#!/usr/bin/env python3
import random
import itertools
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path(__file__).parent.parent.resolve()
DIST_DIR = PROJECT_ROOT / "dist"
LOG_BRIEF = PROJECT_ROOT / "logs" / "dependency_matrix.log"
LOG_FULL = PROJECT_ROOT / "logs" / "dependency_matrix_full.log"

# Dependency matrices to test
DEPENDENCIES = {
    "numpy": ["1.21.*", "2.0.*","2.3.*"][::-1],
    "networkx": ["2.5.*", "3.0.*", "3.5.*"][::-1],
    "tqdm": ["4.44.*", "4.60.*","4.67.*"][::-1],
    "matplotlib": ["3.3.*", "3.7.*", "3.10.*"][::-1],
}
for v in DEPENDENCIES.values():
    random.shuffle(v)
OPTIONAL_DEPENDENCIES = {
    "rasterio": ["1.3.*", "1.4.*"][::-1],
    "xarray": ["2023.*", "2024.*", "2025.*"][::-1],
}
for v in OPTIONAL_DEPENDENCIES.values():
    random.shuffle(v)

# Python executables to try (filter to ones that exist)
PYTHONS = [
    "/Users/alex/mambaforge/envs/p3.9/bin/python",
    "/Users/alex/mambaforge/envs/p3.10/bin/python",
    "/Users/alex/mambaforge/envs/p3.11/bin/python",
    "/Users/alex/mambaforge/envs/p3.12/bin/python",
    "/Users/alex/mambaforge/envs/p3.13/bin/python",
]
random.shuffle(PYTHONS)

def log(write_to, msg):
    write_to.write(msg + "\n")
    write_to.flush()
    print(msg)

def run(cmd, env=None, cwd=None, log_file=None):
    if log_file:
        log(log_file, f"[CMD] {' '.join(cmd)}")
    subprocess.check_call(cmd, env=env, cwd=cwd, stdout=log_file, stderr=subprocess.STDOUT)

def py_tag(python_exe):
    # e.g., cp310, cp313
    out = subprocess.check_output([python_exe, "-c", "import sys;print(f'cp{sys.version_info.major}{sys.version_info.minor}')"], text=True)
    return out.strip()

def build_wheel(python_exe, longlog):
    # Clean build artifacts
    for p in (PROJECT_ROOT / "build", PROJECT_ROOT / "PyOCN.egg-info"):
        shutil.rmtree(p, ignore_errors=True)
    DIST_DIR.mkdir(exist_ok=True)

    # Build a wheel for this Python
    run([python_exe, "-m", "pip", "install", "-U", "pip", "build", "wheel", "setuptools"], log_file=longlog)
    run([python_exe, "-m", "build", "--wheel", "--outdir", str(DIST_DIR)], log_file=longlog)

    # Pick the wheel matching the Python tag
    tag = py_tag(python_exe)
    wheels = sorted(DIST_DIR.glob(f"pyocn-*{tag}*.whl"))
    if not wheels:
        raise RuntimeError(f"No wheel found for {tag} in {DIST_DIR}")
    return wheels[-1]

def combo_iterator():
    all_deps = {**DEPENDENCIES, **OPTIONAL_DEPENDENCIES}
    keys, vals = zip(*all_deps.items())
    for vals_combo in itertools.product(*vals):
        yield dict(zip(keys, vals_combo))

def test_one_combo(python_exe, deps, wheel_path, shortlog, longlog):
    if not Path(python_exe).exists():
        log(longlog, f"[SKIP] Missing Python: {python_exe}")
        return False

    env_dir = Path(tempfile.mkdtemp(prefix="pyocn-test-"))
    work_dir = Path(tempfile.mkdtemp(prefix="pyocn-work-"))
    try:
        # Create venv
        run([python_exe, "-m", "venv", str(env_dir)], log_file=longlog)
        pip = env_dir / ("Scripts" if sys.platform.startswith("win") else "bin") / "pip"
        py = env_dir / ("Scripts" if sys.platform.startswith("win") else "bin") / "python"

        # Upgrade build tooling
        run([str(pip), "install", "-U", "pip", "wheel", "setuptools"], log_file=longlog)

        # Install pinned dependencies (may fail due to incompatibilities)
        pinned = [f"{pkg}=={ver}" for pkg, ver in deps.items()]
        run([str(pip), "install"] + pinned, log_file=longlog)

        # Install our wheel
        run([str(pip), "install", str(wheel_path)], log_file=longlog)

        # Headless test code (no GUI needed)
        test_code = r"""
import os, tempfile, pathlib, sys
os.environ["MPLBACKEND"] = "Agg"
# ensure repo root is not in sys.path
sys.path = [p for p in sys.path if not (p and p == os.getcwd())]
import PyOCN as po
po.SUPPRESS_WARNINGS = True
ocn = po.OCN.from_net_type("I", (10, 10), gamma=0.5, random_state=2126666666, verbosity=0)
res = ocn.fit(n_iterations=1_000, xarray_out=True)
dag = ocn.to_digraph()
ocn2 = po.OCN.from_digraph(dag, gamma=0.5)
po.plot_ocn_raster(ocn2)  # smoke only
out = str(pathlib.Path(tempfile.gettempdir()) / "pyocn_test_output.tif")
ocn.to_gtiff(0, 0, "ESRI:102336", out)
print("OK", out)
"""
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        run([str(py), "-c", test_code], cwd=str(work_dir), env=env, log_file=longlog)

        log(shortlog, f"[PASS] py={python_exe} deps={deps}")
        log(longlog, f"[PASS] py={python_exe} deps={deps}")
        return True

    except subprocess.CalledProcessError:
        log(shortlog, f"[FAIL] py={python_exe} deps={deps}")
        log(longlog, f"[FAIL] py={python_exe} deps={deps}")
        return False
    finally:
        shutil.rmtree(env_dir, ignore_errors=True)
        shutil.rmtree(work_dir, ignore_errors=True)

def main():
    with open(LOG_BRIEF, "w") as shortlog, open(LOG_FULL, "w") as longlog:
        log(longlog, f"== PyOCN dependency matrix test == {datetime.now().isoformat()} ==")
        existing_pythons = [p for p in PYTHONS if Path(p).exists()]
        missing_pythons = [p for p in PYTHONS if not Path(p).exists()]
        if missing_pythons:
            log(longlog, f"[INFO] Missing Python executables (skipped): {missing_pythons}")
        if not existing_pythons:
            log(longlog, "[ERROR] No Python executables found. Aborting.")
            sys.exit(1)

        for py_exe in existing_pythons:
            log(longlog, f"\n=== Building wheel for {py_exe} ===")
            try:
                wheel = build_wheel(py_exe, longlog)
                log(longlog, f"[INFO] Built wheel: {wheel.name}")
            except Exception as e:
                log(longlog, f"[ERROR] Build failed for {py_exe}: {e}")
                continue

            # Iterate all combinations (core + optional)
            for deps in combo_iterator():
                test_one_combo(py_exe, deps, wheel, shortlog, longlog)

        log(longlog, "== Done ==")

if __name__ == "__main__":
    os.chdir(PROJECT_ROOT)
    main()