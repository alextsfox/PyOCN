if __name__ == "__main__":
    import ast
    from collections import defaultdict
    import pandas as pd
    logfile = "/Users/alex/Documents/pet-projects/lang-learn/C/OCN/logs/dependency_matrix.log"
    with open(logfile, "r") as f:
        lines = [ln.strip() for ln in f.readlines()]
        
    results = defaultdict(list)
    for i, ln in enumerate(lines):
        success = ln[1:5] == "PASS"
        python_ver = ln.split("=")[1].strip().split("/")[-3][1:]
        
        deps = ast.literal_eval(ln.split("=")[2])
        depnames = list(deps.keys())
        depvers = list(deps.values())

        results["success"].append(success)
        results["python_ver"].append(python_ver)
        for n, v in deps.items():
            results[n].append(float(v.strip(".*")))
    successes = pd.DataFrame(results).query("success").drop(columns=["success"])
    # find the pareto-minimal combos of dependency versions for each python version
    minimal_combos = []
    for py_ver, group in successes.groupby("python_ver"):
        group = group.reset_index(drop=True)
        dominated = set()
        for i, row_i in group.iterrows():
            for j, row_j in group.iterrows():
                if i == j or j in dominated:
                    continue
                le_all = True
                lt_any = False
                for dep in depnames:
                    vi = row_i[dep]
                    vj = row_j[dep]
                    if vi > vj:
                        le_all = False
                        break
                    if vi < vj:
                        lt_any = True
                if le_all and lt_any:
                    dominated.add(i)
                    break
        minimal_combos.extend(group.drop(index=dominated).to_dict(orient="records"))
    pd.DataFrame(minimal_combos).to_csv("/Users/alex/Documents/pet-projects/lang-learn/C/OCN/logs/dep_test_results.csv", index=False)