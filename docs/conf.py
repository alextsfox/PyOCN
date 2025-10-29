# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Tell Sphinx where your package is
import sys
import os
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyOCN'
copyright = '2025, Alexander S. Fox'
author = 'Alexander S. Fox'
version = "1.4.20251029"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# autodoc_mock_imports = ["matplotlib", "networkx", "numpy", "xarray", "rasterio", "tqdm"]

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode', 
    'sphinx.ext.napoleon',
    # 'sphinx_autodoc_typehints',  # pip install sphinx-autodoc-typehints
    'myst_parser',              # pip install myst-parser
]
# MyST settings
myst_enable_extensions = ["colon_fence", "deflist"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme' # 'alabaster' # pip install sphinx-rtd-theme
html_static_path = ['_static']
