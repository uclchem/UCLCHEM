# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'UCLCHEM'
copyright = '2024, Leiden Observatory'
author = "Vermariën"

release = "0.1"
version = "0.1.0"

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.apidoc'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

apidoc_module_dir = '../../src/uclchem'
apidoc_output_dir = 'reference'
apidoc_excluded_paths = ['tests']
apidoc_separate_modules = True



import os
import sys

sys.path.insert(0, os.path.abspath('../../src/uclchem'))
