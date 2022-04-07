#!/usr/bin/bash
#Create markdown file for UCLCHEM docs from the docstrings in the python module
pydoc-markdown -p uclchem --render-toc > uclchemwrap.md