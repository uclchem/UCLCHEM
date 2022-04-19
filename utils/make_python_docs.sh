#!/usr/bin/bash
#User must supply the path to their uclchem.github.io repo clone as the first argument
#don't forget to go add the id and title to the notebook markdowns

#Create markdown file for UCLCHEM docs from the docstrings in the python module
api_file=Tutorials/pythonapi.md
echo "---" > $api_file
echo "id: pythonapi" >> $api_file
echo "title: Python Reference" >> $api_file
echo "---" >> $api_file
pydoc-markdown -p uclchem --render-toc >> $api_file


#create default parameters file
python utils/generate_param_docs.py src/fortran_src/defaultparameters.f90 Tutorials/parameters.md

#create markdown file from notebooks then change image links to point to img folder
#of uclchem.github.io. Then move images and notebooks to the right place.
jupyter nbconvert --to markdown Tutorials/*.ipynb --output-dir Tutorials
sed -i 's|\[png\](.*\/|[png](img\/|g' Tutorials/*.md
mv Tutorials/*files/*.png $1/website/static/img
mv Tutorials/*.md $1/website/docs
