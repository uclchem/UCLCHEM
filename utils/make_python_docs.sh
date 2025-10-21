#!/usr/bin/bash
#User must supply the path to their uclchem.github.io repo clone as the first argument
#don't forget to go add the id and title to the notebook markdowns

#Create markdown file for UCLCHEM docs from the docstrings in the python module
api_file=notebooks/start-pythonapi.md
echo "---" > $api_file
echo "id: pythonapi" >> $api_file
echo "title: Python Reference" >> $api_file
echo "---" >> $api_file
pydoc-markdown -p uclchem --render-toc >> $api_file
sed -i 's|# Table of Contents|# Python API|g' $api_file
#get rid of version entries in the markdown
sed -i 's|^.*\_version.*||g' $api_file

#create default parameters file
python utils/generate_param_docs.py src/fortran_src/defaultparameters.f90 notebooks/start-parameters.md

# First run all the notebooks
PYTHON_EXEC=$(which python)
$PYTHON_EXEC utils/run_notebooks.py notebooks

#create markdown file from notebooks then change image links to point to img folder
#of uclchem.github.io. Then move images and notebooks to the right place.
jupyter nbconvert --to markdown notebooks/*.ipynb --output-dir notebooks
#point images to paths in website
sed -i 's|\[png\](.*\/|[png](.\/assets\/|g' notebooks/*.md
#remove styles from tables because they throw mdx errors
find notebooks/*.md -exec python utils/remove_table_styles.py  {} +
find notebooks/*.md -exec sed -i 's/{/\\{/g' {} +
find notebooks/*.md -exec sed -i 's/}/\\}/g' {} +

mv notebooks/*files/*.png $1/website/docs/assets
cp notebooks/assets/* $1/website/docs/assets

mv notebooks/*.md $1/website/docs
rm -r notebooks/*files/
