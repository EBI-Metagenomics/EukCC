#jupyter nbconvert --to script treelineage/treelineage.ipynb
python3 setup.py sdist bdist_wheel
pip3 install -e .
