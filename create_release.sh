
# compute coverage
coverage run -m pytest
coverage report -m
rm badges/coverage.svg
coverage-badge -o badges/coverage.svg

# make sure tests are all in order
pytest

# now build package
python3 setup.py sdist bdist_wheel
