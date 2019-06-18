import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eukcc",
    version="0.0.1",
    author="Paul Saary",
    author_email="eukcc@paulsaary.de",
    description="Check eukaryotic genomes or MAGs for completeness and contamination",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    py_modules = ['eukcc', 'base', 'fileoperations'],
    scripts=["scripts/runGMES",
             "scripts/EukCC.py"],
    install_requires=["ete3", "pyfaidx"],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
