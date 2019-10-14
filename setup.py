import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eukcc_dev",
    version="0.0.2",
    author="Paul Saary",
    author_email="eukcc@paulsaary.de",
    description="Check eukaryotic genomes or MAGs for completeness and contamination",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    py_modules=['workflow', 'base', 'fileoperations'],
    scripts=["scripts/runGMES"],
    entry_points={
        'console_scripts': ['eukcc = eukcc.__main__:main']
        },
    install_requires=["ete3", "pyfaidx", "configargparse",
                      "PyQt5"],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
