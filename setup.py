import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = {}
with open("eukcc/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    name="eukcc",
    version=version["__version__"],
    author="Paul Saary",
    author_email="eukcc@paulsaary.de",
    description="Check eukaryotic genomes or MAGs for completeness and contamination",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Finn-Lab/EukCC/",
    py_modules=["workflow", "base", "fileoperations"],
    scripts=["scripts/runGMES", "scripts/filter_euk_bins.py"],
    entry_points={"console_scripts": ["eukcc = eukcc.__main__:main"]},
    install_requires=["ete3", "pyfaidx", "configargparse", "pyyaml", "pygmes"],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Unix",
    ],
)
