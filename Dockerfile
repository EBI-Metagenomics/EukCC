# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:4.9.2

RUN conda config --set ssl_verify no; conda install -y -c bioconda -c conda-forge \
            'metaeuk=4.a0f584d' pplacer  \
            'epa-ng=0.3.8' python=3.8 hmmer=3.3 git \
            minimap2 bwa  pysam biopython \
            samtools=1.12; conda clean --all ; \
            pip install jsonpickle numpy ete3;
ADD https://api.github.com/repos/Finn-Lab/EukCC/git/refs/heads/master version.json
RUN pip install https://github.com/Finn-Lab/EukCC/archive/refs/heads/master.zip

ENTRYPOINT ["eukcc"]
