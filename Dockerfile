# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:4.9.2

RUN conda config --set ssl_verify no; conda install -y -c bioconda -c conda-forge \
            'metaeuk=4.a0f584d' pplacer  \
            'epa-ng=0.3.8' python=3.8 hmmer=3.3 minimap2 bwa  pysam biopython samtools=1.12; conda clean --all 
RUN pip install jsonpickle numpy ete3; wget https://github.com/Finn-Lab/EukCC/archive/refs/tags/v.2.0.RC1.zip;\
        apt update; apt install -y unzip; \
        rm -rf /var/lib/apt/lists/*;\
        unzip v.2.0.RC1.zip; \
        cd EukCC-v.2.0.RC1; \
        pip install .; 

        #pip install  "eukcc>=2,<3"

CMD eukcc --help
