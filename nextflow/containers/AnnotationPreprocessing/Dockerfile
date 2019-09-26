FROM nfcore/base:latest

LABEL \
  author="Mahesh Binzer-Panchal" \
  description="Annotation preprocessing environment for NBIS Annotation service" \
  maintainer="mahesh.panchal@nbis.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/annotation-preprocessing/bin:$PATH

