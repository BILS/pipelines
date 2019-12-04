FROM nfcore/base:latest

LABEL \
  author="Mahesh Binzer-Panchal, Jacques Dainat" \
  description="Functional annotation preparation environment for NBIS Annotation service" \
  maintainer="mahesh.panchal@nbis.se, jacques.dainat@nbis.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/functional-annotation-preparation/bin:$PATH
