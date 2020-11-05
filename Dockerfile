FROM bioconductor/bioconductor_docker:RELEASE_3_11

#Install R packages
RUN R -e "install.packages('ordinal')"
RUN R -e "install.packages('dplyr')"
RUN R -e "install.packages('forcats')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('metafolio')"
RUN R -e "install.packages('patchwork')"
RUN R -e "BiocManager::install('biomaRt', update=T, ask=F)"

RUN apt-get update -y && apt-get install -y \
    python

# Install scripts
RUN mkdir /home/analyst
COPY *.sh /home/analyst/
COPY *.py /home/analyst/
COPY *.R /home/analyst/
RUN chmod a+x /home/analyst/*.sh

