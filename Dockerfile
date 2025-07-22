FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_snap_workflow:0.15.4-2c75e1-wip-f74bbc

WORKDIR /tmp/docker-build/work/

# Install R, from ArchRProject
RUN apt-get update -y && \
    apt-get install -y \
        r-base \
        r-base-dev \
        apt-transport-https \
        aptitude \
        build-essential \
        gdebi-core \
        gfortran \
        libhdf5-dev \
        libatlas-base-dev \
        libbz2-dev \        
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgdal-dev \
        libgit2-dev \
        libgsl-dev \
        libjpeg-dev \
        libicu-dev \
        liblzma-dev \
        libmagick++-dev \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libpcre3-dev \
        libssl-dev \
        libtcl8.6 \
        libtiff5 \
        libtk8.6 \
        libxml2-dev \
        libxt-dev \
        libx11-dev \
        libtiff-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        locales \
        make \
        pandoc \
        r-cran-rjava \
        tzdata \
        vim \
        wget \
        zlib1g-dev

# Upgrade R to version 4.3.0
RUN wget https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
RUN tar zxvf R-4.3.0.tar.gz
RUN cd R-4.3.0 && ./configure --enable-R-shlib
RUN cd R-4.3.0 && make && make install

# Have to install devtools, cairo like this; see https://stackoverflow.com/questions/20923209
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Install java
RUN apt install -y default-jdk
RUN R CMD javareconf

# Fix systemd conflict with timedatectl
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site

# Installation of R packages with renv
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/renv/renv_1.0.5.tar.gz', repos = NULL, type = 'source')"
COPY atx_snap.Rproj /root/atx_snap.Rproj
COPY renv.lock /root/renv.lock
RUN mkdir /root/renv
COPY renv/activate.R /root/renv/activate.R
WORKDIR /root
RUN R -e "renv::restore()"

# Install BPCells, seurat-disk
RUN R -e "remotes::install_github('bnprks/BPCells/r', ref = 'a3096e5', upgrade = 'never')"
RUN R -e "remotes::install_github('mojaveazure/seurat-disk', ref = '877d4e1', upgrade = 'never')"

RUN R -e "BiocManager::install('sparseMatrixStats')"

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/spam64_2.10-0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "remotes::install_github('jpmcga/ArchR', ref = '9583d51')"

COPY . /root/

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
