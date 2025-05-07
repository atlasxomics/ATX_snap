FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base-cuda12:565f-main

WORKDIR /tmp/docker-build/work/

SHELL [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
ENV TZ='Etc/UTC'
ENV LANG='en_US.UTF-8'

ARG DEBIAN_FRONTEND=noninteractive

# Install base packages
RUN apt-get update -y && apt-get install -y libz-dev git

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

# Latch SDK
# DO NOT REMOVE
RUN pip install latch==2.52.2
RUN mkdir /opt/latch

# Install specific version of numpy
RUN pip install numpy==1.25.2

# Install pip dependencies from `requirements.txt`, pychromvar with chunks
COPY requirements.txt /opt/latch/requirements.txt
RUN pip install --requirement /opt/latch/requirements.txt

RUN pip install 'rapids-singlecell[rapids12]' --extra-index-url=https://pypi.nvidia.com
RUN pip install --no-cache-dir git+https://github.com/pinellolab/pychromVAR.git@7fc47cb02ed36e0ce4c53c5c08bfe17b1ee626a7
RUN git clone https://github.com/latchbio-workflows/gmacs.git && \
    cd gmacs && \
    pip install -e .

RUN pip3 uninstall -y aiobotocore botocore awscli s3transfer
RUN pip3 install awscli

RUN apt-get update -y && apt-get install -y bedtools
ENV TQDM_DISABLE=1

# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

# Fix bug in snapatac2.tools.diff_test()
COPY _diff.py /usr/local/lib/python3.11/dist-packages/snapatac2/tools/_diff.py 

# Copy stuff for ArchRProject
COPY atx_snap.Rproj /root/atx_snap.Rproj
COPY .renvignore /root/.renvignore

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
