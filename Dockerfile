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
RUN apt-get update -y && apt-get install -y libz-dev wget git

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

# Install forked version of SnapATAC2
RUN apt-get install -y cmake libclang-dev 
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install git+https://github.com/jpmcga/SnapATAC2.git#subdirectory=snapatac2-python
# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

RUN apt-get update && apt-get install -y libblosc-dev
RUN pip install --no-binary numcodecs numcodecs

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
