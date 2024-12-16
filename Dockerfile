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
RUN pip3 uninstall -y awscli boto3 botocore s3transfer
RUN pip3 install awscli

# Install specific version of numpy
RUN pip install numpy==1.25.2

COPY requirements.txt /opt/latch/requirements.txt
RUN pip install --requirement /opt/latch/requirements.txt
RUN pip install 'rapids-singlecell[rapids12]' --extra-index-url=https://pypi.nvidia.com
RUN pip install zstandard
RUN pip install --no-cache-dir git+https://github.com/pinellolab/pychromVAR.git@7fc47cb02ed36e0ce4c53c5c08bfe17b1ee626a7

# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
