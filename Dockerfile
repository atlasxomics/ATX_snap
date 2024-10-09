FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

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
RUN apt-get install libz-dev wget

# Latch SDK
# DO NOT REMOVE
RUN pip install latch==2.52.2
RUN mkdir /opt/latch

# Install pip dependencies from `requirements.txt`
RUN pip install numpy==1.25.2

COPY requirements.txt /opt/latch/requirements.txt
RUN pip install --requirement /opt/latch/requirements.txt
RUN pip install 'rapids-singlecell[rapids12]' --extra-index-url=https://pypi.nvidia.com
RUN wget https://developer.download.nvidia.com/compute/cuda/12.6.2/local_installers/cuda-repo-debian11-12-6-local_12.6.2-560.35.03-1_amd64.deb
RUN dpkg -i cuda-repo-debian11-12-6-local_12.6.2-560.35.03-1_amd64.deb
RUN cp /var/cuda-repo-debian11-12-6-local/cuda-*-keyring.gpg /usr/share/keyrings/
RUN add-apt-repository contrib
RUN apt-get update
RUN apt-get -y install cuda-toolkit-12-6
RUN apt-get install -y nvidia-open
RUN apt-get install -y cuda-drivers
# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
