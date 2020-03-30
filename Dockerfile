# Set the base image (in this case debian)
FROM debian:stable-slim

# Set the working directory
WORKDIR /usr/src/UCLCHEM/

# Copy the current directory to the working directory in the volume
ADD . .

# This hack is widely applied to avoid python printing issues in docker containers.
# See: https://github.com/Docker-Hub-frolvlad/docker-alpine-python3/pull/13
ENV PYTHONUNBUFFERED=1

# Install apk dependencies and packages
RUN apt-get update && apt-get install -y apt-utils python3 python3-pip gcc gfortran

# Alias python3 to python (and pip3 to pip)
RUN ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip

# Install pip (python -m is a workaround to avoid invoking using old wrapper)
# and install packages
RUN python -m pip install --upgrade --force-reinstall pip \
    && pip install numpy \
    && python -mpip install matplotlib