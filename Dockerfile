# Set the base image (in this case debian)
FROM alpine:3.4

# Set the working directory
WORKDIR /usr/src/UCLCHEM

# Copy the current directory to the working directory in the volume
ADD . .

# This hack is widely applied to avoid python printing issues in docker containers.
# See: https://github.com/Docker-Hub-frolvlad/docker-alpine-python3/pull/13
ENV PYTHONUNBUFFERED=1

# Install apk dependencies and packages
RUN apk update && apk add --no-cache --update-cache openssl python3-dev build-base gcc gfortran freetype-dev libpng-dev \
    && rm -rf /var/lib/apt/lists/* \ 
    && ln -s /usr/include/locale.h /usr/include/xlocale.h

# Alias python3 to python (and pip3 to pip)
RUN ln -s /usr/bin/python3 /usr/bin/python \
    && ln -s /usr/bin/pip3 /usr/bin/pip

# Install packages using pip
RUN pip install --upgrade pip \
    # && pip install -r requirements.txt
    && pip install numpy \
    && python -mpip install matplotlib

# Set the default shell
SHELL [ "bash" ]
