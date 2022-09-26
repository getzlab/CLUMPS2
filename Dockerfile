FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.5
MAINTAINER Shankara Anand

WORKDIR /build

RUN apt-get update && apt-get install software-properties-common -y &&\
apt-get update &&  apt-get install build-essential git zlib1g-dev \
libgsl-dev libboost-iostreams-dev -y

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jre-headless && \
    apt-get clean

COPY . /build

# Install clumps
RUN python3 -m pip install -e .

# Test
RUN clumps -h
