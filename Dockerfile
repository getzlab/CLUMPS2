FROM ubuntu:18.04
MAINTAINER Shankara Anand

RUN apt-get update && apt-get install software-properties-common -y &&\
apt-get update &&  apt-get install build-essential git zlib1g-dev \
libgsl-dev libboost-iostreams-dev -y

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jre-headless && \
    apt-get clean

# Install python3
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update
RUN apt-get install -y python3.6 python3-pip python3.6-dev
RUN python3 -m pip install --upgrade pip setuptools

# Copy github
RUN mkdir clumps
COPY . /clumps/

# Install clumps
RUN python3 -m pip install -e ./clumps/.

# Cheese b/c Aaron sucks
RUN python3 -m pip install --no-cache-dir firecloud-dalmatian>=0.0.16

# Test
RUN clumps -h
