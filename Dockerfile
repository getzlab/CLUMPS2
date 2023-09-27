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

# Install prerequisites 
RUN pip install twobitreader statsmodels scipy pyopenssl prody mkl-random mkl-fft lxml jpype1 canine biopython tqdm

# Install clumps
COPY . /build
RUN python3 -m pip install -e .

#add ipdb for easier debugging
RUN pip install ipdb
RUN export PYTHONBREAKPOINT=ipdb.set_trace

# Test
RUN clumps -h
