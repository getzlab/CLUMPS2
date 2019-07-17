FROM bitnami/minideb
MAINTAINER Shankara Anand

RUN install_packages build-essential

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Install python3
RUN install_packages python3 python3-pip
