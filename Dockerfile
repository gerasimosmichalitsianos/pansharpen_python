# use python2 interpreter here. With some minor changes,
# we can use python3.
FROM ubuntu:latest
FROM python:3
COPY . /bin

# import needed python source files
ADD bin/pansharpen.py /

# Update base container install
RUN apt-get update
RUN apt-get upgrade -y

# Install GDAL dependencies
RUN apt-get install -y python3-pip libgdal-dev locales

# Ensure locales configured correctly
RUN locale-gen en_US.UTF-8
ENV LC_ALL='en_US.utf8'

# Set python aliases for python3
RUN echo 'alias python=python3' >> ~/.bashrc
RUN echo 'alias pip=pip3' >> ~/.bashrc

# Update C env vars so compiler can find gdal
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

# This will install latest version of GDAL
RUN pip3 install numpy
RUN pip3 install GDAL==2.2.3
RUN pip3 install PyWavelets

# these are the possible command-line inputs:
# -h','--help
# -d','--directory
# -v','--version
# -o','--output-directory
# -i','--ignore
ENTRYPOINT [ "python3", "pansharpen.py" ]
