#
# Dockerfile for PyRFQ project
# Jared Hwang, June 2019
#

FROM ubuntu:18.04

# Updating ubuntu packages
RUN apt-get update && yes|apt-get upgrade

RUN apt-get install -y \
    emacs \
    sudo \
    wget \
    bzip2 \
    make \
    git \
    vim \
    gcc \
    gfortran \
    libx11-dev \
    openmpi-bin libopenmpi-dev \
    python3-mpi4py \
    python3-h5py \
    python3-pip \
    python-pip \
    && rm -rf /var/lib/apt/lists*

RUN pip install Forthon \
    && pip3 install Forthon

# Establishing home environment
# RUN useradd --create-home pyrfq
RUN adduser --disabled-password --gecos '' pyrfq
RUN adduser pyrfq sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
WORKDIR /home/pyrfq
RUN mkdir /home/pyrfq/installation/
COPY ./transfer/installation /home/pyrfq/installation

# Installing anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
RUN bash Anaconda3-2019.03-Linux-x86_64.sh -b
RUN rm Anaconda3-2019.03-Linux-x86_64.sh

# Set path to conda
# ENV PATH /home/pyrfq/anaconda3/bin:$PATH
ENV PATH /root/anaconda3/bin:$PATH

# Updating conda packages
RUN conda update conda

# Create (PyRFQ) environment from yml and activate it automatically
RUN conda env create -f ./installation/environment.yml
RUN echo "source activate PyRFQ" > ~/.bashrc

ENV PATH /root/anaconda3/envs/PyRFQ/bin:$PATH
RUN pip install Forthon 
# Install pygist
RUN cd ./installation/pygist \
    && python3 setup.py config \
    && python3 setup.py install \
    && cd ../

RUN pip install Forthon \
    && pip3 install Forthon

# Install warp as root
RUN cd ./installation/warp/pywarp90 \
    && make cleanall \
    && rm -f *local* \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3 \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3.pympi \
    && echo 'FCOMPEXEC= --fcompexec mpifort' >> Makefile.local3.pympi \
    && make install3 \
    && make clean3 \    
    && make pinstall3 \
    && make pclean3

RUN apt-get update \
    && apt-get install -y \
       dh-make \
       gcc \
       g++ \
       gfortran \
       libeigen3-dev \
       python-dev python3-dev \
       python-numpy python3-numpy \
       patchelf \
       libtbb-dev \
       zlib1g-dev \
       libboost-all-dev \
       cmake \
       cpio \
       libdune-common-dev libdune-geometry-dev libdune-grid-dev libdune-localfunctions-dev \
       mpi-default-dev \
       cython cython3 \
       python-setuptools python3-setuptools \
       python-mpi4py python3-mpi4py \
       python-scipy python3-scipy \
       dh-python 
       python-all python3-all \
       libpython3.5 libpython2.7 \
    && conda install -c anaconda h5py 

# ENV PATH /root/anaconda3/bin:$PATH
RUN cd ./installation/bempp \
    && python setup.py install

RUN cd ./installation/dans_pymodules \
    && pip install .

RUN cd ./installation/pythonocc-utils \
    && pip install .
