#
# Dockerfile for PyPA related projects
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
    gfortran \
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
    dh-python \
    python-all python3-all \
    libpython3.5 libpython2.7 \
    gmsh \
    && rm -rf /var/lib/apt/lists*

RUN pip install Forthon \
    && pip3 install Forthon

# Establishing home environment
RUN adduser --disabled-password --gecos '' pypa
RUN adduser pypa sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
WORKDIR /home/pypa
RUN mkdir /home/pypa/installation/
COPY ./installation /home/pypa/installation

# Installing anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
RUN bash Anaconda3-2019.03-Linux-x86_64.sh -b
RUN rm Anaconda3-2019.03-Linux-x86_64.sh

# Set path to conda
ENV PATH /root/anaconda3/bin:$PATH

# Create (PyPA) environment from yml and activate it automatically
RUN conda update conda
RUN conda env create -f ./installation/environment.yml
RUN echo "source activate PyPA" > ~/.bashrc
ENV PATH /root/anaconda3/envs/PyPA/bin:$PATH

# Install serial and parallel Warp and related modules
RUN cd ./installation/pygist \
    && python3 setup.py config \
    && python3 setup.py install \
    && cd ../

RUN pip install Forthon \
    && pip3 install Forthon

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

# Install bempp and related modules
RUN conda install -c anaconda h5py 
RUN cd ./installation/bempp \
    && python setup.py install

# dans_pymodules
RUN cd ./installation/dans_pymodules \
    && pip install .

# pythonocc-utils required with pythonocc
RUN cd ./installation/pythonocc-utils \
    && pip install .