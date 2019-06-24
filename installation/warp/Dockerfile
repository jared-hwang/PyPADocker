FROM ubuntu:18.10

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    sudo \
    wget \
    make \
    git \
    vim \
    gcc \
    gfortran \
    libx11-dev \
    openmpi-bin libopenmpi-dev \
    python3 \
    python3-pip \
    python3-numpy \
    python3-scipy \
    python3-mpi4py \
    python3-h5py \
    && rm -rf /var/lib/apt/lists/*

# This is a bit of a hack so that the name "python" is defined.
RUN ln -sf /usr/bin/python3 /usr/bin/python

# openPMD-viewer is installed mainly for tests
# Note: matplotlib is installed with pip since the apt-get install matplotlib
#       needs the time zone to be set.
RUN pip3 --no-cache-dir install matplotlib \
    openPMD-viewer \
    Forthon

# Install pygist
ENV GIT_SSL_NO_VERIFY 1
RUN git clone https://bitbucket.org/dpgrote/pygist.git \
    && cd pygist \
    && python3 setup.py config \
    && python3 setup.py install \
    && cd ../ \
    && rm -rf pygist

# Create a new user and copy the current branch of Warp
# into the Docker container
RUN useradd --create-home warp_user
RUN mkdir /home/warp_user/warp/
COPY ./ /home/warp_user/warp/

# Compile warp
RUN cd /home/warp_user/warp/pywarp90 \
    && make cleanall \
    && rm -f *local* \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3 \
    && echo 'FCOMP= -F gfortran' >> Makefile.local3.pympi \
    && echo 'FCOMPEXEC= --fcompexec mpifort' >> Makefile.local3.pympi \
    && make install3 \
    && make clean3 \
    && make pinstall3 \
    && make pclean3

RUN chown -R warp_user /home/warp_user/warp/
RUN chgrp -R warp_user /home/warp_user/warp/

# Grant sudo access without password
RUN echo 'warp_user ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

# Switch to the new user
WORKDIR /home/warp_user
USER warp_user

# This is needed to get around a bug in openmpi that would print copious error messages
# Unfortunately, this turns off CMA and uses shared memory for communication.
# An alternative is to do "docker run --cap-add SYS_PTRACE ...", which keeps CMA.
ENV OMPI_MCA_btl_vader_single_copy_mechanism none

# Prepare the run directory
RUN mkdir run/
WORKDIR /home/warp_user/run/
