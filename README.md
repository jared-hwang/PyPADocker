
# PyPA Projects Docker

### Overview

This is the Dockerfile for a PyPA image that contains most modules necessary for operation of packages in the PyPA (PyParticleAccelerator) network. This includes:

*   anaconda
*   bempp
*   warp
*   pythonocc

and a variety of packages that can be seen in the dockerfile, as well as the conda `environment.yml`.  For more information on Docker in general or to understand the principles behind containerization, see [the Docker documentation](https://docs.docker.com/get-started/) for details. 

### Install Docker

*   Visit [the Docker documentation](https://docs.docker.com/install/) for your given distribution and installation instructions.

### Build the PyPA image

*   Run the docker build command in the dockerfile directory (if using linux, `sudo` must be used with every docker command).

    `docker build -t pypa .`

    This may take a considerable amount of time.

### Run a PyPA image

*   Running the following will open a shell from the previous `pypa`and remove it once exited (in the interest of reduction of bloat).

    `docker run -it --rm pypa`

*   To create a volume that the container can access on the host computer, add the following parameter before pypa:

    `-v [absolute path to host directory]:[desired absolute path in image]`

*   Enable the use of Gui Applications 
    * Linux
        * Add the following parameters before pypa:

        `-v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY`

        <sub>Note: you may have to run <code>xhost +local:docker</code> on the host machine to give docker X11 display permissions.</sub>

    * Windows
        * Download and install [VcXsrv Windows X Server](https://sourceforge.net/projects/vcxsrv/) with all settings allowed (Xming would also work)
        * Get your IP address with `ipconfig` 
        * Add `-e DISPLAY=[IP ADDRESS]:0.0` to the run command
    * Mac
        * Install [XQuartz](https://www.xquartz.org/)
        * Get your IP address with `ipconfig` and run `xhost +[IP ADDRESS]`
        * Add `-e DISPLAY=[IP ADDRESS]:0 -v /tmp/.X11-unix:/tmp/.X11-unix` to the run command
    * For details on any of the above, see [this webpage](https://cuneyt.aliustaoglu.biz/en/running-gui-applications-in-docker-on-windows-linux-mac-hosts)

    
### Customizing the image

The build file is cached, so any changes or added build instructions can be added onto the end without rebuilding the entire image. 

### Note for Windows

Linux containers require a "minimal linux kernal" to be run on Windows, so installation may require a few extra steps. See [this page](https://tutorials.ubuntu.com/tutorial/tutorial-windows-ubuntu-hyperv-containers) for more details.
