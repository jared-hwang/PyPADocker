<h1>PyPA Projects Docker</h1>

<h3>Overview</h3>
    <p>This is the Dockerfile for a PyPE image that contains most modules necessary for operation of packages in the PyPA (PyParticleAccelerator) network. This includes: </p>
        <ul>
            <li>anaconda</li>
            <li>bempp</li>
            <li>warp</li>
            <li>pythonocc</li>
        </ul>
    <p>and a variety of packages that can be seen in the dockerfile, as well as the conda <code>environment.yml</code>.</p>

<h3>Install Docker</h3>
    <ul>
        <li>Visit [the Docker Documentation](https://docs.docker.com/install/) for your given distribution and installation instructions.</li>
    </ul>
    
<h3>Build the PyPA image</h3>
    <ul>
        <li>Run the docker build command (if using linux, <code>sudo</code> must be used with every docker command).</li>
        <br>
        ```<code>docker build -t pypa .</code>```
        <br><br>
        This creates a docker image named <code>pypa</code> from the current directory. The process may take a considerable amount of time, as it is downloading and installing several external modules.
    </ul>

<h3>Run a PyPA image</h3>
    <ul>
        <li>Running the following will open a shell from the previous <code>pypa</code>and remove it once exited (in the interest of reduction of bloat.</li>
            <!-- <br> -->
            <code>docker run -it --rm pyrfq</code><br>
        <li>To create a volume that the container can access on the host computer, add the following to the run command: </li>
            <!-- <br> -->
            <code>-v [absolute path to host directory]:[absolute path to image directory]</code>
            <br>
        <li>To enable the use of GUI applications add the following to the run command (on linux):</li>
            <!-- <br> -->
            <code>
                -v /tmp/.X11-unix:/tmp/.X11-unix \
                -e DISPLAY \
            </code>
            <br>
            Note: You may have to run <code>xhost +local:docker</code> on the host machine to give docker X11 display permissions.
    </ul>


sudo docker run -it --rm \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -e DISPLAY \
    -v /home/rfq-dip-student/jared/docker/transfer/files/:/home/pyrfq/testfiles \
    pyrfq


may have to run xhost +local:docker on host machine to give docker display permissions