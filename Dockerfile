# How to run Xenomapper in a docker container

# Docker File Author / Maintainer
# MAINTAINER Matthew Wakefield <matthew.wakefield@unimelb.edu.au>

# First you will need install docker.  On MacOS you do
# brew install Caskroom/cask/virtualbox
# brew install docker
# brew install boot2docker
# boot2docker download
# boot2docker init
# boot2docker up
# $(boot2docker shellinit)

# You then need to create a docker container.  In the same directory as this file you type

# docker build -t xenomapper .  

# This will run the following executable part of this file

FROM ubuntu:trusty-20190425
MAINTAINER Matthew Wakefield <matthew.wakefield@unimelb.edu.au>
RUN apt-get update && apt-get install -y \
	samtools \
	python3-pip \
	git
RUN pip3 install git+https://github.com/genomematt/xenomapper.git

# congratulations - you now have a docker container with xenomapper installed!

# To test it run
# docker run -i xenomapper xenomapper --version
# (note the first xenomapper is the container name, the second the program name)

# To do anything useful with it you will need to mount some data in the container
# The simplest version of this is to mount a directory from the host.
# This must be an absolute path, and something like $HOME/my/data/directory:/data
# will mount your host directory in the container as /data

# docker run -it -v $HOME/repos/xenomapper/xenomapper/tests/data/:/data/ xenomapper xenomapper --primary_sam /data/paired_end_testdata_human.sam --secondary_sam /data/paired_end_testdata_mouse.sam --unresolve /data/docker_unresolved.sam

# Note that you will get the primary_specific on std_out, the read count category summary on std_err and unresolved in a host file $HOME/repos/xenomapper/xenomapper/tests/data/docker_unresolved.sam
