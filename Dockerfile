#brew install Caskroom/cask/virtualbox
#brew install docker
#brew install boot2docker
#boot2docker download
#boot2docker up
#$(boot2docker shellinit)




FROM ubuntu:14.04
RUN apt-get update
RUN apt-get install -y python3-pip
RUN apt-get install -y git
RUN pip3 install git+https://github.com/genomematt/xenomapper.git



# docker build -t xenomapper_container .  

# docker run -i xenomapper_container