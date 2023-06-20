FROM ubuntu:22.04
# ENV PYTHONUNBUFFERED 1
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git apt-utils  wget gnuplot  python3-pip \
    libmysqlclient-dev vim mysql-client apache2-dev

RUN mkdir /opt/interop
WORKDIR /opt/interop

RUN wget https://github.com/Illumina/interop/releases/download/v1.1.15/InterOp-1.1.15-Linux-GNU.tar.gz
RUN tar -xvf  InterOp-1.1.15-Linux-GNU.tar.gz
RUN ln -s InterOp-1.1.15-Linux-GNU interop
RUN rm InterOp-1.1.15-Linux-GNU.tar.gz


RUN mkdir /opt/iskylims
WORKDIR /opt/iskylims

# RUN git clone https://github.com/BU-ISCIII/iSkyLIMS.git .
RUN git clone https://github.com/smonzon/iSkyLIMS.git .
RUN git checkout develop
# RUN git checkout develop

RUN mkdir -p /opt/iskylims/documents/wetlab/tmp
RUN mkdir -p /opt/iskylims/documents/drylab
RUN mkdir -p /opt/iskylims/logs



# Starting iSkyLIMS
# for develop
RUN python3 -m pip install -r conf/requirements.txt

# for main
# RUN python3 -m pip install -r conf/pythonPackagesRequired.txt
RUN django-admin startproject iskylims .
RUN /bin/bash -c 'grep ^SECRET iskylims/settings.py > ~/.secret'


# Copying config files and script
RUN cp conf/docker_settings.py /opt/iskylims/iskylims/settings.py
RUN cp conf/urls.py /opt/iskylims/iskylims/

RUN sed -i "/^SECRET/c\\$(cat ~/.secret)" iskylims/settings.py
ENV PATH="usr/bin:$PATH"
# Expose and run
EXPOSE 8001
# CMD python3 manage.py runserver 0:8000
