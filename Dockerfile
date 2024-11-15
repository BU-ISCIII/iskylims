FROM ubuntu:24.04
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git wget lsb-release \
    libmysqlclient-dev \
    python3-pip libpq-dev python3-venv python3-wheel \
    apache2-dev \
    gnuplot pkg-config rsync

# Set MYSQLCLIENT_CFLAGS and MYSQLCLIENT_LDFLAGS using pkg-config
RUN export MYSQLCLIENT_CFLAGS="$(pkg-config --libs mysqlclient)" && \
    export MYSQLCLIENT_LDFLAGS="$(pkg-config --cflags mysqlclient)"

# Set git repository
RUN mkdir /srv/iskylims 
WORKDIR /srv/iskylims

# Copy the local git repository to docker image directory
COPY . /srv/iskylims

# Create and activate a virtual environment
RUN python3 -m venv /srv/iskylims/venv

# Install dependencies within the virtual environment
RUN /srv/iskylims/venv/bin/pip install -r conf/requirements.txt

# Set default install type
ARG INSTALL_TYPE=app
ARG GIT_REVISION=main

# Execute the installation script
RUN /bin/bash install.sh --install $INSTALL_TYPE --git_revision $GIT_REVISION --conf conf/docker_install_settings.txt --docker
WORKDIR /opt/iskylims

# Expose
EXPOSE 8001

# Start the application
CMD ["python", "/opt/iskylims/manage.py", "runserver", "0:8001"]