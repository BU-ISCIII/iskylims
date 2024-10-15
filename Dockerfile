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
<<<<<<< HEAD
    python3-pip libpq-dev python3-venv python3-wheel \
    apache2-dev cron \
    gnuplot pkg-config rsync
=======
    python3-pip libpq-dev \
    python3-wheel apache2-dev \
    gnuplot pkg-config
>>>>>>> ff357111 (hotfix numpy dependency and docker fixes)

# Set MYSQLCLIENT_CFLAGS and MYSQLCLIENT_LDFLAGS using pkg-config
RUN export MYSQLCLIENT_CFLAGS="$(pkg-config --libs mysqlclient)" && \
    export MYSQLCLIENT_LDFLAGS="$(pkg-config --cflags mysqlclient)"

# Set git repository
RUN mkdir /srv/iskylims 
WORKDIR /srv/iskylims
<<<<<<< HEAD

# Copy the local git repository to docker image directory
COPY . /srv/iskylims

# Create and activate a virtual environment
RUN python3 -m venv /srv/iskylims/venv
ENV PATH="/srv/iskylims/venv/bin:$PATH"
ENV PATH="/usr/sbin/cron:$PATH"
=======
RUN pip install -r conf/requirements.txt 

RUN bash install.sh --install app --conf conf/docker_install_settings.txt --docker
>>>>>>> ff357111 (hotfix numpy dependency and docker fixes)

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