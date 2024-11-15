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
<<<<<<< HEAD
    python3-pip libpq-dev python3-venv python3-wheel \
    apache2-dev cron \
    gnuplot pkg-config rsync
=======
    python3-pip libpq-dev \
    python3-wheel apache2-dev \
    gnuplot pkg-config
>>>>>>> ff357111 (hotfix numpy dependency and docker fixes)
=======
    python3-pip libpq-dev python3-venv python3-wheel \
    apache2-dev \
    gnuplot pkg-config rsync
>>>>>>> 5a7fbaef (update ubuntu version to 24.04)

# Set MYSQLCLIENT_CFLAGS and MYSQLCLIENT_LDFLAGS using pkg-config
RUN export MYSQLCLIENT_CFLAGS="$(pkg-config --libs mysqlclient)" && \
    export MYSQLCLIENT_LDFLAGS="$(pkg-config --cflags mysqlclient)"

# Set git repository
<<<<<<< HEAD
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

<<<<<<< HEAD
RUN bash install.sh --install app --conf conf/docker_install_settings.txt --docker
<<<<<<< HEAD
>>>>>>> ff357111 (hotfix numpy dependency and docker fixes)
=======
=======
RUN bash install.sh --install app --git_revision main --conf conf/docker_install_settings.txt --docker
>>>>>>> 1f1de18d (Solved issue #274. Docker installation fails)
>>>>>>> 8f943a0e (Solved issue #274. Docker installation fails)
=======
RUN mkdir /srv/iskylims
WORKDIR /srv/iskylims
RUN pip install -r conf/requirements.txt

RUN bash install.sh --install app --git_revision main --conf conf/docker_install_settings.txt --docker
>>>>>>> 5a7fbaef (update ubuntu version to 24.04)

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
<<<<<<< HEAD
CMD ["python", "/opt/iskylims/manage.py", "runserver", "0:8001"]
=======
CMD ["python3", "/opt/iskylims/manage.py", "runserver", "0:8001"]
>>>>>>> 5a7fbaef (update ubuntu version to 24.04)
