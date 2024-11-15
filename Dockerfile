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
RUN git checkout develop

RUN pip install -r conf/requirements.txt 

RUN bash install.sh --install app --git_revision main --conf conf/docker_install_settings.txt --docker

WORKDIR /opt/iskylims

# Expose
EXPOSE 8001
# Start the application
CMD ["python3", "/opt/iskylims/manage.py", "runserver", "0:8001"]