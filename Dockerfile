FROM ubuntu:22.04
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git wget lsb-core lsb-release \
    libmysqlclient-dev \
    python3-pip libpq-dev \
    python3-wheel apache2-dev \
    gnuplot

RUN git clone https://github.com/bu-isciii/iskylims.git /srv/iskylims
WORKDIR /srv/iskylims
RUN git checkout develop

RUN pip install -r conf/requirements.txt 

RUN bash install.sh --install app --dev --conf conf/docker_install_settings.txt --docker

WORKDIR /opt/iskylims

# Expose
EXPOSE 8001
# Start the application
CMD ["python3", "/opt/iskylims/manage.py", "runserver", "0:8001"]