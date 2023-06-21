FROM ubuntu:22.04
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git lsb-core lsb-release

RUN git clone https://github.com/saramonzon/iskylims.git /srv/iskylims
WORKDIR /srv/iskylims
RUN git checkout develop

RUN bash install.sh --install full --dev --conf conf/docker_install_settings.txt --docker

ENV PATH="usr/bin:$PATH"

# Expose
EXPOSE 8001
