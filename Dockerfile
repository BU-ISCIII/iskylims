FROM ubuntu:22.04
ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Updates
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y

# Essential software
RUN apt-get install -y \
    git wget lsb-core lsb-release

RUN git clone https://github.com/saramonzon/iskylims.git /srv/iskylims
WORKDIR /srv/iskylims
RUN git checkout develop

RUN bash install.sh --install full --tables --dev --conf conf/docker_install_settings.txt --docker

WORKDIR /opt/iskylims
ENV PATH="/usr/bin:/opt/iskylims/virtualenv/bin:$PATH"
ENV PYTHONPATH="/opt/iskylims/virtualenv/lib/python3.10/site-packages:${PYTHONPATH}"

# Expose
EXPOSE 8001
# Start the application
CMD ["python3", "/opt/iskylims/manage.py", "runserver", "0:8001"]