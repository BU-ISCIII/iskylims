#!/usr/bin/bash
#iSkyLIMS version 2.0.0 Docker installation

docker build -t iskylimsv2 .

docker-compose up -d
echo "Waiting 20 seconds for starting dabase and web services..."

sleep 20
echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylimsv2_app python3 manage.py migrate
docker exec -it iskylimsv2_app python3 manage.py makemigrations django_utils core iSkyLIMS_wetlab drylab clinic
docker exec -it iskylimsv2_app python3 manage.py migrate

echo "Loading in database initial data"
docker exec -it iskylimsv2_app python3 manage.py loaddata conf/new_installation_loading_tables.json

echo "Creating super user "
docker exec -it iskylimsv2_app python3 manage.py createsuperuser
echo "Starting iSkyLIMS"
docker exec -it iskylimsv2_app python3 manage.py runserver 0:8000 &

