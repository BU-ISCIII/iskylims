#!/usr/bin/bash

echo "Deploying test containers..."
docker compose build --no-cache 
docker compose up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20

echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylimsv2_app python3 manage.py migrate
docker exec -it iskylimsv2_app python3 manage.py makemigrations django_utils core wetlab drylab
docker exec -it iskylimsv2_app python3 manage.py migrate

echo "Loading in database initial data"
docker exec -it iskylimsv2_app python3 manage.py loaddata conf/first_install_tables.json

echo "Creating super user "
docker exec -it iskylimsv2_app python3 manage.py createsuperuser

echo "You can now access iskylims via: http://localhost:8001"