#!/usr/bin/bash

echo "Deploying test containers..."
docker compose build --no-cache 
docker compose up -d

echo "Waiting 20 seconds for starting database and web services..."
sleep 20

echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylims_app python3 manage.py migrate
docker exec -it iskylims_app python3 manage.py makemigrations django_utils core wetlab drylab
docker exec -it iskylims_app python3 manage.py migrate

echo "Loading in database initial data"
docker exec -it iskylims_app python3 manage.py loaddata conf/first_install_tables.json

echo "Creating super user "
docker exec -it iskylims_app python3 manage.py createsuperuser

echo "Download testing files and copy it to samba container"
wget https://zenodo.org/record/8091169/files/iskylims_demo_data.tar.gz
docker cp iskylims_demo_data.tar.gz samba:/mnt
docker exec -it samba tar -xf  /mnt/iskylims_demo_data.tar.gz -C /mnt

echo "deleting compress testing file"
docker exec -it samba rm /mnt/iskylims_demo_data.tar.gz
rm  iskylims_demo_data.tar.gz


echo "Running crontab"
docker exec -it iskylims_app python3 manage.py crontab add
echo "You can now access iskylims via: http://localhost:8001"