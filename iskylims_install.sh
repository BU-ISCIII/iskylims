#!/usr/bin/bash
#iSkyLIMS version 2.0.0 installation instruction

wget -L https://raw.githubusercontent.com/BU-ISCIII/iSkyLIMS/develop/Dockerfile
wget -L https://raw.githubusercontent.com/BU-ISCIII/iSkyLIMS/develop/docker-compose.yml

docker build -t iskylimsv2 .

docker-compose up -d 
echo "Waiting 30 seconds for starting dabase and web services..."

sleep 30
echo "Creating the database structure for iSkyLIMS"
docker exec -it iskylims_installation_web1_1 python3 manage.py migrate
docker exec -it iskylims_installation_web1_1 python3 manage.py makemigrations django_utils iSkyLIMS_core iSkyLIMS_wetlab iSkyLIMS_drylab iSkyLIMS_clinic
docker exec -it iskylims_installation_web1_1 python3 manage.py migrate

echo "Loading in database initial data"
docker exec -it iskylims_installation_web1_1 python3 manage.py loaddata conf/new_installation_loading_tables.json

echo "Creating super user "
docker exec -it iskylims_installation_web1_1 python3 manage.py createsuperuser
chown -R $USER:$USER .


docker-compose up -d
echo "Deleting dockr files"
rm DockerFile
rm docker-compose.yml

