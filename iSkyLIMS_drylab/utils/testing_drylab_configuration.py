# Generic imports
import os

import django.contrib.auth.models
from django.conf import settings

# Local imports
import iSkyLIMS_drylab.models


def get_config_file(config_file):
    c_file = []
    try:
        with open(config_file, "r") as fh:
            for line in fh:
                if "PASSWORD" in line:
                    hide_passwd = line.split("=")
                    hide_passwd[1] = "XXXXXXXXXXXXXXXXX"
                    line = " = ".join(hide_passwd)
                line = line.replace("\n", "")
                c_file.append(line)
    except Exception:
        return
    return c_file


def get_iSkyLIMS_settings():
    s_file = []
    settings_file = os.path.join(settings.BASE_DIR, "iSkyLIMS", "settings.py")
    try:
        with open(settings_file, "r") as fh:
            for line in fh:
                if "PASSWORD" in line or "SECRET_KEY" in line:
                    if "AUTH_PASSWORD_VALIDATORS" not in line:
                        if "=" in line:
                            split_separator = "="
                        else:
                            split_separator = ":"
                        hide_passwd = line.split(split_separator)
                        hide_passwd[1] = "XXXXXXXXXXXXXXXXX"
                        line = " = ".join(hide_passwd)

                line = line.replace("\n", "")
                s_file.append(line)
    except Exception:
        return

    return s_file


def create_service_test(service_requested):
    service_results = []
    # Check if test service exists
    if iSkyLIMS_drylab.models.Service.objects.filter(
        service_request_number__exact=service_requested
    ).exists():
        delete_service = iSkyLIMS_drylab.models.Service.objects.get(
            service_request_number__exact=service_requested
        )
        delete_service.delete()

    # Check user is defined in database
    if not django.contrib.auth.models.User.objects.filter(
        username__exact="test_userDrylab"
    ).exists():
        django.contrib.auth.models.User.objects.create_user(
            username="test_userDrylab",
            email="test_userDrylab@iSkyLIMS.com",
            password="test_userD",
        )

    try:
        user_name = django.contrib.auth.models.User.objects.get(
            username__exact="test_userDrylab"
        )
        service_results.append(("User defined", "OK"))
    except Exception:
        service_results.append(("User defined", "NOK"))
    try:
        service_platform = iSkyLIMS_drylab.models.Platform.objects.first()
        service_results.append(("Platform defined", "OK"))
    except Exception:
        service_results.append(("Platform defined", "NOK"))
    try:
        service_file_ext = iSkyLIMS_drylab.models.FileExt.objects.first()
        service_results.append(("File extension defined", "OK"))
    except Exception:
        service_results.append(("File extension defined", "NOK"))

    for i in range(len(service_results)):
        if "NOK" in service_results[i]:
            return service_results, "NOK"
    try:
        new_test_service = iSkyLIMS_drylab.models.Service(
            service_request_number=service_requested,
            serviceUserId=user_name,
            servicePlatform=service_platform,
            serviceFileExt=service_file_ext,
            serviceStatus="recorded",
        )
        new_test_service.save()

        service_results.append(("Service Test creation", "OK"))
        return service_results, "OK"
    except Exception:
        service_results.append(("Service Test creation", "NOK"))
        return service_results, "NOK"


def create_resolution_test(resolution_number, service_requested):
    resolution_test = []
    # get service object
    service = iSkyLIMS_drylab.models.Service.objects.get(
        service_request_number=service_requested
    )

    # Create resolution object
    try:
        test_resolution = iSkyLIMS_drylab.models.Resolution(
            resolution_serviceID=service,
            resolutionNumber=resolution_number,
            resolution_full_number=str("Test_" + resolution_number),
        )
        resolution_test.append(("Resolution creation", "OK"))
    except Exception:
        resolution_test.append(("Resolution creation", "NOK"))

    test_resolution.get_resolution_information()

    resolution_test.append(("Folder structure creation", "OK"))

    return resolution_test


def delete_test_service(service_name):
    try:
        d_service = iSkyLIMS_drylab.models.Service.objects.get(
            service_request_number__exact=service_name
        )
        d_service.delete()
    except Exception:
        return

    return
