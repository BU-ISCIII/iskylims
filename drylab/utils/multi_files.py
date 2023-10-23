# encoding: utf-8
# Generic imports
import json

from django.http import HttpResponse

# Local imports
import drylab.config
import drylab.models


class JSONResponse(HttpResponse):
    """JSONResponse -- Extends HTTPResponse to handle JSON format response.

    This response can be used in any view that should return a json stream of
    data.

    Usage:

        def a_iew(request):
            content = {'key': 'value'}
            return JSONResponse(content, mimetype=response_mimetype(request))

    """

    def __init__(
        self, obj="", json_opts=None, mimetype="application/json", *args, **kwargs
    ):
        json_opts = json_opts if isinstance(json_opts, dict) else {}
        content = json.dumps(obj, **json_opts)
        super(JSONResponse, self).__init__(content, mimetype, *args, **kwargs)


def get_and_save_service_file(request):
    """
    Description:
        The function get the sevice files and stored in database
    Input:
        request      user request with the files
    Return:
        data having the information of the uploaded file
    """
    file_data = {}
    file_data["file"] = request.FILES["file"]
    file_data["file_name"] = file_data["file"].name
    files = [
        {
            "name": file_data["file"].name,
            "type": file_data["file"].content_type,
            "size": file_data["file"].size,
            "deleteType": "DELETE",
        }
    ]
    if request.FILES["file"].size > drylab.config.MAX_UPLOAD_SIZE:
        files[0]["errors"] = "maxFileSize"
        files[0]["error_detail"] = drylab.config.ERROR_FILE_TOO_BIG
    else:
        # store the file

        new_upload_file_obj = (
            drylab.models.UploadServiceFile.objects.create_upload_file(file_data)
        )
        files[0]["file_id"] = new_upload_file_obj.get_upload_file_id()
        files[0]["deleteUrl"] = str(
            "upload-file-delete=" + new_upload_file_obj.get_upload_file_id()
        )

    data = {"files": files}
    return data


def get_uploaded_files(service_obj):
    """
    Description:
        The function return the user uploaded files in a list
    Input:
        service_obj      # service instance
    Return:
        with file and file name stored  on database
    """
    file_list = []
    if drylab.models.UploadServiceFile.objects.filter(
        upload_service=service_obj
    ).exists():
        file_objs = drylab.models.UploadServiceFile.objects.filter(
            upload_service=service_obj
        )
        for file_obj in file_objs:
            file_list.append(
                [
                    file_obj.get_upload_file_full_path_and_name(),
                    file_obj.get_upload_file_name(),
                ]
            )
    return file_list


def update_upload_file_with_service(file_id, service_obj):
    """
    Description:
        The function check if service id exists
    Input:
        service_id      # id of the service
    Return:
        None
    """
    if drylab.models.UploadServiceFile.objects.filter(pk__exact=file_id).exists():
        drylab.models.UploadServiceFile.objects.get(
            pk__exact=file_id
        ).update_service_id(service_obj)
    return


def check_if_file_is_linked_to_service(file_id):
    """
    Description:
        The function check if file is used in a defined service
    Input:
        file_id      # id of the file
    Return:
        True if file is used by a service
    """
    if (
        drylab.models.UploadServiceFile.objects.filter(pk__exact=file_id)
        .exclude(upload_service=None)
        .exists()
    ):
        return True
    return False


def delete_service_file(file_id):
    """
    Description:
        The function delete the file from database
    Input:
        file_id      # id of the file
    Return:
        None
    """
    if drylab.models.UploadServiceFile.objects.filter(pk__exact=file_id).exists():
        drylab.models.UploadServiceFile.objects.filter(pk__exact=file_id).delete()
    return None
