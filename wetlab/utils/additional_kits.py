import json

from core.models import CommercialKits, Protocols, ProtocolType, UserLotCommercialKits
from core.utils.protocols import get_protocol_obj_from_id
from wetlab.models import *
from wetlab.utils.library import get_lib_prep_obj_from_id
from wetlab.config import *


def analyze_and_store_input_additional_kits(form_data):
    """
    Description:
        The function get the user input  for the library preparation additional kits
        and store them in database.
    Input:
        form_data   # User form
    Constant:
        HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS
        ERROR_EMPTY_VALUES
    Functions:
        check_empty_fields  # located at utils.library_preparation file
    Return:
        ERROR message if some of the data are missing, or a list of recorded library preparation obj
    """

    lib_prep_ids = form_data["libPrepIds"].split(",")
    lib_prep_code_ids = form_data["libPrepCodeIds"].split(",")
    # additional_kit_ids = form_data['additional_kit_ids'].split(',')
    json_data = json.loads(form_data["used_kits"])
    fixed_heading_length = len(HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS)

    heading_in_excel = form_data["headings"].split(",")
    full_heading_length = len(heading_in_excel)
    stored_additional_kits = {}
    sample_names = []

    """
    if check_empty_fields(json_data) :
        stored_additional_kits['ERROR'] = ERROR_EMPTY_VALUES
        return stored_additional_kits
    """
    protocol_obj = get_lib_prep_obj_from_id(lib_prep_ids[0]).get_protocol_obj()
    for row_index in range(len(json_data)):
        right_id = lib_prep_ids[lib_prep_code_ids.index(json_data[row_index][1])]

        library_prep_obj = get_lib_prep_obj_from_id(right_id)
        sample_names.append(library_prep_obj.get_sample_name())
        for c_index in range(fixed_heading_length, full_heading_length):
            kit_name = heading_in_excel[c_index]
            additional_kit_lib_obj = AdditionaKitsLibPrepare.objects.filter(
                kit_name__exact=kit_name, protocol_id=protocol_obj
            ).last()
            commercial_kit_obj = additional_kit_lib_obj.get_commercial_kit_obj()
            if json_data[row_index][c_index] == "":
                user_lot_commercial_obj = None
            else:
                user_lot_commercial_obj = UserLotCommercialKits.objects.filter(
                    based_commercial=commercial_kit_obj,
                    chip_lot__exact=json_data[row_index][c_index],
                ).last()
                # increase the user Lot Kit use
                user_lot_commercial_obj.set_increase_use()
            user_additional_kit = {}
            user_additional_kit["lib_prep_id"] = library_prep_obj

            user_additional_kit["additionalLotKits"] = additional_kit_lib_obj
            user_additional_kit["userLotKit_id"] = user_lot_commercial_obj

            AdditionalUserLotKit.objects.create_additional_user_lot_kit(
                user_additional_kit
            )

        # update the library prepareation state
        library_prep_obj.set_state("Updated additional kits")

    stored_additional_kits["stored_lib_ids"] = list(
        zip(sample_names, lib_prep_code_ids)
    )

    return stored_additional_kits


def define_table_for_additional_kits(protocol_id):
    """
    Description:
        The function get the commercial protocols defined for the library preparation
        protocol.
    Input:
        protocol_id__exact  # Id of the protocol
    Constant:
        HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    Return:
        kit_data
    """
    kit_data = {}
    protocol_obj = get_protocol_obj_from_id(protocol_id)
    if CommercialKits.objects.filter(protocol_kits=protocol_obj).exists():
        c_kits = CommercialKits.objects.filter(protocol_kits=protocol_obj)
        c_kit_names = []
        for c_kit in c_kits:
            c_kit_names.append(c_kit.get_name())
        kit_data["kits"] = c_kit_names
        kit_data["protocol_name"] = protocol_obj.get_name()
        kit_data["protocol_id"] = protocol_id
        kit_data["heading"] = HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    return kit_data


def get_additional_kits_data_to_modify(protocol_id):
    """
    Description:
        The function get additional kits information for the given protocol
        to allow to modify.
    Input:
        protocol_id    # id of the protocol that are linked to the additional kits
    Constamt:
        HEADING_FOR_MODIFYING_ADDITIONAL_KITS
    Return:
        add_kit_info
    """
    add_kit_info = {}
    protocol_obj = Protocols.objects.get(pk__exact=protocol_id)
    if AdditionaKitsLibPrepare.objects.filter(protocol_id=protocol_obj).exists():
        add_kit_objs = AdditionaKitsLibPrepare.objects.filter(
            protocol_id=protocol_obj
        ).order_by("kit_order")
        add_kit_info["add_kit_data"] = []
        for add_kit_obj in add_kit_objs:
            data = add_kit_obj.get_add_kit_data_for_javascript()
            data.insert(1, "")
            data.append(add_kit_obj.get_add_kit_id())
            add_kit_info["add_kit_data"].append(data)
        if CommercialKits.objects.filter(protocol_kits=protocol_obj).exists():
            c_kits = CommercialKits.objects.filter(protocol_kits=protocol_obj)
            add_kit_info["c_kit_names"] = []
            for c_kit in c_kits:
                add_kit_info["c_kit_names"].append(c_kit.get_name())

        add_kit_info["protocol_name"] = protocol_obj.get_name()
        add_kit_info["protocol_id"] = protocol_id
        add_kit_info["heading"] = HEADING_FOR_MODIFYING_ADDITIONAL_KITS

    return add_kit_info


def get_additional_kits_from_lib_prep(lib_prep_ids):
    """
    Description:
        The function get additional kits (and the user Lot kits) for the  list of
        library preparation ids
    Input:
        lib_prep_ids    # ID List for the library preparation
    Constamt:
        HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS
    Return:
        additional_kits
    """
    additional_kits = {}
    additional_kits["data"] = []
    lib_prep_code_ids = []
    kit_name_list = []
    protocol_obj = LibPrepare.objects.get(pk__exact=lib_prep_ids[0]).get_protocol_obj()
    additional_kits["kit_heading"] = []
    if AdditionaKitsLibPrepare.objects.filter(
        protocol_id=protocol_obj, kitUsed=True
    ).exists():
        additional_kits_objs = AdditionaKitsLibPrepare.objects.filter(
            protocol_id=protocol_obj, kitUsed=True
        ).order_by("kit_order")

        for additional_kit_obj in additional_kits_objs:
            kit_name = additional_kit_obj.get_kit_name()
            kit_name_list.append(kit_name)
            kit_commercial_obj = additional_kit_obj.get_commercial_kit_obj()
            user_lot = []
            if UserLotCommercialKits.objects.filter(
                based_commercial=kit_commercial_obj, run_out=False
            ).exists():
                user_lot_kit_objs = UserLotCommercialKits.objects.filter(
                    based_commercial=kit_commercial_obj, run_out=False
                ).order_by("expiration_date")
                for user_lot_kit_obj in user_lot_kit_objs:
                    user_lot.append(user_lot_kit_obj.get_lot_number())
            additional_kits["kit_heading"].append([kit_name, user_lot])
    additional_kits["fix_heading"] = HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS
    data_length = len(HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS) + len(
        additional_kits_objs
    )

    for lib_prep_id in lib_prep_ids:
        lib_prep_obj = get_lib_prep_obj_from_id(lib_prep_id)
        lib_prep_code_ids.append(lib_prep_obj.get_lib_prep_code())
        row_data = [""] * data_length
        row_data[0] = lib_prep_obj.get_sample_name()
        row_data[1] = lib_prep_obj.get_lib_prep_code()
        additional_kits["data"].append(row_data)
    additional_kits["lib_prep_ids"] = ",".join(lib_prep_ids)
    additional_kits["lib_prep_code_ids"] = ",".join(lib_prep_code_ids)
    additional_kits["full_heading"] = ",".join(
        HEADING_FIX_FOR_ASSING_ADDITIONAL_KITS + kit_name_list
    )

    return additional_kits


def get_additional_kits_list(app_name):
    """
    Description:
        The function get protocol list with information if additional kits are:
        Defined; Pending; or Not requested
    Return:
        data_prot with "protocol_name" "protocol_id" and status of the
        additional index
    """
    additional_kits = []
    if ProtocolType.objects.filter(
        molecule=None, protocol_type__icontains="additional", apps_name__exact=app_name
    ).exists():
        protocol_types = ProtocolType.objects.filter(
            molecule=None, apps_name__exact=app_name
        )
        for protocol_type in protocol_types:
            if Protocols.objects.filter(type=protocol_type).exists():
                protocols = Protocols.objects.filter(type=protocol_type).order_by(
                    "type"
                )
                for protocol in protocols:
                    data_prot = []
                    data_prot.append(protocol.get_name())
                    data_prot.append(protocol.pk)
                    if AdditionaKitsLibPrepare.objects.filter(
                        protocol_id=protocol
                    ).exists():
                        data_prot.append("Defined")
                    else:
                        data_prot.append("")
                    additional_kits.append(data_prot)
    return additional_kits


def get_all_additional_kit_info(protocol_id):
    """
    Description:
        The function get all information of the additional kits for the protocol in
        the input request
    Input:
        protocol_id   # id of the protocol
    Constant:
        HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    Return:
        kit_info
    """
    kit_info = {}
    protocol_obj = get_protocol_obj_from_id(protocol_id)
    if AdditionaKitsLibPrepare.objects.filter(protocol_id=protocol_obj).exists():
        kit_info["kit_heading"] = HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
        kit_info["kit_data"] = []
        kits = AdditionaKitsLibPrepare.objects.filter(
            protocol_id=protocol_obj
        ).order_by("kit_order")
        for kit in kits:
            kit_info["kit_data"].append(kit.get_all_kit_info())

    return kit_info


def get_additional_kit_obj_form_id(add_kit_id):
    """
    Description:
        The function get the information from the additional kits used in the
        library preparation
    Input:
        add_kit_id   # id of the  AdditionaKitsLibraryPreparation
    Return:
        object of the additional kit or None
    """
    if AdditionaKitsLibPrepare.objects.filter(pk__exact=add_kit_id).exists():
        return AdditionaKitsLibPrepare.objects.filter(pk__exact=add_kit_id).last()
    else:
        return None


def get_additional_kits_used_in_sample(sample_id):
    """
    Description:
        The function get the information from the additional kits used in the
        library preparation
    Input:
        sample_id   # sample id
    Constamt:
        HEADING_FOR_DISPLAY_ADDITIONAL_KIT_LIBRARY_PREPARATION
    Return:
        kit_data
    """
    kit_data = {}
    kit_data["protocols_add_kits"] = {}
    if LibPrepare.objects.filter(sample_id__pk__exact=sample_id).exists():
        kit_data[
            "heading_add_kits"
        ] = HEADING_FOR_DISPLAY_ADDITIONAL_KIT_LIBRARY_PREPARATION
        library_preparation_items = LibPrepare.objects.filter(
            sample_id__pk__exact=sample_id
        ).order_by("protocol_id")
        for lib_prep in library_preparation_items:
            protocol_name = lib_prep.get_protocol_used()
            if protocol_name not in kit_data["protocols_add_kits"]:
                kit_data["protocols_add_kits"][protocol_name] = []

            kit_used_objs = AdditionalUserLotKit.objects.filter(lib_prep_id=lib_prep)
            for kit_used_obj in kit_used_objs:
                data = kit_used_obj.get_additional_kit_info()
                data.append(lib_prep.get_lib_prep_code())
                kit_data["protocols_add_kits"][protocol_name].append(data)

    return kit_data


def modify_fields_in_additional_kits(form_data, user):
    """
    Description:
        The function get additional kits used in the library preparation and save
        them in database
    Input:
        form_data   # user form data
        user        # requested user
    Constants:
        HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL
    Return:
        saved_fields
    """
    saved_fields = {}
    protocol_id = form_data["protocol_id"]
    protocol_obj = get_protocol_obj_from_id(protocol_id)
    saved_fields["protocol_name"] = protocol_obj.get_name()
    saved_fields["heading"] = HEADING_ADDING_COMMERCIAL_KITS_TO_PROTOCOL

    saved_fields["fields"] = []
    json_data = json.loads(form_data["add_kit_data"])
    for row_data in json_data:
        if row_data[0] == "" and row_data[1] == "":
            continue
        kit_data = {}

        kit_data["commercialKit_id"] = CommercialKits.objects.filter(
            name__exact=row_data[4]
        ).last()
        kit_data["description"] = row_data[5]
        kit_data["kitOrder"] = row_data[2]
        kit_data["kitUsed"] = row_data[3]
        kit_data["user"] = user

        if row_data[0] == "" and row_data[1] != "":
            # new field
            kit_data["kitName"] = row_data[1]
            kit_data["protocol_id"] = protocol_obj
            saved_fields["fields"].append(
                AdditionaKitsLibPrepare.objects.create_additional_kit(
                    kit_data
                ).get_all_kit_info()
            )
            continue
        if row_data[0] != "" and row_data[1] != "":
            # rename field name
            kit_data["kitName"] = row_data[1]
        else:
            kit_data["kitName"] = row_data[0]
        # Update  Field
        add_kit_obj = get_additional_kit_obj_form_id(row_data[-1])
        if not add_kit_obj:
            # Unable to find the object class. Skipping this change
            continue
        add_kit_obj.update_add_kit_fields(kit_data)
        saved_fields["fields"].append(add_kit_obj.get_all_kit_info())
    return saved_fields


def set_additional_kits(form_data, user):
    """
    Description:
        The function get additional kits used in the library preparation and save
        them in database
    Input:
        form_data   # user form data
        user        # requested user
    Return:
        kit_names
    """
    kit_names = []
    protocol_id = form_data["protocol_id"]
    kit_heading_names = [
        "kitName",
        "kitOrder",
        "kitUsed",
        "commercial_kit",
        "description",
    ]
    json_data = json.loads(form_data["kits_data"])
    for row_index in range(len(json_data)):
        if json_data[row_index][0] == "":
            continue
        kit_data = {}
        for i in range(len(kit_heading_names)):
            kit_data[kit_heading_names[i]] = json_data[row_index][i]
        kit_data["protocol_id"] = get_protocol_obj_from_id(protocol_id)
        kit_data["commercialKit_id"] = CommercialKits.objects.filter(
            name__exact=kit_data["commercial_kit"]
        ).last()
        kit_data["user"] = user
        AdditionaKitsLibPrepare.objects.create_additional_kit(kit_data)
        kit_names.append(kit_data["kitName"])
    return kit_names
