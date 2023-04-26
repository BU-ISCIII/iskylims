from iSkyLIMS_drylab.models import Service, ServiceState


"""
    The script is applicable for the upgrade from 2.3.0 to 2.3.1.
    Service state that was defined as option choice in  models,is replaced in
    version 2.3.1  as a primary key in SampleState table, given in this way
    more flexibility for adding more states in the future.

    Script reads the existing service state value reorded in field serviceStatus
    and create the relation to the same state in SampleState table.
"""


def run():
    service_states = [
        "recorded",
        "approved",
        "rejected",
        "queued",
        "in_progress",
        "delivered",
        "archived",
    ]
    s_state_obj = {}
    for s_state in service_states:
        if ServiceState.objects.filter(state_value__iexact=s_state).exists():
            s_state_obj[s_state] = ServiceState.objects.filter(
                state_value__iexact=s_state
            ).last()
    for s_obj in Service.objects.all():
        try:
            s_obj.update_service_state(s_state_obj[s_obj.serviceStatus.lower()])
        except AttributeError:
            print(
                "Service ", s_obj, " contains a state that was not defined in database"
            )
    return
