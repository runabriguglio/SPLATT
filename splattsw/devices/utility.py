###
# Dictionary for JobStatus values.
###
JobStatus = {"queued": 1, "started": 2, "waitingForTrig": 3, "acquiring": 4, "completed": 5, "stopped": 6,
             "canceled": 7, "jumped": 8, "error": 9, "timeout": 10}
###
# Dictionary for ScheduleStatus values.
###
ScheduleStatus = {"empty": 0, "waiting": 1, "running": 2, "completed": 3, "stopped": 4, "error": 5, "timeout": 6}


def get_job_status_key(value):
    for status, code in JobStatus.items():
        if code == value:
            return status
    return ''

def get_schedule_status_key(value):
    for status, code in ScheduleStatus.items():
        if code == value:
            return status
    return ''


def is_schedule_valid(schedule_descriptor_json):
    """
    Verifies whether or not the schedule is valid.

    Args:
        schedule_descriptor_json: A JSON object containing the current schedule descriptor.

    Returns:
        A flag to indicate whether or not the schedule is valid.
    """

    # Assume the schedule is valid
    valid = True

    # Does the schedule have jobs?
    jobs = schedule_descriptor_json['jobs']
    if len(jobs) == 0:
        valid = False

    return valid


def create_schedule_waiting_msg(schedule_descriptor_json):
    """
    Creates a message to describe what the schedule is waiting for.

    Args:
        schedule_descriptor_json: A JSON object containing the current schedule descriptor.

    Returns:
        A string indicating why the schedule is waiting.
    """
    start = schedule_descriptor_json['start']
    repeat = schedule_descriptor_json['repeat']
    msg = '\n  The schedule has been started, but is waiting for the\n' \
          '  following condition to run a job.\n'
    if start['type'] in 'immediate':
        msg += '\n     Repeat = ' + str(repeat['interval'])
    else:
        msg += '\n     Waiting for ' + start['type']

    return msg


def create_schedule_idle_msg(schedule_status_json):
    """
    Creates a message indicating that the schedule is idle.

    Args:
        schedule_status_json: A JSON object containing the current schedule status.

    Returns:
        A string indicating that the schedule is idle.
    """
    msg = '\n  The schedule is: ' + schedule_status_json['status'] + ' ... there is no available job status.\n'
    msg += '\n  The schedule can be started from a browser or another example'
    msg += '\n  program. This example can be run repeatedly to monitor status.\n\n'

    return msg


def create_job_waiting_msg(job_descriptor_json):
    """
    Creates a message to describe what the job is waiting for.

    Args:
        job_descriptor_json: A JSON object containing the current job descriptor.

    Returns:
        A string indicating why the job is waiting.
    """
    acquisition = job_descriptor_json['acquisition']
    start_trigger = acquisition['startTrigger']
    msg = '\n  The job has been started, but is waiting for the\n' \
          '  following condition.\n'

    msg += '\n      Waiting for ' + start_trigger['type'] + '\n'

    return msg

def get_job_status_json(s, base_url, job_name):
    job_status_response = s.get(base_url + '/schedule/jobs/' + job_name + '/status')
    job_status_response.raise_for_status()
    job_status_json = job_status_response.json()

    return job_status_json


# def get_schedule_descriptor_json(s, base_url):
#     schedule_descriptor_endpoint = base_url + '/schedule/descriptor'
#     schedule_descriptor_response = s.get(schedule_descriptor_endpoint)
#     schedule_descriptor_response.raise_for_status()
#
#     return schedule_descriptor_response


def get_schedule_status_json(s, base_url):
    schedule_status_response = s.get(base_url + '/schedule/status')
    schedule_status_response.raise_for_status()
    schedule_status_json = schedule_status_response.json()

    return schedule_status_json


