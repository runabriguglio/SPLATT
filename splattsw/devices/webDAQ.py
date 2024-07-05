from __future__ import print_function
import socket
import sys
import json
import struct
from struct import unpack, calcsize
from time import sleep
from requests import HTTPError, Session, RequestException
import numpy as np
from numpy import array


from SPLATT.splattsw.devices.utility import *
#from IPython.core.release import keywords


# def openwdd(fname):
#     hdr = {}
#     with open(fname,"rb") as wdf:
#         rawData = wdf.read(564)
#         hdr['version'] = int.from_bytes(rawData[0:4], 'little')
#         hdr['size'] = int.from_bytes(rawData[4:8], 'little')
#         hdr['nchannels']= int.from_bytes(rawData[8:12], 'little')
#         hdr['scan_rate']= int.from_bytes(rawData[12:20], 'little')
#         hdr['start_time']= int.from_bytes(rawData[20:28], 'little')
#         hdr['timezone']= rawData[28:44].decode("utf-8")
#         json_hdr_size =  int.from_bytes(rawData[560:564], 'little')
#         jsonRaw = wdf.read(json_hdr_size)
#         hdr['json_hdr']=json.loads(jsonRaw)
#         ndata = hdr['json_hdr']['jobDescriptor']['acquisition']['stopTrigger']['sampleCount']
#         data_it = struct.iter_unpack('<d', wdf.read(ndata*hdr['nchannels']*8)) #4 because double precision 64 bit\n",
#         tmp = np.asarray(list(data_it), dtype='double')
#         data=tmp.reshape(int(tmp.size/4), 4)
#         data = data.T
#     return data

class WebDAQ(object):
    
    def __init__(self, host_address="193.206.154.195"): #, **keywords):
        # Start a session so the password (if any) will be used throughout
        self._s = Session()
        # user = keywords.pop('user','admin')
        # pword = keywords.pop('pwd','admin')
        user = 'admin'
        pword = 'admin'
        self._s.auth = (user, pword)
        self._host_address = host_address
        
        #self._max_read_sample_count = keywords.pop('samples',10000)
        #self._sample_size_in_bytes = keywords.pop('size',8)
        
        self._base_url = ''
        self._initial_schedule_status_code = ''
        
    def connect(self):
        try:
            # Verify the host address is valid.
            host_address = socket.gethostbyname(self._host_address)

            # Host address is valid.
            api_version_response = self._s.get('http://' + host_address + '/api/version')
            api_version_response.raise_for_status()
            api_version_json = api_version_response.json()

            # Get the version of the REST API and create the base URL address.
            version = api_version_json['apiVersion']
            self._base_url = 'http://' + host_address + '/api/' + version

            # Get the device name from the system information.
            system_info_response = self._s.get(self._base_url + '/system/info')
            system_info_response.raise_for_status()
            system_info_json = system_info_response.json()
            dev_name = system_info_json['name']

            print('\nSelected WebDAQ device: ' + dev_name + ' (' + host_address + ')\n')
            print('  Endpoints demonstrated:')
            print('      /schedule/status/, json={"run": True}')
            print('      /schedule/jobs/{job_name}/samples/{sample_index}/{sample_count}/bin')
            print('\n')

            # Get the descriptor for the current schedule.
            schedule_descriptor_endpoint = self._base_url + '/schedule/descriptor'
            schedule_descriptor_response = self._s.get(schedule_descriptor_endpoint)
            schedule_descriptor_response.raise_for_status()
            schedule_descriptor_json = schedule_descriptor_response.json()

            # Get the initial status of the schedule.
            schedule_status_json = get_schedule_status_json(self._s, self._base_url)
            initial_schedule_status_code = schedule_status_json['statusCode']

            # Is the schedule valid?
            valid = is_schedule_valid(schedule_descriptor_json)
            if valid is False:
                # Schedule is not valid.
                print('  No schedule has been created.\n\n'
                      '      All endpoints require a schedule to be created by the user.\n\n'
                      '  Connect to the WebDAQ using a browser to create\n'
                      '  a schedule and then re-run this example.\n')
            else:
                # Schedule is valid - list the jobs
                print('  There is a valid schedule containing the following jobs:\n')
                self._jobs = schedule_descriptor_json['jobs']
                number_of_jobs_in_schedule = len(self._jobs)
                for job_number in range(number_of_jobs_in_schedule):
                    print('      ', job_number, ': ', self._jobs[job_number], sep='')

        except HTTPError as e:
            print('\n', e)

            content = e.response.json()
            print('\n', json.dumps(content, indent=4), sep='')

            # leave the schedule the way we found it ... if it was already running, then do not stop it
            if initial_schedule_status_code != ScheduleStatus['running'] and \
                    initial_schedule_status_code != ScheduleStatus['waiting']:
                s.post(self._base_url + '/schedule/status/', json={"run": False})


    def get_schedule_status(self):
        schedule_status_json = get_schedule_status_json(self._s, self._base_url)
        print('Schedule status', ':',  get_job_status_key(schedule_status_json['statusCode']),sep='')
    
    def get_jobs_status(self):
        number_of_jobs_in_schedule = len(self._jobs)
        for job_number in range(number_of_jobs_in_schedule):
            job_status_json = get_job_status_json(self._s, self._base_url, self._jobs[job_number])
            print('      ', job_number, ': ', self._jobs[job_number], ':',  get_job_status_key(job_status_json['statusCode']),sep='')

    # def stop_all_jobs(self):
    #     jobs = schedule_descriptor_json['jobs'] # was self._
    #     number_of_jobs_in_schedule = len(jobs)
    #     for job_number in range(number_of_jobs_in_schedule):
    #         stop_res = self._s.post(self._base_url + '/schedule/jobs/'+jobs[job_number]+'/status/', json={"stop":True})
    
    def start_schedule(self):
        start_schedule_response = self._s.post(self._base_url + '/schedule/status/', json={"run": True})
        start_schedule_response.raise_for_status()
    
    def stop_schedule(self):
        stop_schedule_response = self._s.post(self._base_url + '/schedule/status/', json={"run": False})
        stop_schedule_response.raise_for_status()

    # def start_job(self,jobname):
    #     start_job_response = self._s.post(self._base_url + '/schedule/jobs/'+jobname+'/status/', json={"run": True})
    #     start_job_response.raise_for_status()
    #     job_status_json = get_job_status_json(self._s, self._base_url, jobname)
        
    def stop_job(self):
        schedule_status_json = get_schedule_status_json(wd._s, wd._base_url)
        jobname = schedule_status_json['currentJobName']
        stop_res = self._s.post(self._base_url + '/schedule/jobs/'+jobname+'/status/', json={"stop":True})
        stop_res.raise_for_status()
    
    def read_data_current_job(self, **keywords):
        
        max_read_sample_count = keywords.pop('samples',10000)
        sample_size_in_bytes = keywords.pop('size',8)
        
        schedule_status_json = get_schedule_status_json(self._s, self._base_url)
        current_job_name = schedule_status_json['currentJobName']
        if current_job_name=='':
            raise RuntimeError('No Job running')
        # Display the current job name.
        print(' Job: ', current_job_name, 'data')

        # Use the job descriptor to get the number of channels in the job.
        job_descriptor_response = self._s.get(self._base_url + '/jobs/' + current_job_name + '/descriptor')
        job_descriptor_response.raise_for_status()
        job_descriptor_json = job_descriptor_response.json()
        channels = job_descriptor_json['channels']
        number_of_channels = len(channels)

        # Display the header for the data table
        print('  Sample number    ', end='')
        for chan, item in enumerate(channels):
            print('      Channel ', item['number'], ' ', sep='', end='')
        print('')

        last_iteration_index = -1
        read_index = 0

        # Read and display the data
    

        # Get the status for the current job
        job_status_json = get_job_status_json(self._s, self._base_url, current_job_name)
        job_status_code = job_status_json['statusCode']

        # Calculate the number of available samples to read.
        samples_acquired = job_status_json['samplesAcquired']
        current_iteration_index = job_status_json['iterationIndex']

        # Reset the samples already read since this is a new instance of the job.
        if current_iteration_index != last_iteration_index:
            read_index = 0

        last_iteration_index = current_iteration_index

        # Read the data.
        data = self._s. get(self._base_url + '/schedule/jobs/' + current_job_name + '/samples/' +
                      str(int(read_index)) + '/' + str(max_read_sample_count) + '/bin')
        data.raise_for_status()

        samples_read = len(data.content) / (sample_size_in_bytes * number_of_channels)

        if samples_read > 0:
            # Convert the array of binary byte data to an array of doubles.
            read_format = str(int(samples_read * number_of_channels)) + 'd'
            read_size = calcsize(read_format)
            numpy_data = array(unpack(read_format, data.content[0: read_size]))

            # Calculate the index to the last sample.
            display_index = read_index + samples_read - 1
            print('\r{:15}     '.format(display_index), end='')

            offset = len(numpy_data) - number_of_channels
            for chan, item in enumerate(channels):
                print('{:12.5f}'.format(numpy_data[int(offset+chan)]), item['unit'], ' ',
                      end='')
                sys.stdout.flush()

            read_index += samples_read
        else:
            print("No samples read.")
            return None

        schedule_status_json = get_schedule_status_json(self._s, self._base_url)

        sleep(0.01)

        if job_status_code == JobStatus['completed'] or \
                job_status_code == JobStatus['stopped'] or \
                job_status_code == JobStatus['jumped']:
            if samples_acquired == read_index:
                current_job_name = schedule_status_json['currentJobName']
                print('')
        elif job_status_code == JobStatus['queued']:
            current_job_name = schedule_status_json['currentJobName']
            print('')

        return numpy_data
    
    
    
