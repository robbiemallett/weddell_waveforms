from netCDF4 import Dataset
import datetime
import pandas as pd
import numpy as np

def get_code_from_name(x):
    
    s = x.split(' ')
    
    code = f'D{s[1]} SP{x[5]}'
    
    return(code)

def speed(density,form):
    """ Return the radar speed from dry snow density
    density can be in g/cm3 or kg/m3
    form can be 'Ulaby','Hallikainen','Tiuri'
    string can be 'speed', 'factor' or 'wrongfactor'
    """
    c = 3e8
    if density > 10:
        density = density/1000
    speed_dict = {"Ulaby" : c*((1+0.51*density)**(-1.5)),
                  "Hallikainen":c*((1+1.91*density)**(-0.5)),
                  "Tiuri":c*((1+1.7*density + 0.7*density**2)**(-0.5))}
    return(speed_dict[form])

def get_time_ticks(times):
    seconds = [t.second for t in times]
    hours = [t.hour for t in times]
    minutes = [t.minute for t in times]
    time_ticks = [f'{h}:{str(m).zfill(2)}:{str(s).zfill(2)}' for h,m,s in zip(hours, minutes, seconds)]
    return time_ticks
    

def prepare_data(data_dict):
    all_data = {'ku':{},'ka':{}}

    for freq in ['ku','ka']:

        data_list = []
        time_ticks = []
        times = []

        for key, item in data_dict[freq].items():

            data_list.append(item['data'])

            time_ticks += item['time_ticks']
            times += item['times']

        full = np.concatenate(data_list,axis=1)

        all_data[freq] = {'full_data':full,
                          'time_ticks':time_ticks,
                          'times':times,
                         }
    
    return all_data
    
def get_time_index(time,times):
    
    deltas = np.array(times) - time

    secs = np.array([d.seconds for d in deltas])
    ms = np.array([d.microseconds for d in deltas])
    
    secs = secs + (ms * 1e-6)
    
    index = np.argmin(secs)
    
    return index


def get_range_index(input_range,ranges):
    
    deltas = np.array(ranges) - input_range
    
    index = np.argmin(np.abs(deltas))
    
    return index

def prepare_dicts(pit,pol,data_dir):

    data_dict = {'ku':{},'ka':{}}

    freqranges = {}

    for fnames, freq in zip([pit['ku_f_names'], pit['ka_f_names']],
                            ['ku','ka']):
        for f in fnames:

#             d = Dataset(pit['dir']+'/'+f)

            with Dataset(data_dir+'/'+f) as d:
                data = np.asarray(d[f'{pol}_power_decon0'])
                times = [datetime.datetime(1970,1,1) + datetime.timedelta(seconds=float(s)) for s in d['start_time']]
                ranges = np.asarray(d['range'])
                time_ticks = get_time_ticks(times)

            data_dict[freq][f] = {'data':data,
                            'times':times,
                            'ranges':ranges,
                            'time_ticks':time_ticks}

            freqranges[freq]=ranges
            
    return (freqranges, data_dict)

