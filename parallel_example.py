'''
Created ab4810 on 23/01/2020

Python script to generate raypath masks for AARM derived datasets.

Read AARM outfiles in form: 

EVT-FILE
1            2             3           4            5           6           7            8           9           10          11          12         13
$DMT_EV_YEAR $DMT_EV_MONTH $DMT_EV_DAY $DMT_EV_HOUR $DMT_EV_MIN $DMT_EV_SEC $DMT_EVNRORG $DMT_EV_LAT $DMT_EV_LON $DMT_EV_DEP $DMT_EV_MMW $DMT_EV_MS $EVNR >> $EVENT_OUT

PHASE-FILE
1            2           3           4           5        6        7    8     9       10     11 12     13 14    15
$DMT_EVNRORG $DMT_EV_LAT $DMT_EV_LON $DMT_EV_DEP $stlat_d $stlon_d $sth $stnr $azim_d $arc_d $p $prett $d $prec $EVNR >> $PHASE_OUT

'''
### Importing various python libraries
# numpy is a useful toolkit for scientific computations
import numpy as np
# matplotlib is a plotting toolkit
# import matplotlib.pyplot as plt

# Obspy is a seismic toolkit
# import obspy
# from obspy.taup import TauPyModel
# from obspy.taup import plot_travel_times
# from obspy.taup import plot_ray_paths
# model = TauPyModel(model='ak135')
# from subprocess import call
# import subprocess
# # import matplotlib
# import sys,glob
# import os.path
# import math
import time
# import matplotlib.pylab as pylab
# params = {'legend.fontsize': 'x-large',
#           'figure.figsize': (16, 10),
#          'xtick.labelsize':'16',
#          'ytick.labelsize':'16'}
# pylab.rcParams.update(params)

from pathlib import Path
home = str(Path.home())


# import multiprocessing
import concurrent.futures


start = time.perf_counter()

def do_something(seconds):
    print(f'Sleeping {seconds} seconds (s)....')
    time.sleep(seconds)
    return f'Done sleeping...{seconds}'



def do_something2(ind):
    global dep
    mask_grid=ind*np.ones(60)+dep
    mask_grid[mask_grid<=53]=0
    do_something(ind)
    return mask_grid

dep_mask_grids=[]
ind = [5, 4, 3, 2, 1] # needs to be iterable
dep=50

# Max workers sets the max number of cores.
with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
    # map returns results in order started...
    results = executor.map(do_something2, ind)
    for res in results:
        dep_mask_grids.append(res)

    print(np.shape(dep_mask_grids))
    print(dep_mask_grids)
    
    
finish= time.perf_counter()

print(f'Finished in {round(finish-start, 2)} seconds(s)')

# change the for loop to def function













# Option 3
#
#
# with concurrent.futures.ProcessPoolExecutor() as executor:
#
#     secs = [5, 4, 3, 2, 1] # needs to be iterable
#     Can also use regular loop rather tahn list comprehension.
#     results = [executor.submit(do_something, sec) for sec in secs]
#
#     for f in concurrent.futures.as_completed(results):
#         print(f.result())
    





#
# # Option 2
# processes=[]
# for _ in range(10): # we are not actualy using the integer in the loop so the underscore is a throw away variable.
#     p=multiprocessing.Process(target=do_something, args=[1.5])
#     p.start()
#     processes.append(p)
#
# for process in processes:
#     process.join()
    

