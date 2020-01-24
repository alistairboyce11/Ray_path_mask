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
import matplotlib.pyplot as plt

# Obspy is a seismic toolkit
import obspy
from obspy.taup import TauPyModel
from obspy.taup import plot_travel_times
from obspy.taup import plot_ray_paths
model = TauPyModel(model='ak135')
from subprocess import call
import subprocess
import matplotlib
import sys,glob
import os.path
import math
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16, 10),
         'xtick.labelsize':'16',
         'ytick.labelsize':'16'}
pylab.rcParams.update(params)

from pathlib import Path
home = str(Path.home())

print('Imports complete.....')

def calc_mantle_fresnel_zones(dep_ints=np.array([0,50]), P_wave_period=1):

    #calculate vp velocity and fzhw at each depth interval
    dep_vel_fzhw=np.zeros((len(dep_ints),3))
    for i in range(len(dep_ints)):
        vpref=float(model.model.s_mod.v_mod.evaluate_below(dep_ints[i],'p'))
        fzhw=float(np.sqrt(((P_wave_period*vpref)*(1/3)+dep_ints[i])**2.-dep_ints[i]**2.))
        
        dep_vel_fzhw[i,0]=float(dep_ints[i])
        dep_vel_fzhw[i,1]=vpref
        dep_vel_fzhw[i,2]=fzhw

    return dep_vel_fzhw

def calc_pierce_points(file_in_loc=home+'/Google_Drive/GITHUB_AB/3D_model_info', filename='BBAFRP20_phase_SUMMARY.txt',out_file='PP_FZ_OUTFILE.txt',turn_depth_factor=0.8):

    ######## Check input parameters ###############################
    if not isinstance(file_in_loc, str):
        print('Bad specification of file_in_loc - must be string')
        sys.exit()
    if not os.path.exists(file_in_loc):
        print('Bad specification of file_in_loc - does not exist!')
        sys.exit()
    if not isinstance(filename, str):
        print('Bad specification of filename - must be string')
        sys.exit()
    if not isinstance(out_file, str):
        print('Bad specification of out_file - must be string')
        sys.exit()
    if os.path.exists(file_in_loc + '/' + out_file):
        print('Bad specification of out_file :' + str(file_in_loc + '/' + out_file))
        try:
            print('deleting old out-file: '+str(file_in_loc + '/' + out_file)+'...')
            os.remove(file_in_loc + '/' + out_file)
        except:
            print('Failed to delete '+str(file_in_loc + '/' + out_file))
            sys.exit()
    pp_fzhw_file = open(file_in_loc + '/' + out_file, 'w')
        
    ##############################################################

    # Make an array with interogation depths then convert to string for passing to tauP
    dep_ints_string=''
    dep_ints=np.arange(0,2891,50)
    for i in range(len(dep_ints)):
        dep_ints_string=dep_ints_string + str(dep_ints[i]) + ','
    dep_ints_string=dep_ints_string[:-1]
    
    
    # Calculate the fresnel zone widths
    dep_vel_fzhw=calc_mantle_fresnel_zones(dep_ints=dep_ints, P_wave_period=1)
    
    print('Starting ray calcs...')
    # param_list=[]
    data_file = open(file_in_loc + '/' + filename,'r')
    data_list=data_file.readlines()   #.strip().split()
    for line in data_list:
        params=line.split()
        if len(params) != 15:
            print('Something wring with input file')
            sys.exit()
        evlat=float(params[1])
        evlon=float(params[2])
        evdep=float(params[3])
        stlat=float(params[4])
        stlon=float(params[5])
        arc_d=float(params[9])
        print('Calculating pierce-points for ray: evlat: '+str(evlat)+', evlon: '+str(evlon)+', evdep: '+str(evdep)+', stlat: '+str(stlat)+', stlon: '+str(stlon)+', arc_d: '+str(arc_d)
        # Calculate ray turning point
        test1 = ['taup_pierce -mod ak135 -h '+str(evdep)+' -ph P -nodiscon -deg '+str(arc_d)+' -turn']
        # Run test1 in terminal
        out1 = subprocess.check_output(test1,shell=True,universal_newlines=True)
        turn = out1.split('\n')[1]
        turn_dep=float(turn.split()[1])
        # print('Turning depth is: '+str(turn_dep))
        # calc PP in all depth intervals in every 50km intervals
        
        test2 = ['taup_pierce -mod ak135 -h '+str(evdep)+' -ph P -pierce '+str(dep_ints_string)+' -nodiscon -sta '+str(stlat)+' '+str(stlon)+' -evt '+str(evlat)+' '+str(evlon)]
        # Run test in terminal
        out2 = subprocess.check_output(test2,shell=True,universal_newlines=True)
        t = out2.split('\n')
        
        if len(t) > 1:
            for i in range(1,len(t)-1):
                u=t[i].split()
                
                if len(u) == 5:
                    pp_dist=float(u[0])
                    pp_dep=float(u[1])
                    pp_time=float(u[2])
                    pp_lat=float(u[3])
                    pp_lon=float(u[4])
                    
                    # Find the correct FZHW for the depth using index
                    index=int(np.where(dep_vel_fzhw[:,0] == pp_dep)[0])
                    pp_fzhw=float(dep_vel_fzhw[index,2])
                    
                    # Only print pp_deps less than some factor of the turning depth.
                    if pp_dep < turn_dep*turn_depth_factor:
                        pp_fzhw_file.write("%f %f %f %f %f %f %f %f %f %f\n" %(evlat, evlon, evdep, stlat, stlon, arc_d, pp_dep, pp_lat, pp_lon, pp_fzhw))
    
    pp_fzhw_file.close()
    print ('Finished writing pierce - fresnel zone file')
    return


def haversine(lat1, long1, lats2, longs2, depth):
    """
    Input depth in km
    Calculate the distance between two points in earth output in km
    """
    d=[]

    for i in range(len(lats2)):
        lat2=lats2[i]
        long2=longs2[i]
        radius = 6371.e3 - depth*1.e3  # m
        dLat = math.radians(lat2 - lat1)
        dLong = math.radians(long2 - long1)

        a = (math.sin(dLat / 2) ** 2 +
             math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dLong / 2) ** 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d.append(radius * c/1000) # div 1000 for kms
    return np.array(d)

def create_mask_grids(file_in_loc=home+'/Google_Drive/GITHUB_AB/3D_model_info', PP_FZ_file_in='PP_FZ_OUTFILE.txt',mask_grid_file='TEST_MASK',lon_min=0,lon_max=10,lat_min=0,lat_max=10,grid_d_lat_lon=0.25,fzmf=1):
    '''
    Here we want to calculate the mask grid
    specify file-in_loc, file-in filename, output mask_grid_file, lat/lon bounds and spacial sampling interval
    fzmf = 'Fresnel zone multiplication factor' Used to broaden the footprint of each fresnel zone based mask.
    '''
    
    ######## Check input parameters ###############################
    if not isinstance(file_in_loc, str):
        print('Bad specification of file_in_loc - must be string')
        sys.exit()
    if not isinstance(PP_FZ_file_in, str):
        print('Bad specification of filename - must be string')
        sys.exit()
    if not os.path.exists(file_in_loc+'/'+PP_FZ_file_in):
        print('Bad specification of file_in_loc + PP_FZ_file_in - does not exist!')
        sys.exit()
    if not isinstance(mask_grid_file, str):
        print('Bad specification of out_file - must be string')
        sys.exit()
        
    ##############################################################

    # First make the outfile
    dep_mask_filename=file_in_loc + '/' + mask_grid_file+'.txt'
    if os.path.exists(dep_mask_filename):
        print('Bad specification of out_file :' + str(dep_mask_filename))
        try:
            print('deleting old out-file: '+str(dep_mask_filename)+'...')
            os.remove(dep_mask_filename)
        except:
            print('Failed to delete '+str(dep_mask_filename))
            sys.exit()
    dep_mask_file = open(dep_mask_filename, 'w')

    ###########################################################

    print('Starting grid calcs...')
    
    data_file = open(file_in_loc+'/'+PP_FZ_file_in,'r')
    d_vector=data_file.read().strip().split()
    d_vector = [float(x) for x in d_vector]
    data_list=np.reshape(d_vector,(-1,10)) # 10 is number of columns in infile
    
    # Make an array with interogation depths, lats and lons
    dep_ints=np.arange(0,2891,50)
    g_lats=np.arange(lat_min,lat_max+grid_d_lat_lon,grid_d_lat_lon)
    g_lons=np.arange(lon_min,lon_max+grid_d_lat_lon,grid_d_lat_lon)

    grid_points=[]
    lat_points=[]
    lon_points=[]
    for i in range(len(g_lats)):
        for j in range(len(g_lons)):
            lat_points.append(g_lats[i])
            lon_points.append(g_lons[j])
            # grid_points.append([g_lats[i],g_lons[j]])
    # grid_points=np.array(grid_points)
    lon_points=np.array(lon_points)
    lat_points=np.array(lat_points)

    
    #Loop over the dep interval array
    for dep in dep_ints:
        
        # Make dep array to paste with each mask
        dep_points=dep*np.ones(len(lat_points))
        
        # Find rows of data_list that correspond to correct depth.
        dep_index=np.where(data_list[:,6]==dep)[0]
        num_pp=len(dep_index)
        print('Found '+str(num_pp)+' pierce points at '+str(dep)+'km depth..')
        # For each selected row of data_list
        if num_pp>1:
            dep_mask_grids=[]
            for ind in dep_index:
                # use pierce point and Fresnel zone
                pp_lat=data_list[ind,7]
                pp_lon=data_list[ind,8]
                fzhw  =data_list[ind,9]
                # Calculate distances to all points on the grid
                grid_to_pp_dist=haversine(pp_lat, pp_lon, lat_points, lon_points, dep)
                # Make a mask grid "on" everywhere - i.e. equal to one
                mask_grid=np.ones(len(lat_points))
                # Convert all points in mask grid less than specified fresnel zone buffer to zero - mask off.
                mask_grid[grid_to_pp_dist<=(fzhw*fzmf)]=0
                # Add these mask grids to the array dep_mask_grid for each pp
                dep_mask_grids.append(mask_grid)
            # Convert the output into an array with rows: len(lat_points), columns: num_pp
            dep_mask_grids=np.transpose(np.array(dep_mask_grids))
            # Sum over all rows to flatten mask for each depth:
            dep_mask_grid_out=dep_mask_grids.sum(axis=1)
            dep_mask_grid_out[dep_mask_grid_out<num_pp]=0
            dep_mask_grid_out[dep_mask_grid_out==num_pp]=1
        
            print('Writing to file: depth: '+str(dep)+'km\n')
            for i in range(len(lat_points)):
                dep_mask_file.write('{0:2f} {1:3f} {2:4f} {3:5f}\n'.format(dep_points[i], lat_points[i], lon_points[i], dep_mask_grid_out[i]))
    dep_mask_file.close()



# calc_pierce_points(file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info',
#                     filename='BBAFRP20_phase_SUMMARY.txt',out_file='BBAFRP20_phase_SUMMARY_PP_FZ.txt',
#                     turn_depth_factor=0.8)
# create_mask_grids(file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info',
#                     PP_FZ_file_in='BBAFRP20_phase_SUMMARY_PP_FZ.txt',mask_grid_file='BBAFRP20_MASK',
#                     lon_min=-24,lon_max=64,lat_min=-44,lat_max=44,grid_d_lat_lon=0.25,fzmf=1)

# Test implementation
calc_pierce_points(file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info',
                    filename='ATS_50_RAYS.txt',out_file='ATS_50_RAYS_PIERCE.txt',
                    turn_depth_factor=0.8)
                    
                    
create_mask_grids(file_in_loc='/Users/ab4810/Google_Drive/GITHUB_AB/3D_model_info',
                    PP_FZ_file_in='ATS_50_RAYS_PIERCE.txt',mask_grid_file='TEST_ATS_50_MASK',
                    lon_min=-24,lon_max=64,lat_min=-44,lat_max=44,grid_d_lat_lon=0.25,fzmf=10)


