#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 11:31:26 2017

Script written to pre-process radar data to remove duplicate elevation angles
and sails scans

Written By: Sarah Stough, University of Alabama in Huntsville

@author: Sarah Stough
"""
import pyart
import numpy as np
import numpy.ma as ma
import copy

#To run these functions, call:$ new_radar = cut_sails(original_radar)

#Define a function that will remove individual 1-D or 2-D arrays of data corresponding to each removed elevation angle
#This can be thought of as a helper function for the "cut_sails" function
#This function takes the index of the starting and ending values of the data in each elevation angle to be cut, the 
#array to operate on, and indication of the number of dimensions in the array
#Note: truc_concat short for truncate then concatenate remaining data
def trunc_concat(start_ind, end_ind, data, array_dim=2):
    if array_dim==2:
        h1 = data[0:start_ind,:]
        h2 = data[end_ind:,:]
    if array_dim==1:
        h1 = data[0:start_ind]
        h2 = data[end_ind:]
    new = ma.concatenate((h1,h2),axis=0)
    return new

#Function to remove all data corresponding to SAILS or undesired hi-res scans - pass a radar object to this function
#For SAILS, tilts are usually (0-based) 1, 3, 5, 8, 9 but it's necessary to remove
#data in reverse order so the count/index list remains correct and can be adjusted on the fly
def cut_sails(original_radar):
    #Make a hard copy of the radar object being passed - ensures that if desired, this original import object won't be altered 
    oradar = copy.deepcopy(original_radar)
    #Find super-res tilts (have 720 elements per sweep with 0.5° beamwidth)
    start_inds = oradar.sweep_start_ray_index['data']
    num_azs = start_inds[1:]-start_inds[0:len(start_inds)-1]
    #SAILS inds have 720 azimuth angles per sweep, isolate those with 720
    sails_inds = np.where(num_azs==720)
    #Keep list of sweeps to drop (every other of bottom 6 tilts; all others above)
    drop_swps = []
    for i in np.arange(len(sails_inds[0])):
        if sails_inds[0][i]<6:
            if sails_inds[0][i]%2 == 1:
                drop_swps.append(sails_inds[0][i])
            if sails_inds[0][i]>6:
                drop_swps.append(sails_inds[0][i])
    drop_swps = np.array(list(reversed(drop_swps)))
	#Extract and reconstruct data around each elevation angle - start in reverse to preserve earlier index values
    for a in drop_swps:
        #Get azimuth index values to cut for the specified elevation angle index
        si = oradar.sweep_start_ray_index['data'][a]
        ei = oradar.sweep_end_ray_index['data'][a]+1
        #1 Truncate the data for reflectivity, phi, zdr, and rho based on azimuth. 
        oradar.fields['reflectivity']['data'] = trunc_concat(si,ei,oradar.fields['reflectivity']['data'])
        oradar.fields['differential_phase']['data'] = trunc_concat(si,ei,oradar.fields['differential_phase']['data'])
        oradar.fields['cross_correlation_ratio']['data'] = trunc_concat(si,ei,oradar.fields['cross_correlation_ratio']['data'])
        oradar.fields['differential_reflectivity']['data'] = trunc_concat(si,ei,oradar.fields['differential_reflectivity']['data'])
        #2 Truncate data based only on the azimuth array (1D, 9720-element arrays)
        oradar.azimuth['data'] = trunc_concat(si,ei,oradar.azimuth['data'],array_dim=1)
        oradar.elevation['data'] = trunc_concat(si,ei,oradar.elevation['data'],array_dim=1)
        oradar.time['data'] = trunc_concat(si,ei,oradar.time['data'],array_dim=1)
        oradar.instrument_parameters['unambiguous_range']['data']= trunc_concat(si,ei,oradar.instrument_parameters['unambiguous_range']['data'], array_dim=1)
        oradar.instrument_parameters['nyquist_velocity']['data']= trunc_concat(si,ei,oradar.instrument_parameters['nyquist_velocity']['data'], array_dim=1)    
        #3 Truncate data based on elevation index (1D, 19-element arrays)
        oradar.fixed_angle['data'] = ma.concatenate((oradar.fixed_angle['data'][0:a],oradar.fixed_angle['data'][a+1:]),axis=0)
        oradar.sweep_end_ray_index['data'] = ma.concatenate((oradar.sweep_end_ray_index['data'][0:a],oradar.sweep_end_ray_index['data'][a+1:]),axis=0)
        oradar.sweep_start_ray_index['data'] = ma.concatenate((oradar.sweep_start_ray_index['data'][0:a],oradar.sweep_start_ray_index['data'][a+1:]),axis=0)
        oradar.sweep_mode['data'] = ma.concatenate((oradar.sweep_mode['data'][0:a],oradar.sweep_mode['data'][a+1:]),axis=0)
        oradar.sweep_number['data'] = ma.concatenate((oradar.sweep_number['data'][0:a],oradar.sweep_number['data'][a+1:]),axis=0)
        #4 Adjust sweep start, end arrays based on number of elements removed
        diff = ei-si
        oradar.sweep_start_ray_index['data'][a:]=oradar.sweep_start_ray_index['data'][a:]-diff
        oradar.sweep_end_ray_index['data'][a:]=oradar.sweep_end_ray_index['data'][a:]-diff
    #5 Adjust the data that relies on count of azimuth angles or count of elevation angles
    oradar.nsweeps = len(oradar.fixed_angle['data'])
    oradar.nrays = len(oradar.azimuth['data'])
    #Pass a new, adjusted radar object back from the function
    return oradar