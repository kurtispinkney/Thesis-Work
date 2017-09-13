#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:26:00 2017

Script to do xlma processing

@author: kurtispinkney
"""
#import packages
import glob
import numpy as np
import sys
import collections
import datetime

#generate a list of filenames
xlmafiles = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/xlma/2016/*')

#identify origin location of Hytop domain
o_lat = 34.93055725
o_lon = -86.08361053

FMT = '%H%M%S'

#loop through xlma files
for file in xlmafiles:
    #generate outfile name time
#    outfiletime = file[83:89]
    outfiletime = file[69:75]
    print(file)
    print(outfiletime)

    #read in xlma data
    data = np.genfromtxt(file,skip_header=68,usecols=(0,1,2,3,7))

    #assign data to variables
    time = data[:,0]
    lat = data[:,1]
    lon = data[:,2]
    flashnum = data[:,4] 

    #identify flashes by looking for data in flashnum variable with counts >~ 10
    flashes = [item for item, count in collections.Counter(flashnum).items() if count > 9]

    #initialize an array that holds the times of the last source associated with each flash
    last_time = np.zeros((np.size(flashes),2))

    #loop through flashes to determine if in domain and determine last flash times
    count = 0
    for i in flashes:
        #find indices of sources associated with flash
        flash_indices = np.where(data[:,4] == i)
        flash_sources = data[flash_indices]

        #find time and lat/lon data of those sources
        time = flash_sources[:,0]
        lat = flash_sources[:,1]
        lon = flash_sources[:,2]      
        
        #convert lat/lon coordinates to km
        lat_km = (lat-o_lat)*110.926
        lon_km = (lon-o_lon)*92.15
        
        #find time and lat/lon data within 125 km x 125 km domain
        time1 = time[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        lat1 = lat[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        lon1 = lon[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        
        #output last source time to second column of last_time array if at least 10 sources are within domain
        if np.size(lat1) >= 10:
            last_time[count][0] = i
            last_time[count][1] = time1[-1]
        count = count + 1
        
    last_time = last_time[~(last_time==0).all(1)]
    timerange = np.size(last_time[:,0]) 
    
    if np.size(flashes) == 1:
        cessation_time = last_time
        cessation_time_flashnum = cessation_time[0,0]
        
    else:
        
        #initiate counter number to loop through flashes to determine if they're 30 min apart 
        #loop through all the flashes
        for j in range(timerange-1):
    #        t1 = str(round(last_time[j][1]).astype(int))
    #        t2 = str(round(last_time[j+1][1]).astype(int))
            t1 = int(round(last_time[j][1]))
            t2 = int(round(last_time[j+1][1]))
    

            #determine the time difference between the flashes
#            tdelta = datetime.datetime.strptime(t2, FMT) - datetime.datetime.strptime(t1, FMT)
            tdelta = datetime.timedelta(seconds=t2)-datetime.timedelta(seconds=t1)
#            print(str(datetime.timedelta(seconds=t2)), last_time[j+1][0])

            #if the time difference between the flashes is 30 min output that data
            if tdelta.seconds > 1800:
                cessation_time = last_time[j+1]
                cessation_time_flashnum = cessation_time[0]
            else:
                cessation_time = last_time[j+1]
                cessation_time_flashnum = cessation_time[0]
    print(cessation_time)
    

    #find all the data associated with the cessation_time variable
    cessation_indices = np.where(data[:,4] == cessation_time_flashnum)
    cessation_sources = data[cessation_indices]
    
    #find time and lat/lon data of those sources
    cessation_time = cessation_sources[:,0]
    cessation_lat = cessation_sources[:,1]
    cessation_lon = cessation_sources[:,2]
    cessation_alt = cessation_sources[:,3]
    cessation_flashnum = cessation_sources[:,4]    
    
    #convert lat/lon coordinates to km
    cessation_lat_km = (cessation_lat-o_lat)*110.926
    cessation_lon_km = (cessation_lon-o_lon)*92.15
    cessation_alt_km = cessation_alt/1000

    #find time and lat/lon data within 125 km x 125 km domain
    cessation_time1 = cessation_time[(np.abs(cessation_lon_km) < 90) & (np.abs(cessation_lat_km) < 90)]
    cessation_lat1 = cessation_lat[(np.abs(cessation_lon_km) < 90) & (np.abs(cessation_lat_km) < 90)]
    cessation_lon1 = cessation_lon[(np.abs(cessation_lon_km) < 90) & (np.abs(cessation_lat_km) < 90)]
    cessation_alt1 = cessation_alt[(np.abs(cessation_lon_km) < 90) & (np.abs(cessation_lat_km) < 90)]
    cessation_flashnum1 = cessation_flashnum[(np.abs(cessation_lon_km) < 90) & (np.abs(cessation_lat_km) < 90)]
    
    #create an array of times displayed in UTC format
    cessation_UTCtime1 = []
    for j in range(np.size(cessation_time1)):
        cessation_UTCtime1 = np.append(cessation_UTCtime1,str(datetime.timedelta(seconds=cessation_time1[j])))
    
    #combine cessation data into a multi-dimensional array
    cessation_ndim = np.zeros(cessation_time1.size, dtype=[('var1', float), ('var2', float),
                                                           ('var3',float),('var4',float),('var5','U32'),('var6',float)])
    #fill multi-dimensional array
    cessation_ndim['var1'] = cessation_time1
    cessation_ndim['var2'] = cessation_lat1
    cessation_ndim['var3'] = cessation_lon1
    cessation_ndim['var4'] = cessation_alt1
    cessation_ndim['var5'] = cessation_UTCtime1
    cessation_ndim['var6'] = cessation_flashnum1
                  
    #output sources to a text file
    np.savetxt('/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_xlma_data/2016/xlma_'+outfiletime+'.txt',
               cessation_ndim, fmt='%f %f %f %f %s %f')