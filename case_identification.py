# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 15:39:32 2017

@author: kurtispinkney
"""

#import packages
import numpy as np
import glob
import sys
import datetime
import time

#read in parsed flash data
files = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_flash_data/2014/flash_*qc.dat')

data = []
FMT = '%H%M%S'

#loop through list of files
casefile = open('/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/2014cases.txt','w')
for file in files:
    #casefile.write('date t1 t2'  + "\n")
    date = file[-12:-6]
    data = np.genfromtxt(file,dtype=[('mystring','S6'),('myfloat','f8'),('yourfloat','f8'),('everyonesfloat','f8')],
                         names=['time','lat','lon','alt'], usecols = (0,1,2,3))
    
    #assign lightning data to variables                     
    timedata = data['time']
    lat = data['lat']
    lon = data['lon']
    alt = data['alt']
    
    #determine values within HYTOP domain
    origin_lat = 34.93055725 #taken from grid origin
    origin_lon = -86.08361053 #taken from grid origin
    lat_km = (lat-origin_lat)*110.926
    lon_km = (lon-origin_lon)*92.15
    alt_km = (alt)/1000.
    
    #just switched from 125 square km radius to 90 square km radius... 
    #year 2012 and 2013 identify cases using 125 square km radius
    lat_parse = lat_km[np.abs(lat_km) < 90]
    lon_parse = lon_km[np.abs(lon_km) < 90]
    
    timedata = timedata[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
    timedata = sorted(timedata)
    lat = lat[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
    lon = lon[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
    alt = alt[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
    
    timerange = np.size(timedata)   
    
    for i in range(timerange-1):
        t1 = timedata[i].decode('utf-8')
        t2 = timedata[i+1].decode('utf-8')
        
#        print(t1, t2)
#        time.sleep(10)
#        t1 = int(round(timedata[i][1]))
#        t2 = int(round(timedata[i+1][1]))
        sys.exit()
        tdelta = datetime.datetime.strptime(t2, FMT) - datetime.datetime.strptime(t1, FMT)
#        tdelta = datetime.timedelta(seconds=t2)-datetime.timedelta(seconds=t1)
        if tdelta.seconds > 1800:
            print(date)
            print(tdelta.seconds)
            print("{} {} {}\n".format(date, t1, t2), file=casefile)
            #casefile.write('{} {} {}\n'.format(date, t1, t2))
#    time.sleep(3)       
    
    
