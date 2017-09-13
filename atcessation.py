#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 18:02:47 2017

Script to do analysis at time of cessation

@author: kurtispinkney
"""

#==============================================================================
#Import packages
#==============================================================================
import pyart 
import numpy as np
import numpy.ma as ma
import glob
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as ticker
import mpl_toolkits.axisartist as AA
import sys
import datetime
from scipy.spatial import ConvexHull
import scipy.stats as stats
from skewt import SkewT


#==============================================================================
#Calculation flags
#==============================================================================
xlma = 1
nalma = 0
box = 1
prop = 0
convexhull = 0
cross = 0
vppr = 0
composite = 1
gradient = 0

#==============================================================================
#Read in data 
#==============================================================================
vertfiles = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/atcessation/2015/*')
#vertfiles = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/postpriorcessation/2012/0518test20120519_000129_def.nc')

#read in bounding box info
boxfile = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/atcessation/15cessation.csv'
boxpoints = np.genfromtxt(boxfile,skip_header=1,delimiter=',')
count = 0 

#initialize lists for gradients
dZh_xlist = []
dZh_ylist = []
dZh_zlist = []
magnitudelist = []

for vertfile in vertfiles:
    
    print(vertfile)
    #assign strings to variables for file in/out processes
    flash_date = vertfile[-20:-14]
    file_time = vertfile[-13:-7]
    date_time = vertfile[-22:-7]
    snd_date = vertfile[-22:-14]
    
    #read in higher vertical resolution gridded data
    vertgrid = pyart.io.read_grid(vertfile,exclude_fields=None)
    
    #assign reflectivity data to a variable
    vertgrid.fields['reflectivity']['data'] = ma.masked_less_equal(vertgrid.fields['reflectivity']['data'],10)
    vertrefl = vertgrid.fields['reflectivity']['data']
    
    #assign HID data to a variable
    vertHID = vertgrid.fields['FH']['data']
    
    #read in distances from higher resolution grid center and convert to km
    vert_x_array_km = vertgrid.x['data']/1000
    vert_y_array_km = vertgrid.y['data']/1000
    vert_z_array_km = vertgrid.z['data']/1000
    
    #assign grid lat/lon data to a variable
    vertlat = vertgrid.point_latitude['data']
    vertlon = vertgrid.point_longitude['data']
    vertalt = vertgrid.point_altitude['data']/1000 #in km

    #==============================================================================
    #Sounding data   
    #==============================================================================
    #read in sounding data
#    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2016/'+snd_date+'_snd.txt'
    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2015/'+snd_date+'.txt'
    sounding = SkewT.Sounding(sndfile)
    
    #assign temp and height data to variables
    stemp = sounding.soundingdata['temp']
    stemp = ma.getdata(stemp)
    shght = sounding.soundingdata['hght']
    shght = ma.getdata(shght)
    
    #find values below 15 km (size of grid in vertical)
    below15hght = shght[np.where(shght < 15000)]
    below15temp = stemp[np.where(shght < 15000)]
    
    #find important isotherms
    #0°C
    freezelayeridx = np.abs(0-below15temp).argmin()
    freezelayerhght = below15hght[freezelayeridx]/1000 #km
    freezegrididx = np.abs(freezelayerhght-vert_z_array_km).argmin()
    #-10°C
    neg10idx = np.abs(-10-below15temp).argmin()
    neg10hght = below15hght[neg10idx]/1000 #km
    neg10grididx = np.abs(neg10hght-vert_z_array_km).argmin()
    #-20°C
    neg20idx = np.abs(-20-below15temp).argmin()
    neg20hght = below15hght[neg20idx]/1000 #km
    neg20grididx = np.abs(neg20hght-vert_z_array_km).argmin()
    #-30°C
    neg30idx = np.abs(-30-below15temp).argmin()
    neg30hght = below15hght[neg20idx]/1000 #km
    neg30grididx = np.abs(neg30hght-vert_z_array_km).argmin()
    #-40°C
    neg40idx = np.abs(-40-below15temp).argmin()
    neg40hght = below15hght[neg40idx]/1000 #km
    neg40grididx = np.abs(neg40hght-vert_z_array_km).argmin()

    #==============================================================================
    #Lightning data  
    #==============================================================================
    if xlma == 1:
        #read in lightning data
        xlmafile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_xlma_data/2015/xlma_'+flash_date+'.txt'
#        lxnfile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_xlma_data_test/xlma_120518.txt'
        data = np.genfromtxt(xlmafile,dtype=None,names=['time','lat','lon','alt'], usecols = (0,1,2,3))
   
        #assign data columns to variables
        time = data['time']
        lat = data['lat']
        lon = data['lon']
        alt = data['alt']
        
        #convert lat/lon from ° to km and altitude to km
        lat_km = (lat-vertgrid.origin_latitude['data'])*110.926
        lon_km = (lon-vertgrid.origin_longitude['data'])*92.15
        alt_km = (alt)/1000.
    
    if nalma == 1:
        nalmafile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_flash_data/2016/flash_'+flash_date+'qc.dat'       
#        nalmadata = np.genfromtxt(nalmafile,dtype=None,names=['time','lat','lon','alt'], usecols = (0,1,2,3))
        nalmadata = np.genfromtxt(nalmafile,dtype=[('mystring','S6'),('myfloat','f8'),('yourfloat','f8'),('everyonesfloat','f8')],
                         names=['time','lat','lon','alt'], usecols = (0,1,2,3))
        #assign data columns to variables
        timedata = nalmadata['time']
        lat = nalmadata['lat']
        lon = nalmadata['lon']
        alt = nalmadata['alt']
        
        #convert lat/lon from ° to km and altitude to km
        lat_km = (lat-vertgrid.origin_latitude['data'])*110.926
        lon_km = (lon-vertgrid.origin_longitude['data'])*92.15
        alt_km = (alt)/1000.
        
        #just switched from 125 square km radius to 90 square km radius... 
        #year 2012 and 2013 identify cases using 125 square km radius
        lat_parse = lat_km[np.abs(lat_km) < 90]
        lon_parse = lon_km[np.abs(lon_km) < 90]
        
        timedata = timedata[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        timedata2 = sorted(timedata)
        lat = lat[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        lon = lon[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        alt = alt[(np.abs(lon_km) < 90) & (np.abs(lat_km) < 90)]
        sorteddata = sorted(zip(timedata, lon, lat), key=lambda x: x[0])
#        sys.exit()
        timerange = np.size(timedata2)   
        FMT = '%H%M%S'
        cesstime = []
        new_lat = []
        new_lon = []
        for i in range(timerange-1):
            t1 = timedata2[i].decode('utf-8')
            t2 = timedata2[i+1].decode('utf-8')
            tdelta = datetime.datetime.strptime(t2, FMT) - datetime.datetime.strptime(t1, FMT)
            if tdelta.seconds > 1800:
                new_lat.append(sorteddata[i][2])
                new_lon.append(sorteddata[i][1])
                cesstime.append(timedata2[i])
        new_lat_km = (new_lat-vertgrid.origin_latitude['data'])*110.926
        new_lon_km = (new_lon-vertgrid.origin_longitude['data'])*92.15
#        sys.exit()
        
    #==============================================================================
    #Bounding box tracking   
    #==============================================================================
    if box == 1:
        #generate an array of possible x,y grid point combinations
        grid_points= np.vstack(np.meshgrid(vert_x_array_km, vert_y_array_km)).reshape(2,-1).T
        
        #identify the lower-left and upper-right points of the box
        ll = np.array([boxpoints[count][0], boxpoints[count][1]])  # lower-left
        ur = np.array([boxpoints[count][2], boxpoints[count][3]])  # upper-right
        
        #identify points within the box
        inidx = np.all(np.logical_and(ll <= grid_points, grid_points <= ur), axis=1)
        inbox = grid_points[inidx]
        
        #find indices of grid x,y arrays within bounding box
        xsort = vert_x_array_km.searchsorted(inbox[:,0])
        ysort = vert_y_array_km.searchsorted(inbox[:,1])
     
    #==============================================================================
    #Composite Zh testing
    #==============================================================================
    if composite == 1:
        Zhmax = np.zeros((np.size(vert_x_array_km),np.size(vert_y_array_km)))
        for aa in range(np.size(vert_x_array_km)):
            for bb in range(np.size(vert_x_array_km)):
                Zhmax[aa][bb] = np.amax(vertrefl[:,aa,bb])
        
        #find grid points closest to NSSTC
#        x = (-86.38 - vertgrid.origin_longitude['data'])*92.15
#        y = (34.43 - vertgrid.origin_latitude['data'])*110.926
        
        #plot composite reflectivities
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig, ax = plt.subplots()
        
        #plot reflectivity
        ax.imshow(Zhmax, origin='lower', vmin=00, vmax=65,
                   extent=(-125, 125, -125, 125),cmap = pyart.graph.cm.NWSRef)
        colors = plt.cm.jet(np.linspace(0, 1, len(time)))
        #plot xlma source points
        plt.scatter(lon_km,lat_km,c=colors)
        plt.title(date_time+'Composite Reflectivity')
        plt.xlabel('East-West Distance from Radar')
        plt.ylabel('North-South Distance from Radar')
        
        #set grid intervals for major and minor x-axis grids
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        
        #set grid intervals for major and minor y-axis grids
        ax.yaxis.set_major_locator(MultipleLocator(50))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        
        #plot the grids
        ax.xaxis.grid(True,'minor')
        ax.yaxis.grid(True,'minor')
        ax.xaxis.grid(True,'major',linewidth=2)
        ax.yaxis.grid(True,'major',linewidth=2)
        ax.add_patch(patches.Rectangle((-90, -90),180,180,fill=False,color='purple',linewidth=2))
        
        #plot circle showing NALMA domain
#        circle1 = plt.Circle((x, y), 160, color='r', fill=False)
#        ax.add_artist(circle1)
        
        #plot rectangle showing radar domain
        ax.add_patch(patches.Rectangle((-90, -90),180,180,fill=False))
        
        #plot bounding box
        if box == 1:
            ax.add_patch(patches.Rectangle((ll[0], ll[1]), ur[0]-ll[0], ur[1]-ll[1], fill=False,color='red',linewidth=2))
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/atcessation/PPI/2015/'+date_time+'_composite.png')
        plt.close(fig)
        plt.clf()
        
    #==============================================================================
    #Convex hull analysis    
    #==============================================================================
    if convexhull == 1:
        #identify points to be used in convex hull calculations
        points = np.column_stack((lon_km,lat_km))
        #determine convex hull
        hull = ConvexHull(points)
        
        #find max reflectivity for each grid column
        Zhmax = np.zeros((np.size(vert_x_array_km),np.size(vert_y_array_km)))
        for aa in range(np.size(vert_x_array_km)):
            for bb in range(np.size(vert_x_array_km)):
                Zhmax[aa][bb] = np.amax(vertrefl[:,aa,bb])
        
        #plot composite reflectivity
        ax = plt.subplot(111)
        ax.imshow(Zhmax,
                  origin='lower', vmin=10, vmax=65,extent=(-90, 90, -90, 90)
                  ,cmap = pyart.graph.cm.NWSRef)
        plt.title(date_time + ' Convex Hull')
        plt.xlabel('Distance from radar (km)')
        plt.ylabel('Distance from radar (km)')
        
        #set grid intervals for major and minor x-axis grids
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        
        #set grid intervals for major and minor y-axis grids
        ax.yaxis.set_major_locator(MultipleLocator(50))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        
        #plot the grids
        ax.xaxis.grid(True,'minor')
        ax.yaxis.grid(True,'minor')
        ax.xaxis.grid(True,'major',linewidth=2)
        ax.yaxis.grid(True,'major',linewidth=2)
        
        #plot the last flash sources
        colors = plt.cm.jet(np.linspace(0, 1, len(time)))
        plt.scatter(lon_km,lat_km,c=colors)
        #plot the convex hull
        for simplex in hull.simplices:
            plt.plot(points[simplex, 0], points[simplex, 1], 'b-')
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/atcessation/ConvexHullTest/test'+date_time+'.png')
    
    #==============================================================================
    #Propagation analysis    
    #==============================================================================
    if prop == 1:
        #identify nearest grid points to sources associated with last flash
        source_xindices = np.zeros(np.size(lon_km))
        source_yindices = np.zeros(np.size(lat_km))
        source_zindices = np.zeros(np.size(alt_km))
        
        #initial size of source array
        pre = np.size(source_zindices)
        
        #using the values associated with the grid poins
        vertalt1d = vertalt[:,0,0]
        pre = np.size(source_zindices)

        #determine altitude of each source point
        for ind in range(np.size(alt_km)):
            source_zindices[ind] = (np.abs(vertalt1d-alt_km[ind])).argmin()
       
        #convert to int
        source_zindices = source_zindices.astype(int)

        #determine lon/lat points associated with each source point
        for i in range(np.size(lat)):
            vertlonarray = vertlon[source_zindices[i],:,:]
            vertlatarray = vertlat[source_zindices[i],:,:]
            source_xindices[i] = (np.abs(vertlonarray-lon[i])).argmin()
            source_yindices[i] = (np.abs(vertlatarray-lat[i])).argmin()
        
        #find the associated values in km from radar origin
        flatpointx = np.ndarray.flatten(vertgrid.point_x['data'][0,:,:])
        flatpointy = np.ndarray.flatten(vertgrid.point_y['data'][0,:,:])
        source_xindices2 = np.zeros(np.size(lon_km))
        source_yindices2 = np.zeros(np.size(lat_km))       
        for k in range(np.size(source_xindices)):
            source_xindices2[k] = flatpointx[source_xindices[k]]
            source_yindices2[k] = flatpointy[source_yindices[k]]
        
        #convert to km
        source_xindices2 = source_xindices2/1000
        source_yindices2 = source_yindices2/1000
        
        #determine indices to be used to find HID points
        final_xindices = np.zeros(np.size(source_xindices2))
        final_yindices = np.zeros(np.size(source_yindices2))
        for l in range(np.size(source_xindices2)):
            final_xindices[l] = np.where(source_xindices2[l] == vert_x_array_km)[0]
            final_yindices[l] = np.where(source_yindices2[l] == vert_y_array_km)[0]
            
        #identify HID values associated with each soure
        source_HID = np.zeros(np.size(source_xindices))
        for j in range(np.size(source_xindices)):
            source_HID[j] = vertHID[source_zindices[j]][final_yindices[j]][final_xindices[j]]
            
        #convert nans to a number and all values to integers
        source_HID = np.nan_to_num(source_HID).astype(int)

        #find all points that don't have data
        nodata = np.where(source_HID == 0)[0]

#==============================================================================
#atrocious code ... rewrite later
#==============================================================================
        deletepoints = []
        for n in range(np.size(nodata)):
            surround = []
            #surrounding grid points at the height of source point
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]-1])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]][final_xindices[nodata[n]]-1])
            surround.append(vertrefl[source_zindices[nodata[n]]][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]-1])
            
            #surrounding grid points at the height above source point
            if source_zindices[nodata[n]] == 30:
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
                surround.append(np.nan)
            else:
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]+1])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]][final_xindices[nodata[n]]+1])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]+1])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]-1])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]][final_xindices[nodata[n]]-1])
                surround.append(vertrefl[source_zindices[nodata[n]]+1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]-1])

            #surrounding grid points at the height below source point
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]+1])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]-1][final_xindices[nodata[n]]-1])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]][final_xindices[nodata[n]]-1])
            surround.append(vertrefl[source_zindices[nodata[n]]-1][final_yindices[nodata[n]]+1][final_xindices[nodata[n]]-1])
            
            #indices associated with each surrounding grid point
            surroundindices = np.zeros((np.size(surround),3))
            surroundindices[0] = [source_zindices[nodata[n]],final_yindices[nodata[n]]+1,final_xindices[nodata[n]]]
            surroundindices[1] = [source_zindices[nodata[n]],final_yindices[nodata[n]]+1,final_xindices[nodata[n]]+1]
            surroundindices[2] = [source_zindices[nodata[n]],final_yindices[nodata[n]],final_xindices[nodata[n]]+1]
            surroundindices[3] = [source_zindices[nodata[n]],final_yindices[nodata[n]]-1,final_xindices[nodata[n]]+1]
            surroundindices[4] = [source_zindices[nodata[n]],final_yindices[nodata[n]]-1,final_xindices[nodata[n]]]
            surroundindices[5] = [source_zindices[nodata[n]],final_yindices[nodata[n]]-1,final_xindices[nodata[n]]-1]
            surroundindices[6] = [source_zindices[nodata[n]],final_yindices[nodata[n]],final_xindices[nodata[n]]-1]
            surroundindices[7] = [source_zindices[nodata[n]],final_yindices[nodata[n]]+1,final_xindices[nodata[n]]-1]
            if source_zindices[nodata[n]] == 30:
                print('surroundindices at the top')
            else:
                surroundindices[8] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]]
                surroundindices[9] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]+1]
                surroundindices[10] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]],final_xindices[nodata[n]]+1]
                surroundindices[11] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]+1]
                surroundindices[12] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]]
                surroundindices[13] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]-1]
                surroundindices[14] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]],final_xindices[nodata[n]]-1]
                surroundindices[15] = [source_zindices[nodata[n]]+1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]-1]
            
            surroundindices[16] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]]
            surroundindices[17] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]+1]
            surroundindices[18] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]],final_xindices[nodata[n]]+1]
            surroundindices[19] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]+1]
            surroundindices[20] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]]
            surroundindices[21] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]-1,final_xindices[nodata[n]]-1]
            surroundindices[22] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]],final_xindices[nodata[n]]-1]
            surroundindices[23] = [source_zindices[nodata[n]]-1,final_yindices[nodata[n]]+1,final_xindices[nodata[n]]-1]
            
            surround = np.array(surround)

#==============================================================================
#end of atrocious code to be rewritten            
#==============================================================================
            #if no data in cube around source delete the point else replace its value
            if np.all(np.isnan(surround)):
                deletepoints.append(nodata[n])
            else:
                replace = np.where(surround == np.nanmax(surround))[0][0]
                source_HID[nodata[n]]=vertHID[surroundindices[replace][0]][surroundindices[replace][1]][surroundindices[replace][2]]

        #delete all the source points that don't have data around them
        source_HID=np.delete(source_HID,deletepoints)
        source_zindices=np.delete(source_zindices,deletepoints)
        
        #final size of source array
        post = np.size(source_zindices)
        
        #determine percent change of source array
        change = ((post-pre)/pre)*100
        print(date_time, change)
        #determine vertical array values
        vertvalues = vert_z_array_km[source_zindices.astype(int)]

        #determine bin edge values
        xedges = np.arange(0,np.max(source_HID)+2,1)
        yedges = np.arange(0,np.max(vertvalues)+1,0.5)
        
        #generate a figure
        fig, ax = plt.subplots()
        ax.set_aspect("equal")
        
        #generate histogram and get bin data
        hist, xbins, ybins, im = ax.hist2d(source_HID,vertvalues,bins=(xedges-0.5,yedges-0.25))
        #set x,y ticks
        ax.set(xticks=range(np.size(xedges)), xlim=[-0.5, np.max(source_HID)+0.5])
        ax.set(yticks=np.arange(0,np.max(vertvalues)+0.5,0.5), ylim=[-0.25, np.max(vertvalues)+0.25])
        #edit x axis labels
        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels = ['N/A','DR','RN','IC','AG','WS','VI','LG','HG','HL','BD']
        ax.set_xticklabels(labels,rotation='-45')
        plt.xlabel('HID Type')
        plt.ylabel('Height (km)')
        plt.title(snd_date+' 2D Histogram of Source HID with Height')
        
        #annotate the bins to show the count of each bin
        for i in range(len(xbins)-1):
            for j in range(len(ybins)-1):
                if hist[i,j] > 0:
                    ax.text(xbins[i]+0.5,ybins[j]+0.25, hist[i,j], 
                            color="w", ha="center", va="center")
#        sys.exit()
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/atcessation/HIDsourcesTest/0518test'+date_time+'.png')
#        plt.close(fig)
#        plt.clf()
#    sys.exit()
    
    #==============================================================================
    #Gradient analysis
    #==============================================================================
    if gradient == 1:
        #find indices of grid point associated with initial last flash source 
        initZh_x = (np.abs(vert_x_array_km-lon_km[0])).argmin()
        initZh_y = (np.abs(vert_y_array_km-lat_km[0])).argmin()
        initZh_z = (np.abs(vert_z_array_km-alt_km[0])).argmin()
        
        #determine reflectivity last flash source
        initZh = vertrefl[initZh_z,initZh_y,initZh_x]

        #determine change in reflectivity for different directions
        dx = [vertrefl[initZh_z,initZh_y,initZh_x-1], vertrefl[initZh_z,initZh_y,initZh_x+1]]
        dy = [vertrefl[initZh_z,initZh_y+1,initZh_x], vertrefl[initZh_z,initZh_y-1,initZh_x]] 
        dz = [vertrefl[initZh_z+1,initZh_y,initZh_x], vertrefl[initZh_z-1,initZh_y,initZh_x]] 
        
        #distance to divide by
        xkm = 2
        ykm = 2
        zkm = 2
        
        #determine if neighbors are masked and alter list if so
        for a in range(np.size(dx)):
            if np.ma.is_masked(dx[a]):
                dx[a] = initZh
                xkm = 1
            if np.ma.is_masked(dy[a]):
                dy[a] = initZh
                ykm = 1
            if np.ma.is_masked(dz[a]):
                dz[a] = initZh
                zkm = 1
        
        #determine gradients of scalar field
        dZh_x = (dx[0]-dx[1])/xkm
        dZh_y = (dy[0]-dy[1])/ykm
        dZh_z = (dz[0]-dz[1])/zkm
        
        #determine magnitude of resulting vector field
        magnitude = np.sqrt(dZh_x**2 + dZh_y**2 + dZh_z**2)
        
        #fill lists
        dZh_xlist.append(dZh_x)
        dZh_ylist.append(dZh_y)
        dZh_zlist.append(dZh_z)
        magnitudelist.append(magnitude)
        
#    sys.exit()  

    #==============================================================================
    #VPPR analysis    
    #==============================================================================
    if vppr == 1:
        #initialize array to hold reflectivity values within bounding box at different temperature heights
        freezeZh = np.zeros(np.size(xsort))
        
        #create arrays with Zh values within bounding box
        for i in range(np.size(xsort)):
            freezeZh[i]=vertrefl[freezegrididx][ysort[i]][xsort[i]]
        
        #identify index of max reflectivity within bounding box at 0°C
        freezeZhmaxidx = np.nanargmax(freezeZh)
        
        #identify max reflectivity valeu at freezing level
        freezeZhmax = freezeZh[freezeZhmaxidx]
        #identify reflectivity value at -20°C 
        neg20Zh = vertrefl[neg20grididx][ysort[freezeZhmaxidx]][xsort[freezeZhmaxidx]]
        
        #test for data at both values 
        freezetest = np.ma.is_masked(freezeZhmax)
        neg20test = np.ma.is_masked(neg20Zh)
        while freezetest:
            freezegrididx = freezegrididx - 1
            freezeZhmax = vertrefl[freezegrididx][ysort[freezeZhmaxidx]][xsort[freezeZhmaxidx]]
            freezetest = np.ma.is_masked(freezeZhmax)
        while neg20test:
            neg20grididx = neg20grididx - 1
            neg20Zh = vertrefl[neg20grididx][ysort[freezeZhmaxidx]][xsort[freezeZhmaxidx]]
            neg20test = np.ma.is_masked(neg20Zh)
        
#        Zhlapse = (neg20Zh - freezeZhmax)/(vert_z_array_km[neg20grididx]-vert_z_array_km[freezegrididx])
        Zhlapse = (neg20Zh - freezeZhmax)/3
        height_change = vert_z_array_km[neg20grididx]-vert_z_array_km[freezegrididx]
        print('height change is', height_change)
        print(date_time+ 'Zh lapse rate is ',Zhlapse)
        
    #==============================================================================
    #Plotting 
    #==============================================================================
    if cross == 1:
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance
        
        #Set up the display for the gridded data
        display = pyart.graph.GridMapDisplay(vertgrid)
        
        #set ranges of x axis for cross sections
        northsouth_range = [ll[1],ur[1]]
        eastwest_range = [ll[0],ur[0]]
        vertical_range = [0,15]
        
        #identify mid points of bounding box
        latmid = round((ur[1]+ll[1])/2)
        lonmid = round((ur[0]+ll[0])/2)
        
        #identify point where cross sections are taken (mid points of bounding box sides)
        cross_lat = latmid
        cross_lon = lonmid
        
        cross_lat = lat[-1]
        cross_lon = lon[-1]
        
        #get values for where lat/lon contours will be taken
        latcontour = np.abs(latmid-vertgrid.x['data']/1000).argmin()
        loncontour = np.abs(lonmid-vertgrid.y['data']/1000).argmin()
#        
        latcoord = vertlat[6,latcontour,loncontour]
        loncoord = vertlon[6,latcontour,loncontour]
        
        #intervals in which contours will be drawn
#        mwlevels = np.arange(0.5,10.5,1)
#        milevels = np.arange(0.5,10.5,1)
#        zdrlevels = np.arange(-1,4.5,.5)
#        kdplevels = np.arange(-.5,3.5,.5)
#        refllevels = np.arange(10,70,5)
        
        #variable min/max
        refl_min = 10
        refl_max = 65
        zdr_min = -6
        zdr_max = 10
        
        colors = plt.cm.jet(np.linspace(0, 1, len(time)))
        
        #generate plots
        ax1 = plt.subplot(1,2,1)
        ax1.set_ylim(vertical_range)
        ax1.set_xlim(northsouth_range)
        #Plot the image
#        contours_1a = ax1.contour(vertkdp[:,:,loncontour],kdplevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='m')
#        ax1.clabel(contours_1a, inline=1, fontsize=10)
#        contours_1b = ax1.contour(vertzdr[:,:,loncontour],zdrlevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='k')
#        ax1.clabel(contours_1b, inline=1, fontsize=10)
#        contours_1c = ax1.contour(Mi[:,:,loncontour],milevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='m')
#        ax1.clabel(contours_1c, inline=1, fontsize=10)
        display.plot_longitude_slice('reflectivity', lon=cross_lon,lat=cross_lat, vmin=refl_min, vmax=refl_max,
                                      axislabels = ('North-South Distance from Radar (km)', 'Height (km)'), 
                                      cmap = pyart.graph.cm.NWSRef)
        ax1.scatter(lat_km,alt_km,c=colors)
#        plt.scatter(lon_km,alt_km,c=colors)
        plt.title(date_time + ' Reflectivity')
#        ax1.plot(lat_km2[cessation_idx],alt1[cessation_idx],'ko')
        ax1.axhline(freezelayerhght,c="blue",linestyle='--')
        ax1.text(ll[1]+1,((freezelayerhght)+0.1),'0C',color='b')
        ax1.axhline(neg10hght,c="red",linestyle=':')
        ax1.text(ll[1]+1,((neg10hght)+0.1),'-10C',color='r')
        ax1.axhline(neg20hght,c="purple",linestyle='-.')
        ax1.text(ll[1]+1,((neg20hght)+0.1),'-20C',color='purple')
    
#        sys.exit()
        ax2 = plt.subplot(1,2,2)
        ax2.set_ylim(vertical_range)
        ax2.set_xlim(eastwest_range)
        #Plot the image
#        contours_2a = ax2.contour(vertkdp[:,latcontour,:],kdplevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='m')
#        ax2.clabel(contours_2a, inline=True, fontsize=10)
#        contours_2b = ax2.contour(vertzdr[:,latcontour,:],zdrlevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='k')
#        ax2.clabel(contours_2b, inline=True, fontsize=10)
#        contours_2c = ax2.contour(Mi[:,latcontour,:],milevels,linewidths=1.5,antialiased=True,extent=(-125,125,0,15),colors='m')
#        ax2.clabel(contours_2c, inline=True, fontsize=10)
        display.plot_latitude_slice('reflectivity', lon=cross_lon,lat=cross_lat, vmin=refl_min, vmax=refl_max,
                                      axislabels = ('East-West Distance from Radar (km)', 'Height (km)'), 
                                      cmap = pyart.graph.cm.NWSRef)
        ax2.scatter(lon_km,alt_km,c=colors)
        plt.title(date_time + ' Reflectivity')
#        ax2.plot(lon_km2[cessation_idx],alt1[cessation_idx],'ko')
        ax2.axhline(freezelayerhght,c="blue",linestyle='--')
        ax2.text(ll[0]+1,((freezelayerhght)+0.1),'0C',color='b')
        ax2.axhline(neg10hght,c="red",linestyle=':')
        ax2.text(ll[0]+1,((neg10hght)+0.1),'-10C',color='r')
        ax2.axhline(neg20hght,c="purple",linestyle='-.')
        ax2.text(ll[0]+1,((neg20hght)+0.1),'-20C',color='purple')
        plt.tight_layout(pad=1)
#        sys.exit()
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/atcessation/CrossSections/'+date_time+'.png')
       
#            count = count+1        
#            sys.exit()
    count = count + 1