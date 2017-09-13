#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 12:57:58 2017

Script to do analysis prior/post cessation

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
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import mpl_toolkits.axisartist as AA
import sys
import datetime
from scipy.spatial import ConvexHull
import scipy
from skewt import SkewT

#==============================================================================
#Calculation flags
#==============================================================================
xlma = 0
ppi = 0
vppr = 0
box = 1
maxzh = 0
icevol = 0
icevolhght = 1

#==============================================================================
#Read in data 
#==============================================================================
vertfiles = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/postpriorcessation/2016/*')
#boxfiles = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/*')

#initialize the box_file_test variable
box_file_test = 'nah fam' 
#count = 0
#boxfile = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/2015/150612.csv'
##        boxfile = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/2015/'+flash_date+'.csv'
#boxpoints = np.genfromtxt(boxfile,skip_header=1,delimiter=',')

#read in data
for vertfile in vertfiles:
    #assign strings to variables for file in/out processes
    flash_date = vertfile[-20:-14]
    file_time = vertfile[-13:-7]
    date_time = vertfile[-22:-7]
    snd_date = vertfile[-22:-14]

    #check to see if the same date is still being used    
    if flash_date != box_file_test:
#        boxfile = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/2015/150612.csv'
        boxfile = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/bounding_boxes/2016/'+flash_date+'.csv'
        boxpoints = np.genfromtxt(boxfile,skip_header=1,delimiter=',')
        count = 0 #to loop through each line in boxpoints array later
        
        #ice volumes with height
        graupbyhght = np.zeros([31,np.size(boxpoints[:,0])])
        icebyhght = np.zeros([31,np.size(boxpoints[:,0])])
        
        #initialize lists to hold values at each temperature height
        freezeZhmax = []
        neg10Zhmax = []
        neg20Zhmax = []
        
        freezeZhavg = []
        neg10Zhavg = []
        neg20Zhavg = []
        
        freezeZdrmax = []
        neg10Zdrmax = []
        neg20Zdrmax = []
        
        freezeZdravg = []
        neg10Zdravg = []
        neg20Zdravg = []
        
        labels = []
        
        graupelvol = []
        nonprecipicevol = []
        
        Zhlapse = []
    else:
        print('Still using '+flash_date+'s boxfile')
    
    #update the name of the date being used for bounding box purposes
    box_file_test = flash_date
    
    #xaxis labels for plotting functions later
    labels.append(file_time)
    
    #read in higher vertical resolution gridded data
    vertgrid = pyart.io.read_grid(vertfile,exclude_fields=None)

    #assign reflectivity data to a variable
    vertgrid.fields['reflectivity']['data'] = ma.masked_less_equal(vertgrid.fields['reflectivity']['data'],10)
    vertrefl = vertgrid.fields['reflectivity']['data']
    
    #assign Zdr data to a variable
    vertzdr = vertgrid.fields['zdr']['data']
    
    #assign Kdp data to a variable
    vertkdp = vertgrid.fields['KDP']['data']
    
    #assign HID data to a variable
    vertHID = vertgrid.fields['FH']['data']
    vertHID = ma.getdata(vertHID)
    vertHID = np.round(vertHID)
    
    #read in distances from higher resolution grid center and convert to km
    vert_x_array_km = vertgrid.x['data']/1000
    vert_y_array_km = vertgrid.y['data']/1000
    vert_z_array_km = vertgrid.z['data']/1000
    
    #assign grid lat/lon data to a variable
    vertlat = vertgrid.point_latitude['data']
    vertlon = vertgrid.point_longitude['data']

    #==============================================================================
    #Sounding data   
    #==============================================================================
    #read in sounding data
#    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2015/20150619.txt'
#    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2015/'+snd_date+'.txt'
    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2016/'+snd_date+'_snd.txt'
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
    #XLMA lightning data  
    #==============================================================================
    if xlma == 1:
        #read in lightning data
        lxnfile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_xlma_data/xlma_'+flash_date+'.txt'
        data = np.genfromtxt(lxnfile,dtype=None,names=['time','lat','lon','alt'], usecols = (0,1,2,3))
   
        #assign data columns to variables
        time = data['time']
        lat = data['lat']
        lon = data['lon']
        alt = data['alt']
        
        #convert lat/lon from ° to km and altitude to km
        lat_km = (lat-vertgrid.origin_latitude['data'])*110.926
        lon_km = (lon-vertgrid.origin_longitude['data'])*92.15
        alt_km = (alt)/1000.

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
    #VPPR    
    #==============================================================================
    if vppr == 1:
        #initialize array to hold reflectivity values within bounding box at different temperature heights
        freezeZh = np.zeros(np.size(xsort))
        neg20Zh = np.zeros(np.size(xsort))
        
        #create arrays with Zh values within bounding box
        for i in range(np.size(xsort)):
            freezeZh[i]=vertrefl[freezegrididx][ysort[i]][xsort[i]]
            neg20Zh[i]=vertrefl[neg20grididx][ysort[i]][xsort[i]]
        
        #identify index of max reflectivity within bounding box at 0°C
        freezeZhmaxidx = np.nanargmax(freezeZh)
        #identify index of max reflectivity within bounding box at -20°C
        neg20Zhmaxidx = np.nanargmax(neg20Zh)
        
        #identify max reflectivity value at freezing level
        freezeZhmax = freezeZh[freezeZhmaxidx]
        
        #identify reflectivity value at -20°C 
        neg20Zhmax = neg20Zh[neg20Zhmaxidx]
        
        #append lapse rate data to list
        Zhlapse.append((neg20Zhmax - freezeZhmax)/3)#3 km is typical distance between freeze height and -20°C height
        height_change = vert_z_array_km[neg20grididx]-vert_z_array_km[freezegrididx]
    #==============================================================================
    #Plotting VPPR results
    #==============================================================================
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance    
    
        #generate plot
        ax1 = plt.subplot(111)
        ax1.plot(Zhlapse)
        ax1.legend()
        
        #labels
        plt.title(snd_date+' Reflectivity Lapse Rate During Cessation Period (0 to -20°C)')
        plt.xlabel('Radar Volume Scan Time (UTC)')
        plt.ylabel('Reflectivity Lapse Rate (dBZ)')
        plt.tight_layout(pad=1)
        
        #format x-axis to note time of each plotted point
        pos_list = np.arange(len(labels))
        ax1.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax1.xaxis.set_major_formatter(ticker.FixedFormatter((labels)))
#        ax1.xaxis.grid(True,'major',linewidth=2)
        
        #format y-axis
#        ax1.set_ylim(25,50)
#        ax1.yaxis.set_major_locator(MultipleLocator(5))
#        ax1.yaxis.set_minor_locator(MultipleLocator(1))
#        ax1.yaxis.grid(True,'major',linewidth=2)
#        ax1.yaxis.grid(True,'minor')
        #save figure
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/priorpostcessation/VPPR/2012/'+flash_date+'vppr.png')
        plt.close(fig)
        plt.clf()
    #==============================================================================
    #Plotting PPI images
    #==============================================================================
    if ppi == 1:
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance    
        
        #composite Zh testing
        Zhmax = np.zeros((np.size(vert_x_array_km),np.size(vert_y_array_km)))
        for aa in range(np.size(vert_x_array_km)):
            for bb in range(np.size(vert_x_array_km)):
                Zhmax[aa][bb] = np.amax(vertrefl[:,aa,bb])
                
        #generate plots
        ax = plt.subplot(111)
        a=ax.imshow(Zhmax,
                  origin='lower', vmin=10, vmax=65,extent=(-125, 125, -125, 125)
                  ,cmap = pyart.graph.cm.NWSRef) 
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
        
        #plot bounding box
        if box == 1:
            ax.add_patch(patches.Rectangle((ll[0], ll[1]), ur[0]-ll[0], ur[1]-ll[1], fill=False,color='red',linewidth=2))

        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/priorpostcessation/PPI/2012/'+date_time+'_priorpost.png')
        plt.close(fig)
        plt.clf()
    
    #==============================================================================
    #Max Zh calculations within bounding box              
    #==============================================================================
    if maxzh == 1:
        #determine grid points closest to relevant temp heights
        freeze_zidx_init = (np.abs(vert_z_array_km-freezelayerhght)).argmin()
        neg10_zidx_init = (np.abs(vert_z_array_km-neg10hght)).argmin()
        neg20_zidx_init = (np.abs(vert_z_array_km-neg20hght)).argmin()
        
        #initialize array to hold reflectivity values within bounding box at different temperature heights
        freezeZh = np.zeros(np.size(xsort))
        neg10Zh = np.zeros(np.size(xsort))
        neg20Zh = np.zeros(np.size(xsort))

        #create arrays with Zh values within bounding box
        for i in range(np.size(xsort)):
            freezeZh[i]=vertrefl[freeze_zidx_init][ysort[i]][xsort[i]]
            neg10Zh[i]=vertrefl[neg10_zidx_init][ysort[i]][xsort[i]]
            neg20Zh[i]=vertrefl[neg20_zidx_init][ysort[i]][xsort[i]]

        #determine max Zh values at each temperature height
        freezeZhmax.append(np.nanmax(freezeZh))
        neg10Zhmax.append(np.nanmax(neg10Zh))
        neg20Zhmax.append(np.nanmax(neg20Zh))
        
        #plot data
        if box_file_test == flash_date:
            plotfreezeZhmax = freezeZhmax
            plotneg10Zhmax = neg10Zhmax
            plotneg20Zhmax = neg20Zhmax
        
        
        #determine avg Zh values at each temperature height
        freezeZhavg.append(np.nanmean(freezeZh))
        neg10Zhavg.append(np.nanmean(neg10Zh))
        neg20Zhavg.append(np.nanmean(neg20Zh))
    #==============================================================================
    #Plotting for max Zh at relevant temperature heights 
    #==============================================================================
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance    
    
        #generate plot
        ax1 = plt.subplot(111)
        ax1.plot(plotfreezeZhmax, label='0°C')
        ax1.plot(plotneg10Zhmax, label='-10°C')
        ax1.plot(plotneg20Zhmax, label='-20°C')
        ax1.legend()
        
        #labels
        plt.title('Maximum Reflectivity on '+snd_date+' During Cessation Period')
        plt.xlabel('Radar Volume Scan Time (UTC)')
        plt.ylabel('Reflectivity (dBZ)')
        plt.tight_layout(pad=1)
        
        #format x-axis to note time of each plotted point
        pos_list = np.arange(len(labels))
        ax1.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax1.xaxis.set_major_formatter(ticker.FixedFormatter((labels)))
        ax1.xaxis.grid(True,'major',linewidth=2)
        
        #format y-axis
        ax1.yaxis.set_major_locator(MultipleLocator(5))
        ax1.yaxis.set_minor_locator(MultipleLocator(1))
        ax1.yaxis.grid(True,'major',linewidth=2)
        ax1.yaxis.grid(True,'minor')
        
        #save figure
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/priorpostcessation/MaxZhTemps/2016/'+flash_date+'Zh.png')
        plt.close(fig)
        plt.clf()
    
    #==============================================================================
    #Ice volume calculations    
    #==============================================================================
    #sum graupel points in mixed phase region in grid boxes identified as convective
    if icevol == 1:
        graupcount = 0
        icecount = 0
        
        #determine closest gridded value to bottom of mixed phase region
        gridmpbotidx = np.abs(freezelayerhght-vert_z_array_km).argmin()
    
        #determine closest gridded value to top of mixed phase region
        gridmptopidx = np.abs(neg20hght-vert_z_array_km).argmin()
        
        #determine grid point associated with negative 10 height
        gridneg10idx = np.abs(neg10hght-vert_z_array_km).argmin()
    
        #identify HID layer corresponding to mixed phase region
        convmixedphase = vertHID[gridmpbotidx:gridmptopidx]

#        #determine ice volumes within mixed phase region   
#        for k in range(gridmpbotidx,gridmptopidx):
#            for i in range(np.size(xsort)):
#                if vertHID[k][ysort[i]][xsort[i]] == 7 or vertHID[k][ysort[i]][xsort[i]] == 8:
#                    graupcount = graupcount + 1
#                if vertHID[k][ysort[i]][xsort[i]] == 3 or vertHID[k][ysort[i]][xsort[i]] == 6 or vertHID[k][ysort[i]][xsort[i]] == 4:
#                        icecount = icecount + 1
        
        #determine total amount of graupel within bounding box at -10°C height
        for i in range(np.size(xsort)):
            if vertHID[gridneg10idx][ysort[i]][xsort[i]] == 7 or vertHID[gridneg10idx][ysort[i]][xsort[i]] == 8:
                graupcount = graupcount + 1
            
        #determine graupel echo volume
        graupelvol = np.append(graupelvol, graupcount * .5)#.5 is grid box volume in km^3
        nonprecipicevol = np.append(nonprecipicevol, icecount * .5)#.5 is grid box volume in km^3
    
    #==============================================================================
    #Plotting for ice volume with time    
    #==============================================================================
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance    
    
        #generate plot
        ax1 = plt.subplot(111)
        ax1.plot(graupelvol, label='Big Ice')
        ax1.plot(nonprecipicevol, label='Small Ice')
        ax1.legend()
        
        #labels
        plt.title('Big Ice Volumes at -10°C')
        plt.xlabel('Radar Volume Scan Time (UTC)')
        plt.ylabel('Ice Volume (km^3)')
        plt.tight_layout(pad=1)
        
        #format x-axis to note time of each plotted point
        pos_list = np.arange(len(labels))
        ax1.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax1.xaxis.set_major_formatter(ticker.FixedFormatter((labels)))
#        ax1.xaxis.grid(True,'major',linewidth=2)
#        
#        #format y-axis
#        ax1.set_ylim(25,50)
#        ax1.yaxis.set_major_locator(MultipleLocator(5))
#        ax1.yaxis.set_minor_locator(MultipleLocator(1))
#        ax1.yaxis.grid(True,'major',linewidth=2)
#        ax1.yaxis.grid(True,'minor')
        
        #save figure
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/priorpostcessation/IceVol/2012/'+flash_date+'.png')
        plt.close(fig)
        plt.clf()
    
    #determine graupel volume with height
    if icevolhght == 1:
        #determine closest gridded value to bottom of mixed phase region
        gridmpbotidx = np.abs(freezelayerhght-vert_z_array_km).argmin()
    
        #determine closest gridded value to top of mixed phase region
        gridmptopidx = np.abs(neg20hght-vert_z_array_km).argmin()
        
        #determine ice volumes with height
        for z in range(np.size(vert_z_array_km)):
            graupcount = 0
            icecount = 0
            for y in range(np.size(xsort)):
                #if HID is equal to low-density or high-density graupel
                if vertHID[z][ysort[y]][xsort[y]] == 7 or vertHID[z][ysort[y]][xsort[y]] == 8:
                    graupcount = graupcount + 1
                #if HID is equal to ice crystals, aggregates, or vertical ice
                if vertHID[z][ysort[y]][xsort[y]] == 3 or vertHID[z][ysort[y]][xsort[y]] == 4 or vertHID[z][ysort[y]][xsort[y]] == 6:
#                if vertHID[z][ysort[y]][xsort[y]] == 4:
#                if vertHID[z][ysort[y]][xsort[y]] == 3 or vertHID[z][ysort[y]][xsort[y]] == 6:
                    icecount = icecount + 1
            
            #determine graupel volume by multiplying by volume of grid box (km^3)
            graupvol = graupcount * .5
            icevol = icecount * .5
            #fill in graup by height array
            graupbyhght[z][count] = graupvol
            icebyhght[z][count] = icevol
            
        #generate figure for plotting
        fig_dimensions = [20.0, 10.0] #Sets plot size for the notebook and figures created (width, height in inches)
        plt.rcParams['figure.figsize'] = fig_dimensions #Establishes plot space in the notebook
        plt.rcParams.update({'font.size': 20}) #Changes the plot font size
        fig = plt.figure(figsize=fig_dimensions) #Creates the matplotlib figure instance    
    
        #generate plot
        ax1 = plt.subplot(111)
        ax1.plot(graupbyhght[:,count], vert_z_array_km, label='Big Ice')
        ax1.plot(icebyhght[:,count],vert_z_array_km, label='Small Ice')
#        ax1.set_xscale('log')
#        ax1.plot(nonprecipicevol, label='nonprecip ice')
        ax1.legend()
        
        #labels
        plt.title(date_time+' Ice Volume With Height')
        plt.xlabel('Ice Volume 'r'($km^3$)')
        plt.ylabel('Height km')
        plt.tight_layout(pad=1)
        #plot lines to represent significant temperature heights
        plt.axhline(freezelayerhght,color='blue')
        ax1.axhline(freezelayerhght,c="blue",linestyle='--')
        ax1.text(0,((freezelayerhght)+0.1),'0C',color='b')
        ax1.axhline(neg10hght,c="red",linestyle=':')
        ax1.text(0,((neg10hght)+0.1),'-10C',color='r')
        ax1.axhline(neg20hght,c="purple",linestyle='-.')
        ax1.text(0,((neg20hght)+0.1),'-20C',color='purple')
        ax1.axhline(neg40hght,c="magenta",linestyle='-')
        ax1.text(0,((neg40hght)+0.1),'-40C',color='magenta')
        
        #format x-axis to note time of each plotted point
#        pos_list = np.arange(len(labels))
#        ax1.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
#        ax1.xaxis.set_major_formatter(ticker.FixedFormatter((labels)))
#        ax1.xaxis.grid(True,'major',linewidth=2)
#        
#        #format y-axis
#        ax1.set_ylim(25,50)
#        ax1.yaxis.set_major_locator(MultipleLocator(5))
#        ax1.yaxis.set_minor_locator(MultipleLocator(1))
#        ax1.yaxis.grid(True,'major',linewidth=2)
#        ax1.yaxis.grid(True,'minor')
        #save figure
        plt.savefig('/Users/kurtispinkney/Desktop/MastersThesis/Images/priorpostcessation/IceVol/Height/2016/'+date_time+'.png')
        plt.close(fig)
        plt.clf()
    count = count + 1
