# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 14:43:32 2016

Grid radar data

@author: kurtispinkney
"""
#==============================================================================
#Import packages
#==============================================================================
import pyart 
import glob
import numpy as np
import numpy.ma as ma
import copy
from time import time
from pyart.core.transforms import antenna_to_cartesian
import sys
from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, 
                            csu_dsd, csu_kdp, csu_misc)
from cut_sails import cut_sails
from skewt import SkewT

#sys.exit()

#Check to make sure that pyart imported correctly (especially if you are concerned about the error)
if hasattr(pyart.graph, 'RadarMapDisplay'):
    print("Py-ART is ready")
else:
    print("ISSUES\n\nMissing\Broken Basemap\n")
    pyart._debug_info()

#==============================================================================
#Gridding flags
#==============================================================================
vert = 1
shy = 0

#==============================================================================
#Read in data
#==============================================================================
files = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/radar_data/2016/*')
#files = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/radar_data/2014/KHTX20120519_000129_V06.gz')

#loop through all radar files
for radar_file in files:
#    sys.exit()
    #extract strings from file names for later operations
    snd_date = radar_file[-22:-14]
    date_time = radar_file[-22:-7]
    
    #read in radar data
#    wsr88d_radar = pyart.io.read(radar_file)
    #apply sails script that Sarah wrote
    o_wsr88d_radar = pyart.io.read(radar_file)
    wsr88d_radar = cut_sails(o_wsr88d_radar)
#    sys.exit()
    #read in sounding data
#    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2015/20150703.txt'
    sndfile = '/Users/kurtispinkney/Desktop/MastersThesis/sounding_data/2016/'+snd_date+'_snd.txt'
    sounding = SkewT.Sounding(sndfile)

    #extract data from the reflectivity and differential reflectivity fields by using deep copy so as to not overwrite original data
    extracted_reflectivity = copy.deepcopy(wsr88d_radar.fields['reflectivity']['data'])
    extracted_zdr = copy.deepcopy(wsr88d_radar.fields['differential_reflectivity']['data'])
    
    #function to extract unmasked data
    def extract_unmasked_data(radar, field, bad=-32768):
        """Simplify getting unmasked radar fields from Py-ART"""
        return radar.fields[field]['data'].filled(fill_value=bad)
     
    #extract unmasked data
    dzN = extract_unmasked_data(wsr88d_radar, 'reflectivity')
    dpN = extract_unmasked_data(wsr88d_radar, 'differential_phase')
    
    #range needs to be supplied as a variable, and it needs to be the same shape as dzN, etc.
    rng2d, az2d = np.meshgrid(wsr88d_radar.range['data'], wsr88d_radar.azimuth['data'])
    
    #calculate kdp, filtered differential phase, and standard deviation of phase
    kdN, fdN, sdN = csu_kdp.calc_kdp_bringi(
        dp=dpN, dz=dzN, rng=rng2d/1000.0, thsd=12, gs=250.0, window=5)
    
    #convert extracted data to linear scale
    linear_reflectivity = ma.power(10.,(extracted_reflectivity/10.))
    linear_zdr = ma.power(10.,(extracted_zdr/10.))

    #define a function to add a field to the radar object (to be able to add new linear-scale data to the radar object)
    def add_field_to_radar_object(field, radar, field_name='linear_reflectivity', units='mm6m-3', 
                              long_name='Linear Reflectivity', standard_name='Linear Reflectivity',
                              dz_field='reflectivity', valid_max=94.5, valid_min=-32.0):
        """
        Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
        make the new field masked the same as reflectivity.
        """
        fill_value = -32768
        masked_field = np.ma.asanyarray(field)
        masked_field.mask = masked_field == fill_value
        if hasattr(radar.fields[dz_field]['data'], 'mask'):
            setattr(masked_field, 'mask', 
                    np.logical_or(masked_field.mask, radar.fields[dz_field]['data'].mask))
            fill_value = radar.fields[dz_field]['_FillValue']
        field_dict = {'data': masked_field,
                      'units': units,
                      'long_name': long_name,
                      'standard_name': standard_name,
                      'valid_max' : valid_max,
                      'valid_min' : valid_min,
                      '_FillValue': fill_value}
        radar.add_field(field_name, field_dict, replace_existing=True)
        return radar
    
    #add the linearized data back to the radar object for gridding
    radar = add_field_to_radar_object(linear_reflectivity, wsr88d_radar, field_name = 'linear_reflectivity', 
                                  units = 'dBZ', long_name = 'Reflectivity in Linear Scale', 
                                  standard_name = 'Linear Reflectivity', dz_field = 'reflectivity',
                                  valid_min = 5.6e-4, valid_max = 2.82e9)

    radar = add_field_to_radar_object(linear_zdr, wsr88d_radar, field_name = 'linear_zdr', units = 'none', 
                                  long_name = 'Differential Reflectivity in Linear Scale',
                                  standard_name = 'Linear ZDR', dz_field = 'reflectivity',
                                  valid_min = 0.16, valid_max = 6.22)
    
    #add kdp data to the radar object for gridding                                
    radar = add_field_to_radar_object(kdN, wsr88d_radar, field_name='KDP', units='deg/km', 
                                       long_name='Specific Differential Phase',
                                       standard_name='Specific Differential Phase', 
                                       dz_field='reflectivity')
    
    radar = add_field_to_radar_object(fdN, wsr88d_radar, field_name='FDP', units='deg', 
                                       long_name='Filtered Differential Phase',
                                       standard_name='Filtered Differential Phase', 
                                       dz_field='reflectivity')
    
    radar = add_field_to_radar_object(sdN, wsr88d_radar, field_name='SDP', units='deg', 
                                       long_name='Standard Deviation of Differential Phase',
                                       standard_name='Standard Deviation of Differential Phase', 
                                       dz_field='reflectivity')
    
    #exclude masked data from gridding operation
    gatefilters = pyart.filters.GateFilter(wsr88d_radar)
    gatefilters.exclude_masked('reflectivity')
    gatefilters.exclude_masked('differential_reflectivity')
#    gatefilters.exclude_masked('linear_reflectivity')
#    gatefilters.exclude_masked('linear_zdr')
    gatefilters.exclude_masked('cross_correlation_ratio')
    gatefilters.exclude_masked('KDP')
    gatefilters.exclude_masked('differential_phase')
    
    #==============================================================================
    #Grid the data
    #==============================================================================
    #for 1km horiz. and 0.5 kmvert grid_shape = 16,251,251; grid_limits = 500, 15500
    #for 2km horiz. and 1.5km vert grid_shape = 2,251,251; grid_limits = 1500, 3000
    #vert grid
    if vert == 1:
        grid = pyart.map.grid_from_radars(wsr88d_radar, #radar data
               grid_shape=(31, 251, 251),   
               grid_limits=((0, 15000.0),(-125000.0, 125000.0), (-125000.0, 125000.0)),  
               fields=['linear_reflectivity', 'linear_zdr', 'reflectivity', 'differential_reflectivity','cross_correlation_ratio','KDP','differential_phase'], 
               gridding_algo="map_gates_to_grid", roi_func = 'constant', constant_roi = 1000.,
               weighting_function='Cressman',gatefilters=gatefilters) 
        print('ive just gridded vert data for '+date_time)
    #shy grid
    if shy == 1:
        grid = pyart.map.grid_from_radars(wsr88d_radar, #radar data
               grid_shape=(2, 126, 126),   
               grid_limits=((1500, 3000.0),(-125000.0, 125000.0), (-125000.0, 125000.0)),   
               fields=['linear_reflectivity', 'linear_zdr', 'reflectivity', 'differential_reflectivity','cross_correlation_ratio','KDP'], 
               gridding_algo="map_gates_to_grid", roi_func = 'constant', constant_roi = 1000.,
               weighting_function='Cressman',gatefilters=gatefilters) 
        print('ive just gridded shy data for '+date_time)
#    sys.exit()
    
    #extract data from the linearized reflectivity and linearized differential reflectivity gridded fields 
    extracted_linear_reflectivity = copy.deepcopy(grid.fields['linear_reflectivity']['data'])
    extracted_linear_zdr = copy.deepcopy(grid.fields['linear_zdr']['data'])
    
    #convert extracted data to log scale
    log_reflectivity = 10.0*ma.log10(extracted_linear_reflectivity)
    log_zdr = 10.0*ma.log10(extracted_linear_zdr)
    
    #define a function to add a field to the grid object
    def add_field_to_grid_object(field, grid, field_name='log_reflectivity', units='dBZ', 
                                  long_name='Reflectivity in dBZ', standard_name='Reflectivity',
                                  dz_field='linear_reflectivity', 
                                  valid_min = -32.0, valid_max = 94.5):
        """
        Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
        make the new field masked the same as reflectivity.
        """
        fill_value = -32768
        #fill_value = 1e20
#        fill_value = -9999.0
        masked_field = np.ma.asanyarray(field)
        masked_field.mask = masked_field == fill_value
        if hasattr(grid.fields[dz_field]['data'], 'mask'):
            setattr(masked_field, 'mask', 
                    np.logical_or(masked_field.mask, grid.fields[dz_field]['data'].mask))
            fill_value = grid.fields[dz_field]['_FillValue']
        field_dict = {'data': masked_field,
                      'units': units,
                      'long_name': long_name,
                      'standard_name': standard_name,
                      'valid_min' : valid_min,
                      'valid_max' : valid_max,
                      '_FillValue': fill_value}
        grid.add_field(field_name, field_dict, replace_existing=True)
        return grid
    
    #add the correctly linear interpolated, gridded data back to the grid object for plotting, etc.
    if vert == 1:
        grid = add_field_to_grid_object(log_reflectivity, grid, field_name = 'reflectivity', units = 'dBZ', long_name = 'Reflectivity in dBZ', standard_name = 'Reflectivity', dz_field='linear_reflectivity', valid_min = -32.0, valid_max = 94.5)
        grid = add_field_to_grid_object(log_zdr, grid, field_name = 'zdr', units = 'dB', long_name = 'Differential Reflectivity in dB', standard_name = 'Differential Reflectivity', dz_field='linear_reflectivity', valid_min = -7.875, valid_max = 7.9375)
    if shy == 1:
        grid = add_field_to_grid_object(log_reflectivity, grid, field_name = 'log_reflectivity', units = 'dBZ', long_name = 'Reflectivity in dBZ', standard_name = 'Reflectivity', dz_field='linear_reflectivity', valid_min = -32.0, valid_max = 94.5)
        grid = add_field_to_grid_object(log_zdr, grid, field_name = 'log_zdr', units = 'dB', long_name = 'Differential Reflectivity in dB', standard_name = 'Differential Reflectivity', dz_field='linear_reflectivity', valid_min = -7.875, valid_max = 7.9375)

    #==============================================================================
    #HID calculations     
    #==============================================================================
    print('im starting hid for '+date_time)
    #define variables needed for HID
    if vert == 1:
        dz = grid.fields['reflectivity']['data']
        dr = grid.fields['zdr']['data']
        kd = grid.fields['KDP']['data']
        rh = grid.fields['cross_correlation_ratio']['data']
    if shy == 1:
        dz = grid.fields['log_reflectivity']['data']
        dr = grid.fields['log_zdr']['data']
        kd = grid.fields['KDP']['data']
        rh = grid.fields['cross_correlation_ratio']['data']

    #function to make sure sounding is monotonic
    def check_sounding_for_montonic(sounding):
        """
        So the sounding interpolation doesn't fail, force the sounding to behave
        monotonically so that z always increases. This eliminates data from
        descending balloons.
        """
        snd_T = sounding.soundingdata['temp']  # In old SkewT, was sounding.data
        snd_z = sounding.soundingdata['hght']  # In old SkewT, was sounding.data
        dummy_z = []
        dummy_T = []
        if not snd_T.mask[0]: #May cause issue for specific soundings
            dummy_z.append(snd_z[0])
            dummy_T.append(snd_T[0])
            for i, height in enumerate(snd_z):
                if i > 0:
                    if snd_z[i] > snd_z[i-1] and not snd_T.mask[i]:
                        dummy_z.append(snd_z[i])
                        dummy_T.append(snd_T[i])
            snd_z = np.array(dummy_z)
            snd_T = np.array(dummy_T)
        return snd_T, snd_z
    
    #match temperature data structure (array) to radar data structure
    def interpolate_data_to_grid(sounding, grid, alt):
        """
        Takes sounding data and interpolates it to gridded radar z-levels. 
        Provide altitude of radar.
        """
        stack = np.empty((len(grid.z['data']),len(grid.x['data']),len(grid.y['data'])))
        i=0
        #Adjust gridded data z for altitude of the radar
        for z in grid.z['data']+alt:
            z_layer = np.zeros((len(grid.x['data']),len(grid.y['data'])))+z
            stack[i,:,:]=z_layer
            i+=1
        grid_z = stack
        grid_T = None
        snd_T, snd_z = check_sounding_for_montonic(sounding)
        shape = np.shape(grid_z)
        grid_z1d = grid_z.ravel()
        grid_T1d = np.interp(grid_z1d, snd_z, snd_T)
        return np.reshape(grid_T1d, shape), grid_z
    
    #set the altitude of the radar for appropriate HID height assignment
    alt = np.array(wsr88d_radar.altitude['data'])
    
    #interpolate sounding data to gridded data
    grid_T, grid_z = interpolate_data_to_grid(sounding,grid,alt)
    
    ##Default weights
    set_weights = {'DZ': 1.5, 'DR': 0.8, 'KD': 1.0, 'RH': 0.8, 'LD': 0.0, 'T': 0.4}
    hid_identifier = 'def'

#    set_weights = {'DZ': 1.5, 'DR': 0.6, 'KD': 1.0, 'RH': 0.8, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'def2'  

#    set_weights = {'DZ': 1.5, 'DR': 0.8, 'KD': 1.0, 'RH': 1.0, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'def3'       

    ##1st adjustment from default
#    set_weights = {'DZ': 1.5, 'DR': 1.0, 'KD': 0.75, 'RH': 0.8, 'LD': 0.0, 'T': 0.4}    
#    hid_identifier = 'new'
        
    ##2nd adjustment from default
#    set_weights = {'DZ': 1.5, 'DR': 0.8, 'KD': 0.7, 'RH': 0.6, 'LD': 0.0, 'T': 0.4}    
#    hid_identifier = 'new2'
            
    ##Adjustment of RHO based on 1st iteration
#    set_weights = {'DZ': 1.5, 'DR': 1.0, 'KD': 0.75, 'RH': 0.6, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'rho'
        
    ##Adjustment of ZDR
#    set_weights = {'DZ': 1.5, 'DR': 1.2, 'KD': 0.8, 'RH': 0.8, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'zdr'
        
    ##Adjustment of ZDR, RHO
#    set_weights = {'DZ': 1.5, 'DR': 1.2, 'KD': 0.8, 'RH': 0.6, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'zdrrho'
        
    ##2nd Adjustment of ZDR, RHO
#    set_weights = {'DZ': 1.5, 'DR': 1.2, 'KD': 0.8, 'RH': 0.4, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'zdrrho2'
        
    ##Adjustment of ZDR, RHO, KDP
#    set_weights = {'DZ': 1.5, 'DR': 1.2, 'KD': 1.0, 'RH': 0.4, 'LD': 0.0, 'T': 0.4}
#    hid_identifier = 'zrk'
    
    #determine HID scores
    scores = csu_fhc.csu_fhc_summer(dz=dz, zdr=dr, rho=rh, kdp=kd, use_temp=True, band='S',
                                    T=grid_T, weights = set_weights)
    fh = np.argmax(scores, axis=0) + 1
#    sys.exit()
    #add HID results to grid
    grid = add_field_to_grid_object(fh, grid, field_name = 'FH', units = 'unitless', long_name = 'Hydrometeor ID', standard_name = 'Hydrometeor ID')

    #==============================================================================
    #Mass calculations   
    #==============================================================================
    #determine mass of water and ice 
    mw, mi = csu_liquid_ice_mass.calc_liquid_ice_mass(dz, dr, grid_z/1000.0, T=grid_T, method=None, fit_a=None, fit_b=None)
     
    #add mass calculation results to grid
    grid = add_field_to_grid_object(mw, grid, field_name = 'Mw', units = 'g/m^3', long_name = 'Liquid Water Mass', standard_name = 'Liquid Water Mass')
    grid = add_field_to_grid_object(mi, grid, field_name = 'Mi', units = 'g/m^3', long_name = 'Ice Water Mass', standard_name = 'Ice Water Mass')
    #==============================================================================
    #Export grid 
    #==============================================================================
    if vert == 1:
        #grid filename
#        grid_filename = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/postpriorcessation/vertgrid/test'+date_time+'_'+hid_identifier+'.nc'
        grid_filename = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/postpriorcessation/2016/'+date_time+'_'+hid_identifier+'.nc'
    
        #write out the netcdf files with the gridded data using the filename and gridded data
        pyart.io.write_grid(grid_filename, grid)
        print('ive printed out '+date_time)
    
    if shy == 1:
        #grid filename
        grid_filename = '/Users/kurtispinkney/Desktop/MastersThesis/gridded_data/postpriorcessation/shygrid/'+date_time+'_'+hid_identifier+'.nc'
    
        #write out the netcdf files with the gridded data using the filename and gridded data
        pyart.io.write_grid(grid_filename, grid)
        print('ive printed out '+date_time)
#    sys.exit()