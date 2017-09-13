# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a script to attempt to qc the flash data. I want to isolate flashes with at least 10 sources.

Written By: Austin Vacek, July 2016
Editted By: Kurtis Pinkney, November 2016
"""

#import packages
import glob

#==============================================================================
# Multiple files at once
#==============================================================================
#read in flash data
files = glob.glob('/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/raw_flash_data/2016/*.dat')

for file in files:
    date = file[-10:-4]
    print(date)
    outfile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_flash_data/2016/flash_'+date+'qc.dat'
    with open(file,'r+') as infile:
        with open(outfile, 'w') as outfile:
            data = infile.readlines()
            goodrows = [row for row in data if not float(row.split()[4]) < 10]
            outfile.writelines(goodrows)

#==============================================================================
# Single files
#==============================================================================
#file = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/raw_flash_data/2flash_120611.dat'
#outfile = '/Users/kurtispinkney/Desktop/MastersThesis/lxn_data/qc_flash_data/2flash_120611qc.dat'
#with open(file,'r+') as infile:
#        with open(outfile, 'w') as outfile:
#            data = infile.readlines()
#            goodrows = [row for row in data if not float(row.split()[4]) < 10]
#            outfile.writelines(goodrows)

#==============================================================================
# xlma files
#==============================================================================
