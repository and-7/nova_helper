#!/usr/bin/env python

## importing libraries

from astropy.table import Table
import numpy as np
from nova_helper import read_data,return_obs_data,combine_obs_data

## using nova_helper to create a smoothed light curve in the V band for the nova V1723 Sco

obs_V_raw,JD_V_raw,mag_V_raw,uncer_V_raw = read_data('ExampleData_V1723_Sco_aavso_2024-06-17.csv','V')
obs_data_V = return_obs_data(obs_V_raw,JD_V_raw,mag_V_raw,uncer_V_raw)
obs_V,JD_V,mag_V,uncer_V,offset_V = combine_obs_data(obs_data_V)

## plotting using matplotlib:
import matplotlib.pyplot as plt

## raw curve
plt.scatter(JD_V_raw,mag_V_raw,s=0.5,c='k')
plt.xlim([np.amin(JD_V)-40,np.amax(JD_V)+40])
plt.xlabel('JD')
plt.ylim([np.amax(mag_V)+2,np.amin(mag_V)-2])
plt.ylabel('Magnitude')
plt.title('V1723 Sco Light Curve (No Offsets, Band=V)')
plt.show()

## smoothed curve
plt.scatter(JD_V,mag_V,s=0.5,c='b')
plt.xlim([np.amin(JD_V)-40,np.amax(JD_V)+40])
plt.xlabel('JD')
plt.ylim([np.amax(mag_V)+2,np.amin(mag_V)-2])
plt.ylabel('Magnitude')
plt.title('V1723 Sco Light Curve (Offsets, Band=V)')
plt.show()

## saving the smoothed data and log of offset information to csv files using Astropy.table

Table([JD_V,mag_V,uncer_V,obs_V],names=['JD','Magnitude','Uncertainty','Observer']).write('Example_V1723_Sco_V_data.csv',format='csv',overwrite=True)
Table([offset_V[0],offset_V[1],offset_V[2]],names=['Observer Code','Number of Observations','Offset']).write('Example_V1723_Sco_V_offset_log.csv',format='csv',overwrite=True)
