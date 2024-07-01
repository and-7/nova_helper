##############################################

from astropy.table import Table
import numpy as np

##############################################

def read_data(file,bandpass):    
   ## reading in data as a table and extracting necessary data
    mytable = Table.read(file, format='csv')
    JD0,magnitude0,uncertainty0,band0,observer0 = np.array(mytable['JD']),list([mytable['Magnitude']])[0],np.array([mytable['Uncertainty']])[0],np.array([mytable['Band']])[0],np.array([mytable['Observer Code']])[0]
    ## some of the magnitudes have a < sign
    ## these are upper limits rather than actual brightness measurements
    ## so we need to mask them out of the lgiht curve data
    for i in range(len(magnitude0)):
        if magnitude0[i][0] == '<':
            magnitude0[i] = None
    magnitude0 = np.array(magnitude0)
    JD,magnitude,uncertainty,band,observer = JD0[magnitude0 != 'None'],magnitude0[magnitude0 != 'None'],uncertainty0[magnitude0 != 'None'],band0[magnitude0 != 'None'],observer0[magnitude0 != 'None']
    magnitude = np.array(magnitude,copy=True,dtype=float)
    ## now we mask out the requested band and return the JD, magnitude, uncertainty, observer
    return observer[band==bandpass],JD[band==bandpass],magnitude[band==bandpass],uncertainty[band==bandpass]

#############################################

def remove_points(observer,JD,magnitude,uncertainty,target_JD):
    ## finding all of the points within 0.001 day of points in the target set
    for i in range(len(target_JD)):
        for j in range(len(JD)):
            if np.abs(target_JD[i]-JD[j]) <= 0.001:
                JD[j] = 0
    ## now we mask out the requested points and return the JD, magnitude, uncertainty, observer as Numpy arrays
    observer_mask,JD_mask,magnitude_mask,uncertainty_mask = observer[JD != 0],JD[JD != 0],magnitude[JD != 0],uncertainty[JD != 0] 
    return observer_mask,JD_mask,magnitude_mask,uncertainty_mask

##############################################


## returning observation data as a nested list:
## [ [[observer1],[JD1],[mag1],[uncer1]] , [[observer2],[JD2],mag2],[uncer2]] , ...]
## sorted in order of the observer code with most to least amount of observations
## if an observer has less than 10 observations, then their data is not included

def return_obs_data(observer,JD,magnitude,uncertainty):
    ## get a list of all the unique observation codes in the V band
    obs_codes_unsorted = []
    for i in observer:
        if i not in obs_codes_unsorted:
            obs_codes_unsorted.append(i)
    ## for each observer, count the associated number of observations
    ## put the observation codes in order of most to least observations
    ## with the minumum number of observations per code being 10
    num_obs = {}
    for j in obs_codes_unsorted:
        if len(JD[observer==j]) >= 10:
            num_obs[j] = len(JD[observer==j])
    obs_codes = dict(sorted(num_obs.items(), key=lambda item: item[1]))
    obs_codes = list(reversed(obs_codes))
    ## assemble a list containing observer codes and all respective observation data
    obs_data = []
    for k in obs_codes:
        temp = [[k for i in range(len(JD[observer==k]))]]
        temp.append(list(JD[observer==k]))
        temp.append(list(magnitude[observer==k]))
        temp.append(list(uncertainty[observer==k]))
        obs_data.append(temp)
    ## return the sorted observation data list
    return obs_data
        
##############################################

## function that combines two sets of data
## into one sorted set of observer, JD, mag, and mag uncertainty

def combined_sorted_data(obs1,JD1,mag1,uncer1,obs2,JD2,mag2,uncer2):
    ## append the new list to the old list
    obs,JD,mag,uncer = obs1.copy(),JD1.copy(),mag1.copy(),uncer1.copy()
    for i in range(len(JD2)):
        obs.append(obs2[i])
        JD.append(JD2[i])
        mag.append(mag2[i])
        uncer.append(uncer2[i])
    ## zip the lists together and sort based on JD
    zipped = zip(JD,obs,mag,uncer)
    zipped = sorted(zipped)
    ## unzip and return the sorted lists
    JD,obs,mag,uncer = zip(*zipped)
    return list(obs),list(JD),list(mag),list(uncer)

##############################################

## we need to find the optimal offset of the second set of data
## to make the combined curve as smooth as possible

## 1/nu parameter calculated for the magnitude
## (standard deviation squared / mean difference between consecutive values squared)

def curve_smoothness(m):
    ## mean value of data set
    mbar = np.mean(m)
    ## number of measurements
    N = len(m)
    ## calculating squared std
    sigma_sq = 0
    for i in range(N):
        sigma_sq += (m[i]-mbar)**2 / (N-1)
    ## calculating mean difference between consecutive values
    delta_sq = 0
    for j in range(N-1):
        delta_sq += (m[j+1]-m[j])**2 / (N-1)
    ## returning the 1/nu value (the ratio of sigma_sq and delta_sq)
    return sigma_sq / delta_sq

##############################################

## defining a function to optimize the offset between two data sets
## args are the reference JD and mag, the new JD and mag, and an array of potential offsets
## default array of potential offsets is np.linspace(-2,2,1500)

## this function runs the smoothness function for every offset in the array
## then it sees which offset produced the optimal output from the smoothness function
## then that offset is applied to combine the data sets into one optimized data set

def optimize_offset(obs_ref,JD_ref,mag_ref,uncer_ref,obs_new,JD_new,mag_new,uncer_new,offsets=np.linspace(-2,2,1500)):
    smoothness = np.zeros(len(offsets))
    ## looping through every entry in the offsets
    for i in range(len(offsets)):
        obs_comb,JD_comb,mag_comb,uncer_comb = combined_sorted_data(obs_ref,JD_ref,mag_ref,uncer_ref,obs_new,JD_new,mag_new+offsets[i],uncer_new)
        smoothness[i] = curve_smoothness(mag_comb)
    ## finding the offset which produced the smoothest curve (minimizing 1/nu)
    optimal_offset = offsets[np.argmax(smoothness)]
    ## combining the reference data and the offset data into one  set
    obs,JD,mag,uncer = combined_sorted_data(obs_ref,JD_ref,mag_ref,uncer_ref,obs_new,JD_new,mag_new+optimal_offset,uncer_new)
    return obs,JD,mag,uncer,optimal_offset

##############################################

## defining a function that takes many observer sets in a band and optimizes them into one light curve
## assumes that the observation data is sorted by observer in the following nested list format:
## [ [[observer1],[JD1],[mag1],[uncer1] ] , [[observer2],[JD2],mag2],[uncer2]] , ...]
## sorted in descending order from most amount of observations to least amount of observations

def combine_obs_data(obs_data):
    ## reference JD and mag is the first in the list (most observation sets)
    obs_ref,JD_ref,mag_ref,uncer_ref = obs_data[0][0],obs_data[0][1],obs_data[0][2],obs_data[0][3]
    ## creating a log of the observer code, the number of observations, and the chosen offset
    offset_names = [obs_data[0][0][0]]
    offset_len = [len(obs_ref)]
    offset_values = [0.]
    ## looping through all observers
    for i in range(len(obs_data)-1):
        ## new JD and mag is the next in the list
        obs_new,JD_new,mag_new,uncer_new = obs_data[i+1][0],obs_data[i+1][1],obs_data[i+1][2],obs_data[i+1][3]
        ## optimizing offsets
        obs_temp,JD_temp,mag_temp,uncer_temp,offset_temp = optimize_offset(obs_ref,JD_ref,mag_ref,uncer_ref,obs_new,JD_new,mag_new,uncer_new)
        offset_names.append(obs_new[0])
        offset_len.append(len(obs_new))
        offset_values.append(offset_temp)
        ## new reference JD and mag is the previous combined data
        obs_ref,JD_ref,mag_ref,uncer_ref = obs_temp,JD_temp,mag_temp,uncer_temp
    ## returning the combined data from all observers
    ## and the offset log data
    return obs_ref,JD_ref,mag_ref,uncer_ref,(offset_names,offset_len,offset_values)

##############################################
