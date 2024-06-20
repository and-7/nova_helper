# nova-helper
## A set of Python functions which read in nova light curve data, optimize zero-point magnitude offsets between observers, and produce a smoothed curve. 
## This program currently supports comma-separated data. Headers for the input data table must include 'JD', 'Magnitude', 'Uncertainty', 'Band', and 'Observer'.
## Defining helper functions for analyzing nova light curves:

-read_data(file,bandpass) reads in the csv data file and selects a specific bandpass. It returns Numpy arrays observer, JD, magnitude, uncertainty.

-return_obs_data(observer,JD,magnitude,uncertainty) creates a nested list with the observation data, sorted by observer code, in descending number of observations. Any observer code with less than 10 observations is not included. The returned list is of the format [ [[observer1],[JD1],[mag1],[uncer1]] , [[observer2],[JD2],mag2],[uncer2]] , ...]. 

-combined_sorted_data(obs1,JD1,mag1,uncer1,obs2,JD2,mag2,uncer2) combines two sets of data consisting of observer, JD, magnitude, and uncertainty, and sorts in order of ascending JD. It returns lists obs, JD, mag, uncer. 

-curve_smoothness(m) returns the 1/nu parameter (standard deviation squared divided by mean difference between consecutive values squared) for a given data set m. 

-optimize_offset(obs_ref,JD_ref,mag_ref,uncer_ref,obs_new,JD_new,mag_new,uncer_new,offsets=np.linspace(-2,2,1500)) combines two data sets into one optimally smooth curve by looping through a range of possible magnitude offsets for the new data, calculating the 1/nu parameter for each. The offset with the largest 1/nu parameter is selected and the offset is applied. The default set of possible offsets is 1500 linearly spaced values ranging from -2 to 2. The returned values are lists obs, JD, mag, uncer and float optimal_offset.

-combine_obs_data(obs_data) combines all data in a sorted nested list into one optimally smooth curve. The input obs_data must be a nested list in the same format returned by return_obs_data. The first observer code in the nested list (the one with the most observations) is taken as the reference with respect to magnitude offsets. Returned values are lists obs_ref, JD_ref, mag_ref, uncer_ref and tuple (offset_names, offset_len, offset_values). This includes smoothed observer, JD, magnitude, and uncertainty data; offset_names, offset_len, and offset_values are lists containing information for each observer code. Offset_names contains the codes, offset_len contains the respective number of observations, and offset_values contains the respective magnitude offset.

### See ExampleScript_V1723_Sco_V.txt for an example of smoothing a light curve using this program.
