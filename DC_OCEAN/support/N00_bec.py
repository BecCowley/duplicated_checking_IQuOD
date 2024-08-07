#!/usr/bin/env python3
"""
This script used for pre-processing profile data (only supports WOD18 netCDF format files).
During the process:
--Numerical metadata are retained.
--String metadata are converted into numerical values by using the ASCII code to convert each letter (including spaces) and then add the ASCII code of each letter to obtain final values.
The results are stored in npz format file.

***For data that is not in WOD18 netCDF format, it needs to be rewritten as WOD18 format firstly.
The variables included in the WOD18 netCDF files and their order are (keep blank if there are no variables):
Access_no
Biology_Accno
center_source
COMS_id
country
crs
dataset
date
Displacemnt_Vol_001
GMT_sample_start_time
lat
Lge_removed
lon
Orig_Stat_Num
originators_cruise_identifier
Preservation
Project
Sampling_duration
Taxon_count_002
Taxon_count_003
temp_QC_COMS
Temperature
Temperature_BOT_cor
Temperature_sigfigs
Temperature_WODflag
Temperature_WODprofileflag
time
WOD_cruise_identifier
wod_unique_cast
WODf
WODfd
WODfp
z
z_BOT_cor
z_sigfig
z_WODflag

Usage:
    Run this script and follow the prompt to enter the directory path containing netCDF files.
"""
import os
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from datetime import datetime, timedelta
import warnings
import sys
import pandas as pd
import xarray as xr

def validate_path(input_path):
    # Normalize the path
    normalized_path = os.path.normpath(input_path)

    # Check if the path exists and is a directory
    if not os.path.exists(normalized_path) or not os.path.isdir(normalized_path):
        return False

    return True


def read_netCDF_formatted_DNA_series(filepath):

    warnings.filterwarnings("ignore")

    # Path where the netCDF files are stored
    filepath = os.path.normpath(os.path.abspath(filepath))
    print(filepath)


    # Get all file names in the directory that are not directories themselves
    filenames = [f for f in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, f))]
    n_prof = len(filenames)

    # print(filenames)
    # Initialize data structures
    meta_names = ['WOD_unique_id', 'accession_number', 'dataset_id', 'lat', 'lon', 'year', 'month', 'day', 'probe type',
                  'recorder', 'hour', 'minute', 'depth_number', 'maximum_depth', 'hasTemp', 'hasSalinity', 'hasOxygen',
                  'hasChlorophyll', 'country_name', 'GMT_time', 'WMO_ID', 'dbase_orig', 'project_name', 'Platform',
                  'ocean_vehicle', 'Institute', 'WOD_cruise_identifier', 'sum_temp', 'sum_salinity', 'sum_depth',
                  'std_depth', 'std_temp', 'std_salinity', 'corr_temp_depth', 'corr_sal_depth']
    # create an empty dataframe DNA_series with meta_names as headers
    DNA_series = pd.DataFrame(columns=meta_names)

    # Process each file
    for idx, filename in enumerate(filenames):
        print(f"Processing file {idx+1}/{n_prof}: {filename}")
        file_path = os.path.join(filepath, filename)

        # Open the netCDF file as xarray dataset
        try:
            f = xr.open_dataset(file_path, decode_times=False)
        except Exception as e:
            print(f"Failed to open {filename}: {str(e)}")
            continue

        # Read 'z' (depth) data
        depth = f['z'].values

        # Apply conditions: set values outside the acceptable range to NaN
        depth = np.where((depth > 12000) | (depth < -10), np.nan, depth)

        # Calculate number of depth measurements
        DNA_series['depth_number'] = np.count_nonzero(~np.isnan(depth), axis=1)  # Using len() since depth is a numpy array

        # Determine maximum depth, taking care only to consider valid (non-NaN) entries
        if np.any(~np.isnan(depth)):
            DNA_series['maximum_depth'] = np.nanmax(depth, axis=1)
        else:
            DNA_series['maximum_depth'] = np.nan

        # Compute the sum and standard deviation of depth, rounding to four decimal places
        DNA_series['sum_depth'] = np.round(np.nansum(depth, axis=1), 4)
        DNA_series['std_depth'] = np.round(np.nanstd(depth, axis =1), 4)

        # Read 'Temperature' data
        if 'Temperature' in f:
            temp = f['Temperature'].values

            # Apply conditions: mask values outside the acceptable range
            temp = np.where((temp > 50) | (temp < -2.5), np.nan, temp)

            # Filter out NaN values for further processing
            mask = np.isnan(temp)
            temp[mask] = np.nan
            depth2 = depth
            depth2[mask] = np.nan # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid temperature readings
            DNA_series['hasTemp'] = np.any(~np.isnan(temp), axis =1)

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            DNA_series['sum_temp'] = np.round(np.nansum(temp, axis=1), 5)
            DNA_series['std_temp'] = np.round(np.nanstd(temp, axis=1), 5)
            # Calculate correlation if both arrays have non-NaN data and at least two data points
            for i in range(temp.shape[0]):
                if len(~np.isnan(temp[i])) >= 2:
                    mask = ~np.isnan(temp[i]) & ~np.isnan(depth2[i])
                    if np.sum(mask) > 1:  # Ensure there are at least two points to correlate
                        DNA_series['corr_temp_depth'][i] = np.corrcoef(temp[i, mask], depth2[i, mask])[0, 1]

            # Handle edge cases explicitly
            DNA_series['sum_temp'] = DNA_series['sum_temp'].replace(0, np.nan)
            DNA_series['std_temp'] = DNA_series['std_temp'].replace(0, np.nan)

       # Read 'Salinity' data
        if 'Salinity' in f:
            sal = f['Salinity'].values

            # Apply conditions: mask values outside the acceptable range
            sal = np.where((sal > 45) | (sal < -1), np.nan, sal)

            # Filter out NaN values for further processing
            mask = np.isnan(sal)
            sal[mask] = np.nan
            depth2 = depth
            depth2[mask] = np.nan  # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid salinity readings
            DNA_series['hasSalinity'] = np.any(~np.isnan(sal), axis=1)

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            DNA_series['sum_salinity'] = np.round(np.nansum(sal, axis=1), 5)
            DNA_series['std_salinity'] = np.round(np.nanstd(sal, axis=1), 5)
            # Calculate correlation if both arrays have non-NaN data and at least two data points
            for i in range(sal.shape[0]):
                if len(~np.isnan(sal[i])) >= 2:
                    mask = ~np.isnan(sal[i]) & ~np.isnan(depth2[i])
                    if np.sum(mask) > 1:  # Ensure there are at least two points to correlate
                        DNA_series['corr_sal_depth'][i] = np.corrcoef(sal[i, mask], depth2[i, mask])[0, 1]

            # Handle edge cases explicitly
            DNA_series['sum_salinity'] = DNA_series['sum_salinity'].replace(0, np.nan)
            DNA_series['std_salinity'] = DNA_series['std_salinity'].replace(0, np.nan)

        if 'Oxygen' in f:
            oxy = f['Oxygen'].values
            DNA_series['hasOxygen'] = np.any(~np.isnan(oxy), axis=1)

        if 'Chlorophyll' in f:
            chl = f['Chlorophyll'].values
            DNA_series['hasChlorophyll'] = np.any(~np.isnan(chl), axis=1)

        # Try to read various attributes and variables
        if 'CODA_id' in f:
            wodid = f['CODA_id'].values
            DNA_series['WOD_unique_id'] = wodid

        # Read geographic coordinates
        if 'lat' in f:
            lat = np.round(f['lat'].values, 4)
            DNA_series['lat'] = lat

        if 'lon' in f:
            lon = np.round(f['lon'].values, 4)
            DNA_series['lon'] = lon

        # Read and process date and time
        time_var = f['time'].values
        dtime = nc.num2date(time_var, f['time'].attrs['units'])
        DNA_series['year'] = dtime.year
        DNA_series['month'] = dtime.month
        DNA_series['day'] = dtime.day
        DNA_series['hour'] = dtime.hour
        DNA_series['minute'] = dtime.minute

        if 'country' in f:
            DNA_series['country_name'] = f['country'].values
            country_name = bytes(country_name[~country_name.mask]).decode('ascii')

        try:
            probe_type = f.variables['Temperature_Instrument'][:]
            probe_type = bytes(probe_type[~probe_type.mask]).decode('ascii')
        except:
            probe_type = ''

        try:
            need_z_fix = f.variables['need_z_fix'][:]
            need_z_fix = bytes(need_z_fix[~need_z_fix.mask]).decode('ascii')
        except:
            need_z_fix = ''

        try:
            recorder = f.variables['Recorder'][:]
            recorder = bytes(recorder[~recorder.mask]).decode('ascii')
        except:
            recorder = ''

        try:
            GMT_time = np.float64(f.variables['GMT_time'][:])
        except:
            GMT_time = np.nan

        try:
            WMO_id = np.float64(f.variables['WMO_ID'][:])
        except:
            WMO_id = np.nan

        # Reading and processing the 'dbase_orig' attribute
        try:
            dbase_orig = f.variables['dbase_orig'][:]
            dbase_orig = bytes(dbase_orig[~dbase_orig.mask]).decode('ascii')
        except:
            dbase_orig = ''

        # Reading and processing the 'Project' attribute
        try:
            project_name = f.variables['Project'][:]
            project_name = bytes(project_name[~project_name.mask]).decode('ascii')
        except:
            project_name = ''

        # Reading and processing the 'Platform' attribute
        try:
            platform = f.variables['platform'][:]
            platform = bytes(platform[~platform.mask]).decode('ascii')
            # print(platform)
        except:
            platform = ''

        # Reading and processing the 'Ocean_Vehicle' attribute
        try:
            ocean_vehicle = f.variables['Ocean_Vehicle'][:]
            ocean_vehicle = bytes(ocean_vehicle[~ocean_vehicle.mask]).decode('ascii')
            # print(ocean_vehicle)
        except:
            ocean_vehicle = ''

        # Reading 'Access_no'
        try:
            accession_number = f.variables['Access_no'][:]
        except:
            accession_number = np.nan

        # Reading and processing the 'Institute' attribute
        try:
            institute = f.variables['Institute'][:]
            institute = bytes(institute[~institute.mask]).decode('ascii') 
            # print('Institute ='+institute)       
        except:
            institute = ''

        # Reading and processing the 'WOD_cruise_identifier' attribute
        try:
            wod_cruise_identifier = f.variables['WOD_cruise_identifier'][:]
            wod_cruise_identifier = bytes(wod_cruise_identifier[~wod_cruise_identifier.mask]).decode('ascii')         
        except:
            wod_cruise_identifier = ''
        # print(wod_cruise_identifier)


        try:
            dataset_name = f.variables['dataset'][:]
            dataset_name = bytes(dataset_name[~dataset_name.mask]).decode('ascii')         
            # dataset_name = ''.join(chr(x) for x in f.variables['dataset'][:]).strip()
            dataset_name=dataset_name.lower()
            if any(sub in dataset_name for sub in ['bod', 'bottle', 'ocean station', 'osd', 'low-resolution', 'low resolution']):
                dataset_id = 1
            elif any(sub in dataset_name for sub in ['towed', 'uor', 'undulating']):
                dataset_id = 10
            elif any(sub in dataset_name for sub in ['ctd', 'xctd']):
                dataset_id = 2
            elif any(sub in dataset_name for sub in ['mbt', 'mechanica', 'mb']):
                dataset_id = 3
            elif any(sub in dataset_name for sub in ['xbt', 'xb', 'expendable']):
                dataset_id = 4
            elif 'sur' in dataset_name or 'surface' in dataset_name:
                dataset_id = 5
            elif any(sub in dataset_name for sub in ['apb', 'autonomous', 'animal']):
                dataset_id = 6
            elif any(sub in dataset_name for sub in ['mrb', 'moored', 'tao']):
                dataset_id = 7
            elif any(sub in dataset_name for sub in ['pfl', 'argo', 'profiling']):
                dataset_id = 8
            elif 'drb' in dataset_name or 'drifting' in dataset_name:
                dataset_id = 9
            elif 'gld' in dataset_name or 'glider' in dataset_name:
                dataset_id = 11
            elif 'dbt' in dataset_name:
                dataset_id = 12
            elif 'std' in dataset_name:
                dataset_id = 13
            elif 'microbt' in dataset_name:
                dataset_id = 14
            else:
                dataset_id = np.nan
        except:
            dataset_id = np.nan

        # Reading 'lat' and 'lon'
        try:
            latitude = np.round(f.variables['lat'][:], 4)
        except:
            latitude = np.nan
        try:
            longitude = np.round(f.variables['lon'][:], 4)
        except:
            ongitude = np.nan


        # Store the processed data
        ######  please make sure the 'position (order)' of each variables are consistent with the order of meta_names
        # txt[idx][0]=str(filename)
        txt[idx][8]=str(probe_type)
        txt[idx][9]=str(recorder)
        txt[idx][18]=str(country_name)
        txt[idx][21]=str(dbase_orig)
        txt[idx][22]=str(project_name)
        txt[idx][23]=str(platform)
        txt[idx][24]=str(ocean_vehicle)
        txt[idx][25]=str(institute)
        txt[idx][26]=str(wod_cruise_identifier)
        strings_columns_order=[8,9,18,21,22,23,24,25,26]

        DNA_series[idx,0:8]=[wod_unique_id, accession_number, dataset_id, latitude, longitude, year, month, day]
        DNA_series[idx,10:18]=[hour, minute, depth_number, maximum_depth, hasTemp, hasSalinity, hasOxygen,hasChlorophyll]
        DNA_series[idx,19:21]=[GMT_time,WMO_id]
        DNA_series[idx,27:35]=[sum_temp, sum_sal, sum_depth, std_depth, std_temp, std_sal, cor_temp_depth, cor_sal_depth]


        # Close the netCDF file
        f.close()


    # Delete WOD_unique_id
    DNA_series=DNA_series[:,1:]
    txt = [row[1:] for row in txt]
    del meta_names[0]   #Delete WOD_unique_id


    # Converts the string to the ASCILL code and sums
    variables_index_to_process = [x - 1 for x in strings_columns_order]
    for i in range(n_prof):
        for j in variables_index_to_process:
            if j < len(txt[i]): 
                # sum of all ACILL for each string varaible
                ASCII_sum = sum(ord(char) for char in txt[i][j] if char != ' ')
                DNA_series[i][j] = ASCII_sum


    ###### check output folders
    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    # Get the parent directory of the current script's directory
    parent_directory = os.path.dirname(script_directory)
    # Define the path for the 'Input_files' directory within the parent directory
    input_files_directory = os.path.join(parent_directory, 'Input_files')
    # Check if the 'Input_files' directory exists, and create it if it does not
    if not os.path.exists(input_files_directory):
        os.makedirs(input_files_directory)

    # Save the data in .npz format
    output_filename=os.path.join(input_files_directory,'DNA_summary.npz')
    np.savez(output_filename, DNA_series=DNA_series, txt=txt,meta_names=meta_names,filenames=filenames)

    print('The DNA formatted file are output to current folder: '+input_files_directory+'\n')
    print('The DNA filename is: '+output_filename)


if __name__ == '__main__':
    # filepath = '../Input_files/WOD18_sample_1995/'
    # take input arguments and convert to file path
    if len(sys.argv) > 1:
        filepath = sys.argv[1].lstrip().rstrip()
    else:
        filepath = input("Please enter the path to your netCDF files: ").lstrip().rstrip()
    
    if validate_path(filepath):
        read_netCDF_formatted_DNA_series(filepath)
    else:
        print("The entered path is not valid. Please ensure the path is correct and try again.")
