#!/usr/bin/env python3

########################################################################
### @copyright Copyright (C) 2024 All Rights Reserved.
### @brief 
### @version
###		Date	|	Author			|	Version		|	Description
### ------------|-------------------|---------------|---------------
### 2024-03-26	|                   |	1.1			|	Create
### 2024-06-01	|                   |	1.2			|	Create
### 2024-07-04	|                   |	1.3			|	Create
######################################################################


"""
This script used for pre-processing profile data (only supports WOD18 netCDF format files).
This script aims at reading the metadata and secondary information from the original netCDF files (we use the WOD18 single netCDF format) and then processing the metadata to create a data-metadata full list

During the process:
--Numerical metadata are retained.
--String metadata are converted into numerical values by using the ASCII code to convert each letter (including spaces) and then add the ASCII code of each letter to obtain final values.
The results are stored in *.npz format file.

***For data that is not in WOD18 netCDF format, it needs to be rewritten as WOD18 format firstly.
The variables included in the WOD18 netCDF files can be checked in Table 3 in the README.md

we used the following metadata and secondly stat information to calculate the Profile Summary Score (PSS) for each profiles
meta_names = ['WOD_unique_id','accession_number', 'dataset_id', 'lat', 'lon', 'year', 'month', 'day', 'probe type',
              'recorder', 'hour', 'minute', 'depth_number', 'maximum_depth', 'hasTemp', 'hasSalinity', 'hasOxygen',
              'hasChlorophyll', 'country_name', 'GMT_time', 'WMO_ID', 'dbase_orig', 'project_name', 'Platform',
              'ocean_vehicle', 'Institute', 'WOD_cruise_identifier', 'sum_temp', 'sum_salinity', 'sum_depth',
              'std_depth', 'std_temp', 'std_salinity', 'corr_temp_depth', 'corr_sal_depth']

The order in 'meta_names' CANNOT be modifed. Please strictly following this order. Keep NaN or set it as '' if this information is missing.

Read the Section 5 of README file to customize your own netCDf file (if necessary).

PSS = Profile Summary Score (See Song et al., 2024;Frontier in Mairne Science)

The calculation of the PSS will be done in the N01_Possible_Duplicate_Check.py

Usage:
    Run this script and follow the prompt to enter the directory path containing netCDF files.
"""
import os
import numpy as np
import netCDF4 as nc
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime, timedelta
import warnings
import sys
import argparse
import xarray as xr

## Used to determine whether the input path is correct
def validate_path(input_path):
    # Normalize the path
    normalized_path = os.path.normpath(input_path)

    # Check if the path exists and is a directory
    if not os.path.exists(normalized_path) or not os.path.isdir(normalized_path):
        return False

    return True


def find_dataset_item(search_string):
    # assign a number based on the dataset name
    datasets_dict = {'bod': 1, 'bottle': 1, 'ocean station': 1, 'osd': 1, 'low-resolution': 1,
                     'low resolution': 1, 'ctd': 2, 'xctd': 2,
                     'mbt': 3, 'mechanica': 3, 'mb': 3, 'xbt': 4, 'xb': 4, 'expendable': 4, 'sur': 5,
                     'surface': 5, 'apb': 6, 'autonomous': 6, 'animal': 6, 'mrb': 7, 'moored': 7, 'tao': 7,
                     'pfl': 8, 'argo': 8, 'profiling': 8, 'drb': 9, 'drifting': 9,
                     'towed': 10, 'uor': 10, 'undulating': 10, 'gld': 11, 'glider': 11,
                     'dbt': 12, 'std': 13, 'microbt': 14}
    matches = [value for key, value in datasets_dict.items() if key in search_string.lower()]
    return matches[0] if matches else 0


# Read NetCDF files and pre-processing metadata and secondary processing data. Retain numerical metadata and convert string metadata into numerical values by using the ASCII code table and then summing these ASCII code values of each letter to obtain a single value.
def read_netCDF_formatted_PSS_series_sub(inputpath, outputpath):
    warnings.filterwarnings("ignore")

    # Path where the netCDF files are stored
    inputpath = os.path.normpath(os.path.abspath(inputpath))
    print(inputpath)

    # Get all file names in the directory that are not directories themselves
    filenames = [f for f in os.listdir(inputpath) if os.path.isfile(os.path.join(inputpath, f))]

    # print(filenames)
    # Initialize data structures
    meta_names = ['filename', 'WOD_id', 'CODA_id', 'Access_no', 'dataset', 'lat', 'lon', 'datetime',
                  'Temperature_Instrument',
                  'Recorder', 'needs_z_fix', 'depth_number', 'maximum_depth', 'hasTemp', 'hasSalinity', 'hasOxygen',
                  'hasChlorophyll', 'country', 'WMO_ID', 'dbase_orig', 'Project', 'Platform',
                  'ocean_vehicle', 'Institute', 'WOD_cruise_identifier', 'sum_temp', 'sum_salinity', 'sum_depth',
                  'std_depth', 'std_temp', 'std_salinity', 'corr_temp_depth', 'corr_sal_depth']
    # create an empty dataframe PSS_series_sub with meta_names as headers
    PSS_series = pd.DataFrame(columns=meta_names)

    # number of files:
    n_files = len(filenames)

    # Process each file
    for idx, filename in enumerate(filenames):
        print(f"Processing file {idx + 1}/{n_files}: {filename}")
        file_path = os.path.join(inputpath, filename)

        # Open the netCDF file as xarray dataset
        try:
            f = xr.open_dataset(file_path, decode_times=False)
        except Exception as e:
            print(f"Failed to open {filename}: {str(e)}")
            continue

        # set up empty dataframe for this file
        PSS_series_sub = pd.DataFrame(columns=meta_names)

        # add the filename
        PSS_series_sub['filename'] = filename

        # number of profiles is number of casts
        n_prof = len(f['cast'].values)

        # Read 'z' (depth) data
        depth = f['z'].values

        # Apply conditions: set values outside the acceptable range to NaN
        depth = np.where((depth > 12000) | (depth < -10), np.nan, depth)

        # Calculate number of depth measurements
        PSS_series_sub['depth_number'] = np.count_nonzero(~np.isnan(depth),
                                                          axis=1)  # Using len() since depth is a numpy array

        # Determine maximum depth, taking care only to consider valid (non-NaN) entries
        if np.any(~np.isnan(depth)):
            PSS_series_sub['maximum_depth'] = np.nanmax(depth, axis=1)
        else:
            PSS_series_sub['maximum_depth'] = np.nan

        # Compute the sum and standard deviation of depth, rounding to four decimal places
        PSS_series_sub['sum_depth'] = np.round(np.nansum(depth, axis=1), 4)
        PSS_series_sub['std_depth'] = np.round(np.nanstd(depth, axis=1), 4)

        # Read 'Temperature' data
        if 'Temperature' in f:
            temp = f['Temperature'].values

            # Apply conditions: mask values outside the acceptable range
            temp = np.where((temp > 50) | (temp < -2.5), np.nan, temp)

            # Filter out NaN values for further processing
            mask = np.isnan(temp)
            depth2 = depth
            depth2[mask] = np.nan  # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid temperature readings
            PSS_series_sub['hasTemp'] = np.any(~np.isnan(temp), axis=1)

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            PSS_series_sub['sum_temp'] = np.round(np.nansum(temp, axis=1), 5)
            PSS_series_sub['std_temp'] = np.round(np.nanstd(temp, axis=1), 5)
            # Calculate correlation if both arrays have non-NaN data and at least two data points
            for i in range(temp.shape[0]):
                if len(~np.isnan(temp[i])) >= 2:
                    mask = ~np.isnan(temp[i]) & ~np.isnan(depth2[i])
                    if np.sum(mask) > 1:  # Ensure there are at least two points to correlate
                        PSS_series_sub['corr_temp_depth'][i] = np.corrcoef(temp[i, mask], depth2[i, mask])[0, 1]

            # Handle edge cases explicitly
            PSS_series_sub['sum_temp'] = PSS_series_sub['sum_temp'].replace(0, np.nan)
            PSS_series_sub['std_temp'] = PSS_series_sub['std_temp'].replace(0, np.nan)

        # Read 'Salinity' data
        if 'Salinity' in f:
            sal = f['Salinity'].values

            # Apply conditions: mask values outside the acceptable range
            sal = np.where((sal > 45) | (sal < -1), np.nan, sal)

            # Filter out NaN values for further processing
            mask = np.isnan(sal)
            depth2 = depth
            depth2[mask] = np.nan  # Assuming 'depth' is already defined and processed similarly

            # Check if there are valid salinity readings
            PSS_series_sub['hasSalinity'] = np.any(~np.isnan(sal), axis=1)

            # Compute the sum, standard deviation, and correlation coefficient, if applicable
            PSS_series_sub['sum_salinity'] = np.round(np.nansum(sal, axis=1), 5)
            PSS_series_sub['std_salinity'] = np.round(np.nanstd(sal, axis=1), 5)
            # Calculate correlation if both arrays have non-NaN data and at least two data points
            for i in range(sal.shape[0]):
                if len(~np.isnan(sal[i])) >= 2:
                    mask = ~np.isnan(sal[i]) & ~np.isnan(depth2[i])
                    if np.sum(mask) > 1:  # Ensure there are at least two points to correlate
                        PSS_series_sub['corr_sal_depth'][i] = np.corrcoef(sal[i, mask], depth2[i, mask])[0, 1]

            # Handle edge cases explicitly
            PSS_series_sub['sum_salinity'] = PSS_series_sub['sum_salinity'].replace(0, np.nan)
            PSS_series_sub['std_salinity'] = PSS_series_sub['std_salinity'].replace(0, np.nan)

        # loop over some of the variables where we just need to know if they exist or not
        for var in ['Oxygen', 'Chlorophyll']:
            if var in f:
                PSS_series_sub[f'has{var}'] = np.any(~np.isnan(f[var].values), axis=1)

        # Try to read various attributes and variables
        for var in ['WMO_id', 'lat', 'lon', 'Access_no', 'WOD_id']:
            if var in f.variables:
                PSS_series_sub[var] = np.float64(f[var].values)

        # Read and process date and time
        time_var = f['time'].values
        PSS_series_sub['datetime'] = nc.num2date(time_var, f['time'].attrs['units'])

        # Loop over variables for ascii sum
        for var in ['CODA_id', 'country', 'Temperature_Instrument', 'needs_z_fix', 'Recorder', 'dbase_orig',
                    'Project', 'Platform', 'WOD_cruise_identifier', 'Institute']:
            if var in f.variables:
                var_char = np.array([x.decode('utf-8') for x in f[var].values], dtype=str)
            else:
                var_char = np.array(['' for _ in range(n_prof)], dtype=str)
            # convert the country_name, needs_z_fix and probe_type etc to ASCII sum
            ASCII_sums = np.array([sum(ord(char) for char in string if char != ' ') for string in var_char])
            PSS_series_sub[var] = ASCII_sums

        if 'dataset' in f.variables:
            dataset_name = np.array([x.decode('utf-8') for x in f['dataset'].values], dtype=str)
        else:
            dataset_name = np.array(['' for _ in range(n_prof)], dtype=str)
        # find the matching value from the dictionary
        dataset_names = np.array([find_dataset_item(item) for item in dataset_name])
        PSS_series_sub['dataset'] = dataset_names

        # append the PSS_series_sub to the PSS_series
        PSS_series = PSS_series.append(PSS_series_sub, ignore_index=True)

        # Close the netCDF file
        f.close()

    # Delete WOD_unique_id column
    PSS_series = PSS_series.drop(columns=['WOD_id'])

    # check output folders
    script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    # Get the parent directory of the current script's directory
    parent_directory = os.path.dirname(script_directory)
    # Define the path for the 'Input_files' directory within the parent directory
    input_files_directory = os.path.join(parent_directory, 'Input_files')
    # Check if the 'Input_files' directory exists, and create it if it does not
    if not os.path.exists(outputpath):
        os.makedirs(outputpath)

    # Save the data in .npz format
    output_filename = os.path.join(outputpath, 'Profile_Summary_Score_list.npz')
    # Convert the DataFrame to a dictionary of arrays
    data_dict = {col: PSS_series[col].values for col in PSS_series.columns}
    # Save the dictionary to an .npz file
    np.savez(output_filename, **data_dict)

    print('*******************FINISHED***********************************\n')
    print('The Profile Summary Score formatted file are output to current folder: ' + outputpath + '\n')
    print('The Profile Summary Score filename is: ' + output_filename)


if __name__ == '__main__':

    OutputDir = os.path.dirname(os.path.abspath(__file__)) + "/../Input_files"
    InputDir = OutputDir + "/WOD18_sample_1995"

    parser = argparse.ArgumentParser(description='Create Profile Summary Score')
    parser.add_argument("-i", "--input", type=str, default=InputDir)
    parser.add_argument("-o", "--output", type=str, default=OutputDir)
    args = parser.parse_args()
    InputDir = args.input
    OutputDir = args.output

    # check Input/Output dir vaild
    isInputOK = validate_path(InputDir)
    isOutputOK = validate_path(OutputDir)
    if (not (isInputOK and isOutputOK)):
        print("The entered path is not valid. Please ensure the path is correct and try again.")
        raise Exception("Invalid InputDir or OutputDir!", InputDir, OutputDir)

    PSS_summary_filename = OutputDir + "/Profile_Summary_Score_list.npz"
    iAct = 1
    if os.path.exists(PSS_summary_filename):
        iAct = input("Update Profile Summary Score list or not(1: Yes (default); 0: No): ") or 1
        print(iAct)

    if (iAct == 1 or iAct == '1'):
        read_netCDF_formatted_PSS_series_sub(InputDir, OutputDir)
        print("Profile_Summary_Score_list.npz Complete !")
