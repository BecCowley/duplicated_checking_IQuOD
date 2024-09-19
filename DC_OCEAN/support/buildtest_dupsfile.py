#!/usr/bin/env python3
import xarray as xr


def extract_subset_from_netcdf():
    """
    Extract a subset of a netCDF file and save it to a new netCDF file.

    Parameters:
    inputfile (str): The path to the input netCDF file.
    outputfile (str): The path to the output netCDF file.
    """
    inputfile = '/oa-decadal-climate/work/observations/CARSv2_ancillary/CODA/CODAv1/2011/MNF_CODA_2011_ctd.nc'
    outputfile = '/oa-decadal-climate/work/observations/CARSv2_ancillary/CODA/CODAv1/testdups/MNF_CODA_2011_ctd.nc'

    print(inputfile)

    # Open the netCDF file
    ds = xr.open_dataset(inputfile)

    # change z to positive values (temporary fix for MNF and RAN data only)
    ds['z'] = ds['z'] * -1

    # Extract a subset of the data where Platform contains 'Southern Surveyor'
    subset = ds.where(ds.Platform.str.lower().str.contains('southern'), drop=True)

    # Save the subset to a new netCDF file
    subset.to_netcdf(outputfile)

    # Close the netCDF file
    ds.close()

extract_subset_from_netcdf()
