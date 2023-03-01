'''
Author: Alex Crawford
Date Created: 10 Mar 2019
Date Modified: 22 Aug 2019 -- Update for Python 3
               01 Apr 2020 -- Switched output to netCDF instead of GeoTIFF;
                               no longer dependent on gdal module (start V5)
               19 Oct 2020 -- pulled the map creation out of the for loop
               06 Oct 2021 -- added a wrap-around for inputs that prevents
                              empty cells from forming along either 180° or
                              360° longitude (start V6)
               15 Nov 2022 -- replaced "np.int" with "int"

Purpose: Reads in netcdf files & reprojects to the NSIDC EASE2 Grid North.

Edited: Lara Tobias-Tarsh Feb 28th 2023
Takes variables from a bash script of running python from script
'''

'''********************
Import Modules
********************'''
print("Loading modules.")
import os
from sys import argv
import numpy as np
from netCDF4 import Dataset
import xesmf as xe
import CycloneModule_13_2 as md

'''********************
Define Variables
********************'''
print("Defining variables")


'''********************
Parse Script Variables
********************'''
print("Parsing Script Variables")
print("ORDER: ra,var,ymin,ymax,mmin,mmax,dmin,dmax,size,verd,vert")
# File Variables:
ra = string(argv[1])
var = string(argv[2])

if ra == "ERA5":
    ncvar = "msl"
else:
    ncvar = "msl"
nctvar = "time"
ncext = '.nc'

# Time Variables
ymin, ymax = int(argv[3]), int(argv[4])
mmin, mmax = int(argv[5]), int(argv[6])
dmin, dmax = int(argv[7]), int(argv[8])

mons = ["01","02","03","04","05","06","07","08","09","10","11","12"]
dpm = [31,28,31,30,31,30,31,31,30,31,30,31] # days per month (non leap year)
timestep = 1 # in hours
startdate = [1900,1,1] # The starting date for the reanalysis time steps

# Inputs for reprojection
xsize, ysize =int(argv[9])*100, -int(argv[9])*100 # in meters
# select number of grid cells based on the grid size:
if argv[9] == 100:
    nx, ny = 180, 180 # number of grid cells; use 180 by 180 for 100 km grid
else if argv[9] == 50:
    nx, ny = 360, 360
else if argv[9] == 25:
    nx, ny = 720, 720
else:
    raise ValueError(f'Invalid Grid Size: {argv[9]}. Choose from 100, 50, 25km')

# Path Variables
path = "/home/laratt/Documents/cycloneTracking/cyclonetracking"
inpath = path+"/"+ra+"/"+var #
outpath = path+"/"+ra+"/"+var+"_EASE2_N0_"+str(int(xsize/1000))+"km" #
suppath = path+"/Projections"


'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Obtain list of nc files:
os.chdir(outpath)
fileList = os.listdir(inpath)
fileList = [f for f in fileList if (f.endswith(ncext) & f.startswith(ra))]

# Identify the time steps:
ref_netcdf = Dataset(inpath+"/"+fileList[-1])

# Create latitude and longitude arrays:
lons = ref_netcdf.variables['longitude'][:]
lats = ref_netcdf.variables['latitude'][:]

print(f'\n')
print(lons)
print(f'\n')
print(lats)
print(f'\n')

outprjnc = Dataset(suppath+'/EASE2_N0_'+str(int(xsize/1000))+'km_Projection.nc')
outlat = outprjnc['lat'][:].data
outlon = outprjnc['lon'][:].data

print(f'\n')
print(outlon)
print(f'\n')
print(outlat)
print(f'\n')

# Close reference netcdf:
ref_netcdf.close()

# Define Grids as Dictionaries
grid_in = {'lon': np.r_[lons,lons[0]], 'lat': lats}
grid_out = {'lon': outlon, 'lat': outlat}


print(f'\n')
print(grid_in)
print(f'\n')
print(grid_out)
print(f'\n')

# Create Regridder
regridder = xe.Regridder(grid_in, grid_out, 'bilinear')

print("Step 2. Set up dates of analysis")
years = range(ymin,ymax+1)
mos = range(mmin,mmax+1)
hrs = [h*timestep for h in range(int(24/timestep))]

ly = md.leapyearBoolean(years) # annual boolean for leap year or not leap year

# Start the reprojection loop
print("Step 3. Load, Reproject, and Save")
for y in years:
    Y = str(y)

    for m in mos:
        M = mons[m-1]

        mlist, hlist = [], []
        ncList = [f for f in fileList if Y+M in f]

        if len(ncList) > 1:
            print("Multiple files with the date " + Y+M + " -- skipping.")
            continue
        if len(ncList) == 0:
            print("No files with the date " + Y+M + " -- skipping.")
        else:
            nc = Dataset(inpath+"/"+ncList[0])
            tlist = nc.variables[nctvar][:]

            # Restrict days to those that exist:
            if m == 2 and ly[y-ymin] == 1 and dmax > dpm[m-1]:
                dmax1 = 29
            elif dmax > dpm[m-1]:
                dmax1 = dpm[m-1]
            else:
                dmax1 = dmax

            # For days that DO exist:
            for d in range(dmin,dmax1+1):
                timeD = md.daysBetweenDates(startdate,[y,m,d])*24

                print(" " + Y + " " + M + " " + str(d))

                for h in hrs:
                    # Establish Time
                    timeH = timeD + h

                    # Read from netcdf array
                    inArr = nc.variables[ncvar][np.where(tlist == timeH)[0][0],:,:]

                    # Transform data
                    outArr = regridder(np.c_[inArr,inArr[:,0]])
                    outArr[outlat < 0] = np.nan # Limits to Northern Hemisphere

                    # Add to list
                    mlist.append(outArr)
                    hlist.append(timeH)

        # Write monthly data to netcdf file
        ncf = Dataset(ra+"_EASE2_N0_"+str(int(xsize/1000))+"km_"+var+"_Hourly_"+Y+M+".nc", 'w')
        ncf.description = 'Mean sea-level pressure from ERA5. Projection specifications\
        for the EASE2 projection (Lambert Azimuthal Equal Area;\
        lat-origin = 90°N, lon-origin=0°, # cols = ' + str(nx) + ',\
        # rows = ' + str(ny) + ', dx = ' + str(xsize) + ', dy = ' + str(ysize) + ', units = meters'
        ncf.source = 'netCDF4 python module'

        ncf.createDimension('time', len(mlist))
        ncf.createDimension('x', nx)
        ncf.createDimension('y', ny)
        ncft = ncf.createVariable('time', int, ('time',))
        ncfx = ncf.createVariable('x', np.float64, ('x',))
        ncfy = ncf.createVariable('y', np.float64, ('y',))
        ncfArr = ncf.createVariable(ncvar, np.float64, ('time','y','x'))

        try:
            ncft.units = nc.variables[nctvar].units
        except:
            ncft.units = 'hours since 1900-01-01 00:00:00.0'

        ncfx.units = 'm'
        ncfy.units = 'm'
        ncfArr.units = 'Pa'

        # For x and y, note that the upper left point is the edge of the grid cell, but
        ## for this we really want the center of the grid cell, hence dividing by 2.
        ncft[:] = np.array(hlist)
        ncfx[:] = np.arange(-xsize*(nx-1)/2, xsize*(nx-1)/2+xsize, xsize)
        ncfy[:] = np.arange(-ysize*(ny-1)/2, ysize*(ny-1)/2+ysize, ysize)
        ncfArr[:] = np.array(mlist)

        ncf.close()

print("Complete.")


"""

'''*******************************************
Main Analysis
*******************************************'''
print("Main Analysis")

# Obtain list of nc files:
os.chdir(outpath)
fileList = os.listdir(inpath)
fileList = [f for f in fileList if (f.endswith(ncext) & f.startswith(ra))]

# Identify the time steps:
ref_netcdf = Dataset(inpath+"/"+fileList[-1])

# Create latitude and longitude arrays:
lons = ref_netcdf.variables['longitude'][:]
lats = ref_netcdf.variables['latitude'][:]

outprjnc = Dataset(suppath+'/EASE2_N0_'+str(int(xsize/1000))+'km_Projection.nc')
outlat = outprjnc['lat'][:].data
outlon = outprjnc['lon'][:].data

# Close reference netcdf:
ref_netcdf.close()

# Define Grids as Dictionaries
grid_in = {'lon': np.r_[lons,lons[0]], 'lat': lats}
grid_out = {'lon': outlon, 'lat': outlat}

# Create Regridder
regridder = xe.Regridder(grid_in, grid_out, 'bilinear')

print("Step 2. Set up dates of analysis")
years = range(ymin,ymax+1)
mos = range(mmin,mmax+1)
hrs = [h*timestep for h in range(int(24/timestep))] 

ly = md.leapyearBoolean(years) # annual boolean for leap year or not leap year

# Start the reprojection loop
print("Step 3. Load, Reproject, and Save")
for y in years:
    Y = str(y)
        
    for m in mos:
        M = mons[m-1]
        
        mlist, hlist = [], []
        ncList = [f for f in fileList if Y+M in f]
    
        if len(ncList) > 1:
            print("Multiple files with the date " + Y+M + " -- skipping.")
            continue
        if len(ncList) == 0:
            print("No files with the date " + Y+M + " -- skipping.")
        else:
            nc = Dataset(inpath+"/"+ncList[0])
            tlist = nc.variables[nctvar][:]
            
            # Restrict days to those that exist:
            if m == 2 and ly[y-ymin] == 1 and dmax > dpm[m-1]:
                dmax1 = 29
            elif dmax > dpm[m-1]:
                dmax1 = dpm[m-1]
            else:
                dmax1 = dmax
                
            # For days that DO exist:
            for d in range(dmin,dmax1+1):
                timeD = md.daysBetweenDates(startdate,[y,m,d])*24
                
                print(" " + Y + " " + M + " " + str(d))
                
                for h in hrs:
                    # Establish Time
                    timeH = timeD + h
                    
                    # Read from netcdf array
                    inArr = nc.variables[ncvar][np.where(tlist == timeH)[0][0],:,:]
                    
                    # Transform data
                    outArr = regridder(np.c_[inArr,inArr[:,0]])
                    outArr[outlat < 0] = np.nan # Limits to Northern Hemisphere
                                        
                    # Add to list
                    mlist.append(outArr)
                    hlist.append(timeH)

        # Write monthly data to netcdf file
        ncf = Dataset(ra+"_EASE2_N0_"+str(int(xsize/1000))+"km_"+var+"_Hourly_"+Y+M+".nc", 'w')
        ncf.description = 'Mean sea-level pressure from ERA5. Projection specifications\
        for the EASE2 projection (Lambert Azimuthal Equal Area;\
        lat-origin = 90°N, lon-origin=0°, # cols = ' + str(nx) + ',\
        # rows = ' + str(ny) + ', dx = ' + str(xsize) + ', dy = ' + str(ysize) + ', units = meters'
        ncf.source = 'netCDF4 python module'
        
        ncf.createDimension('time', len(mlist))
        ncf.createDimension('x', nx)
        ncf.createDimension('y', ny)
        ncft = ncf.createVariable('time', np.int, ('time',))
        ncfx = ncf.createVariable('x', np.float64, ('x',))
        ncfy = ncf.createVariable('y', np.float64, ('y',))
        ncfArr = ncf.createVariable(ncvar, np.float64, ('time','y','x'))
        
        try:
            ncft.units = nc.variables[nctvar].units
        except:
            ncft.units = 'hours since 1900-01-01 00:00:00.0'
        
        ncfx.units = 'm'
        ncfy.units = 'm'
        ncfArr.units = 'Pa'
        
        # For x and y, note that the upper left point is the edge of the grid cell, but
        ## for this we really want the center of the grid cell, hence dividing by 2.
        ncft[:] = np.array(hlist)
        ncfx[:] = np.arange(-xsize*(nx-1)/2, xsize*(nx-1)/2+xsize, xsize)
        ncfy[:] = np.arange(-ysize*(ny-1)/2, ysize*(ny-1)/2+ysize, ysize)
        ncfArr[:] = np.array(mlist)
        
        ncf.close()

print("Complete.")
"""