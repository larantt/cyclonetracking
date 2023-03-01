# Progress for Cyclone Tracking on Linux
## Dependencies
* This doesn't work in newer python. You need to do this in Python 3.8.5 or xesmf will fail.

SSH key to csrwks:
```
ssh laratt@csrwks2019-0067.engin.umich.edu
```

Create the environment that lets you reproject in Linux:
```
conda create -n cycTrackingLinux -c conda-forge python=3.8.5 numpy=1.20.0 netcdf4=1.5.5.1 scipy=1.4.1 matplotlib=3.3.4 basemap=1.2.2 cartopy=0.20.2 xesmf=0.5.1 cdsapi=0.5.1 metpy=1.3.1
```

Run the scripts on this machine and nohup output:\n
Reprojection:\n
```
nohup python3 /home/laratt/Documents/cycloneTracking/cyclonetracking/'Version 13_2 Scripts'/C2_Reprojection6_E5.py &
```
Cyclone Detection:
```
nohup python3 /home/laratt/Documents/cycloneTracking/cyclonetracking/'Version 13_2 Scripts'/C3_CycloneDetection_13_2.py &
```
System Detection:
```
nohup python3 /home/laratt/Documents/cycloneTracking/cyclonetracking/'Version 13_2 Scripts'/C3_SystemDetection_13.py &
```

## Runs
* Feb 28th 2023 - Run in 100km on Linux machine
