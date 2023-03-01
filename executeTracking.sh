#!/bin/bash
## exit on the first error.
set -e

reproject() {

    echo "Reprojecting SLP fields with parameters: "
    echo "ra: $1"
    echo "var: $2"
    echo "start date: $3 / $5 / $7"
    echo "end date: $4 / $6 / $8"
    echo "grid size: $9 km"

    python3 /home/laratt/Documents/cycloneTracking/scriptTracking/Version 13_2 Scripts/C2_Reprojection6_E5.py "$@"
    echo "Completed Reprojection"
}

cycloneDetection() {

    echo "Starting Cyclone Tracking with variables: "
    echo "ra: $1"
    echo "var: $2"
    echo "start date: $3 / $5 / $7"
    echo "end date: $4 / $6 / $8"
    echo "grid size: $9 km"
    echo "kernel size:" ${13}
    
    python3 /home/laratt/Documents/cycloneTracking/scriptTracking/Version 13_2 Scripts/C3_CycloneDetection_13_2.py "$@"
    echo "Completed cyclone detection"
}

systemDetection() {
    echo "Starting System Tracking with variables: "
    echo "ra: $1"
    echo "var: $2"
    echo "start date: $3 / $5 / $7"
    echo "end date: $4 / $6 / $8"
    echo "grid size: $9 km"
    echo "existing tracks: ${12}"

    python3 /home/laratt/Documents/cycloneTracking/scriptTracking/Version 13_2 Scripts/C3_SystemDetection_13.py "$@"
    echo "Completed system detection"
}

subset() {
    echo "Starting System Tracking with variables: "
    echo "min lifespan ${16}"
    echo "min track length ${17}"
    echo "min latitude ${18}"
    echo "max latitude ${19}"
    echo "min longitude ${20}"
    echo "min longitude ${21}"

    python3 /home/laratt/Documents/cycloneTracking/scriptTracking/Version 13_2 Scripts/C4_Subset_Crossing_byLatLon_andLength_13.py "$@"
    echo "Completed subsetting cyclones"

}

exportcsv() {
    echo "Exporting cyclone objects to csv:"
    python3 /home/laratt/Documents/cycloneTracking/scriptTracking/Version 13_2 Scripts/C17_ExportToCSV_V13.py "$@"
    echo "Completed exporting objects"
}


echo "Enter variables:\n Order: ra,var,ymin,ymax,mmin,mmax,dmin,dmax,size, verd,vert,prior,ksizekm,
bboxnum,type,minlifespan,mintracklength,minlat, maxlat, minlon, maxlon, reproject(y/n)\n"

## standard example:
## ./executeTracking.sh ERA5 msl 1979 2021 1 1 1 1 100 13_2 EX 0 200 01 System 1 1000 [add these lol]

if [${22} == "y"]
then
reproject()
echo "Starting Reprojection Sequence..."
else
echo "Starting Reprojection Sequence..."

cycloneDetection()
systemDetection()
subset()
exportcsv()
