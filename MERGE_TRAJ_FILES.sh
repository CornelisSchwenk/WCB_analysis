#!/bin/bash

#--------------------------------------------
#In this script you will be able to merge the trajectory files using CDO.
#
#Trajectory files must be named like so:
#
#	traj_01_tst00000325_p001_dom001.nc
#	        ^9th pos.             ^31st pos. (bash counts from 1)
#
#Note, this script can take a while! Probably more than an hour.
#--------------------------------------------

#write paths to your data and where you want the trajectory files to go
path_NWP_data="path/to/ICON/NWP/data/"
path_traj_data="path/where/you/want/merged/files/"

#go to NWP data directory
cd ${path_NWP_data}

#create list of traj_01 files
traj_list=($(ls traj_01*))

#Merge traj_01.., traj_02.., traj_03.., etc. files
for i in ${traj_list[@]}
do
	echo ${name}
	name=$(echo $i | cut -c 9-31) #this picks out the "tst00000325_p001_dom001" part
        cdo merge ${p_n}*$name* ${p_h}traj_total_${name}.nc
done

echo "------------------------"
echo "Done merging, now summing files"
echo "------------------------"

#Sum the merged files across domains (traj_total_..._dom01 + traj_total_..._dom02 + ...) 
#If you have more or less than 3 domains, then adjust below accordingly
tot_list=$(ls traj_total*dom001*)
for j in $(echo $tot_list)
do
	name=$(echo $j | cut -c 12-27)
	echo ${name}
	cdo enssum traj_total_${name}_dom001.nc traj_total_${name}_dom002.nc traj_total_${name}_dom003.nc traj_sum_${name}.nc
done
