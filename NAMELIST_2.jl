
using NCDatasets
using Dates

#-------------------------------------------------
#In this script you will to define all of the necessary paths and constants
#for your WCB case.
#-------------------------------------------------


#-------------------------------------------------
#       !!!! IMPORTANT !!!! 
#       
#If you are reading this before you have run the merge and sum scripts, then 
#you need to do that first.
#-------------------------------------------------


#-------------------------------------------------
#       Other notes
#
# -The trajectory files for some domain X must write the value "0" if the 
#  trajectory is currently not in that domain. If instead it writes -999 then
#  you have to modify the trajectory files first.
#
# -The "rtime" variable in the traj_* files MUST BE IN SECONDS!
#
# -The "rtime" output MUST be in a time step, such that if one divides 48 by 
#  the time step in hours, it is an integer. 
#  For example: a half hour time step = 0.5, and 48 / 0.5 = 96 is an integer. 
#  However, if for whatever reason, rtime is written for example every 2800 seconds,
#  this this is a 0.7777 hour time step, and 48 / 0.7777 is not an integer. Then
#  this script will NOT work.
#
# -The pressure MUST BE IN Pa! NOT hPa! 
#
# -Latitude and Longitude in the trajectory files MUST be in radians, not degrees
#-------------------------------------------------


#-------------------------------------------------
#       Define paths
#-------------------------------------------------

global TRAJ_PATH        = "/lustre/project/nhr-tpchange/gziyan/online_trajs/case02_2011121106/"
global NWP_PATH         = "/lustre/project/nhr-tpchange/amiltenb/wcb_case_2011121106/"
global CODE_PATH        = "/lustre/project/nhr-tpchange/gziyan/WCB_julia_code/"

#-------------------------------------------------
#       Define filename convention
#-------------------------------------------------

#If your NWP files look like: mod_NWP_global_DOM03_20150825T003000Z.nc, then
global NWP_FILE_FORMAT  = "mod_NWP_global_DOM0X_YYYYMMDDTHHMMSSZ.nc"
#exactly like this (don't replaxe X, YYYY, MM and so on). X stands for domain,
#YYYY for year, MM for month, DD for day and so on. 

#If your traj files look like: traj_01_tst00000325_p001_dom001.nc, then
global TRAJ_FILE_FORMAT = "traj_0N_tst_TTTTTTTT_p00W_dom00X.nc"

#-------------------------------------------------
#       Define variable lists
#-------------------------------------------------

#List of variables in your NWP datasets.

global NWP_VAR_LIST = ["lon", "lat", "height", "height_bnds", "height_2", "plev", "plev_bnds", "plev_2", "plev_2_bnds", "plev_3", "plev_3_bnds", "height_3", "time", "u", "v", "w", "temp", "pres", "rho", "pv_full", "qv", "qc", "qi", "qr", "qs", "qg", "qnc", "qnr", "qni", "qns", "qng", "tqr", "tqc", "tqs", "tqi", "tqg", "pres_sfc", "pres_msl", "cape_ml", "cin_ml", "tot_prec", "prec_con", "lhfl_s", "shfl_s", "clc", "clct", "clch", "clcl", "clcm", "t_2m", "sob_t", "sob_s", "thb_t", "thb_s"]

#List of variables in your trajectory files. 

global TRAJ_VAR_LIST = ["idx", "CellId", "lon", "lat", "alt", "rtime", "w_v", "rho", "t", "p", "pv", "tmin", "tmax", "qv", "qc", "qr", "qi", "qs", "qg", "qh", "qnc", "qnr", "qni", "qns", "qng", "qnh", "qi_in", "qs_in", "qg_in", "qh_in", "qr_in", "qi_out", "qs_out", "qg_out", "qh_out", "qr_out", "r_evap", "f_evap", "evap", "subl", "r_melt", "c_melt", "cond", "depo", "qc_nuc", "qi_hom", "qihh", "wbf", "qx_dep", "satad2", "r_frez", "qx_rim", "ni_hom", "ni_hh", "nc_nuc", "n_sim", "qvturc", "qcturc", "qiturb", "qvconv", "qcconv", "qrconv", "qiconv", "qsconv", "qvturb", "qcturb"]

#Of the variables in TRAJ_VAR_LIST, which ones are "microphysical" variables,
#meaning that they are accumulated along the simulation. 

global MICRO_VAR_LIST = ["qi_in", "qs_in", "qg_in", "qh_in", "qr_in", "qi_out", "qs_out", "qg_out", "qh_out", "qr_out", "r_evap", "f_evap", "evap", "subl", "r_melt", "c_melt", "cond", "depo", "qc_nuc", "qi_hom", "qihh", "wbf", "qx_dep", "satad2", "r_frez", "qx_rim", "ni_hom", "ni_hh", "nc_nuc", "n_sim", "qvturc", "qcturc", "qiturb", "qvconv", "qcconv", "qrconv", "qiconv", "qsconv", "qvturb", "qcturb"]
#-------------------------------------------------
#       Define start and end time and time steps
#-------------------------------------------------

global START_OF_SIM     = DateTime(2011,12,11,6,0,0) #year, month, day, hour, minute, second
global END_OF_SIM       = DateTime(2011,12,15,12,0,0)

#Number of seconds after START_OF_SIM after which trajectory files start
#you can figure this number out by looking at the first (non-zero) "rtime" 
#value in your first trajectory file.  
global START_OF_TRAJS::Int      = 9780 #IN SECONDS

#Time step of trajectory files in seconds. 
global TRAJ_TIME_STEP::Int      = 1800 #IN SECONDS
#this then calculates the time step in hours
global TRAJ_TIME_STEP_HOURS = TRAJ_TIME_STEP / 3600

#How many time steps does the first trajectory file have?
global NUMBER_RTIME_STEPS::Int  = 199

#-------------------------------------------------
#       Cyclone check
#-------------------------------------------------

#Here you need to define two times, where you define a longitude and latitude
#box within which the trajectories must be located to count as part of the WCB
#that you want to look at. For this, you must do some pre-vizualizations with 
#the NWP data.

#The box is defined by lo1,lo2,la1,la2.
#lo1 < lo2, la1 < la2. 
#
#       THIS WILL NOT WORK WHEN THE WCB CROSSES 0°/180° longitude.
#       YOU WILL NEED TO WORK THAT OUR YOURSELF. BLAME EUROCENTRICITY.

#  (lo1,la2)-------------------------------(lo2,la2)
#  |               ooooooooo               |       
#  |                    oooooooooooo       |            
#  |                       ooooooooooooo   |                 
#  |                         ooooooooooooo |             
#  |                   ooooooooooooooooooo |                   
#  |              ooooooooooooooooooooooo  |                       
#  |        oooooooooooooooooooooooooooo   |                            
#  |      oooooooooooooooooooooooo         |                          
#  |    ooooooooooooooo                    |               
#  |  oooooo                               |     
#  |                                       | 
#  (lo1,la1)-------------------------------(lo2,la1)

global TIME_CHECK_1 = DateTime(2011,12,12,7,0) #enter first check time here
global LON_CHECK_1 = (-45,-5) #in degrees, enter values here
global LAT_CHECK_1 = (35,63) #in degrees, enter values here

global TIME_CHECK_2 = DateTime(2011,12,12,18,0) #enter second check time here
global LON_CHECK_2 = (-40,0) #in degrees, enter values here
global LAT_CHECK_2 = (40,63) #in degrees, enter values here

print("Defined all global variables!")



