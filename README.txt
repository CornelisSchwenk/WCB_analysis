This repository is for the analysis of WCB trajectory data. 

Created by Cornelis Schwenk (corny.schwenk@gmail.com) during my work at the Johannes Gutenberg University Mainz in the Group of Annette Miltenberger 

File structure:
	- NAMELIST_X.jl (X = 1,2,3,...)
		Define paths to data, variable names, start of simulation, end of simulation, number of timesteps in simulation, time step intervals in seconds, longitude and latitude of WCB at times used in the selection algorithm

	- CASE.jl
		Define X (1,2,3,...) which denotes which WCB case you are looking at (if you are analysing multiple simulations). Each case needs to have its own path for the trajectory data, which you put into NAMELIST_X.jl

	- analysis_mogon.slurm
		SLURM file which merges the trajectory files (requires the climate data operators (CDO) module) and then runs run_WCB_selection.jl which uses functions from WCB_selection.jl to write netCDF files containing trajectory data and microphysical variables such as the condensation ratio (CR) and precipitation efficiency (PE)

	- run_WCB_selection.jl
		Julia script that runs the WCB selection and writing of NetCDF files using functions from WCB_selection.jl

	- WCB_selection.jl
		Julia script creating functions for writing NetCDF files containing trajectory data

	- WCB_functions.jl
		Julia script defining functions for trajectory data calculation (PE, CR, hydrometeor flux, etc.)

	- use_config.jl
		include this when using Julia REPL for user friendliness

	- statistical_functions.jl
		Defines some functions for the statistical manipulation of data.

	- thermodynamic_functions.jl
		Defines important thermodynamic functions like RH_i and RH_w.

	- install_packages.jl
		Install necessary packages

Instructions

install neccessary packages
	julia --project=. install_packages.jl

define case number in CASE.jl
define paths and variables etc. in NAMELIST_X.jl
modify paths in analysis_mogon.slurm and run:
	sbatch analysis_mogon.slurm
Alternatively run the commands from analysis_mogon.slurm in the command line. 
