#---------------------------------------------------------------------------
#       "Skript for the selection of WCB trajectories from ICON trajectory files"
#       "Created by Cornelis Schwenk 22.11.2023"
#       "c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

#       ----------   INSTRUCTIONS   ----------

#       "to write the metdata files just include this script in julia and      " #
#       "then exectute the function write_all_WCB_meta_files                   " #


#       "to write the large WCB files with variables just include this script  " #
#       "in julia and then exectute the function write_WCB_var_file            " #            

#--------------------------------------------------------
using Pkg
using Revise
using OhMyREPL
using JuMP       
using Glob, NCDatasets, Plots
using Plots.PlotMeasures
using Statistics
using BenchmarkTools, Dates, DataStructures
using DelimitedFiles
using NaNMath; nm=NaNMath
using Interpolations

include("CASE.jl")
include("NAMELIST_$(CASE).jl")

#-----------------------------------
#       "Load necessary files"
#----------------------------------

#summed trajectory files. You need to run MERGE_TRAJ_FILES.sh to have these. 
global TRAJ_FILES = glob("traj_sum*",TRAJ_PATH)

#This will be empty at first. You will write these files using this script.
global WCB_fls  = glob("WCB_tau_tst*",TRAJ_PATH)

#-----------------------------------------------------------------------
#       "define important global values and thigns"
#-----------------------------------------------------------------------

#runtime trajectory output timestep
global RT_STEP  = TRAJ_TIME_STEP / 3600 #rtime time steps in hours
global RT_0     = START_OF_TRAJS / 3600 #rtime beginning in hours

global RT_ARRAY_REALTIME        = collect(RT_0:RT_STEP:(((NUMBER_RTIME_STEPS * RT_STEP) + RT_0 - RT_STEP)))
global RT_ARRAY_TIMESTEPS       = collect(0:RT_STEP:(NUMBER_RTIME_STEPS * RT_STEP - RT_STEP))

#       "Conversion of rtime (in time steps) to DateTime"
#       "rt_to_dt(0) gives you the DateTime of the first time step in the trajectories"
rt_to_dt(rt) = START_OF_SIM + Minute(rt*60) + Second(START_OF_TRAJS)

#Enter a datetime, and it will give you the closest rt IN TIME STEPS
function dt_to_nearest_rt(dt)
    rt_dt_array = rt_to_dt.(RT_ARRAY_TIMESTEPS)
    ix = argmin(abs.(dt .- rt_dt_array))
    return ix
end

#"for writing the names in the netcdf file"
#"This is purely for name writing and if dependent on how trajectory files are named!"
global NAME_IND         = 18
global TITLE_IND        = 27

#Now define some variables that will only be able to exist once the file WCB_tau_vars_realtime.nc
#has been written
#
if length(glob("WCB_tau_vars_realtime.nc",TRAJ_PATH)) > 0
    global dsr = Dataset(TRAJ_PATH*"WCB_tau_vars_realtime.nc")
    global PR = dsr["p"][:,:] ./ 100
end


#----------------------------------------------------------------------
#       "THE THREE FOLLOWING FUNCTIONS ARE THE IMPORTANT ONES           " #
#       "THEY USE ALL OF THE OTHER FUNCTIONS DEFINED BELOW               " #
#----------------------------------------------------------------------

function write_all_WCB_meta_files()
        #"THIS IS A BIG FUNCTION. IT WRITES ALL OF THE METADATA WCB FILES."
	#"Meta WCB files give the index and positions of wcbs within the trajectory dataset."
        
	#"This function goes through each trajectory file and writes a separate ncfile that"
	#"contains the information on:" 
	#"	- index of trajs that are WCB (w)"
	#"	- where they begin their ascent (a)" 
	#"	- how long the ascent takes in time steps (t6)."
	#"This is done by the "write_WCB_ncfile" function, see for more detail on selection process"

	println("length TRAJ_FILES = $(length(TRAJ_FILES))")
        for i in 1:length(TRAJ_FILES) 	
	    f = TRAJ_FILES[i]
	    print("$(i), ")
            write_WCB_ncfile(f)
        end
	println("! done !")
end

function write_WCB_var_file_realtime()
        #"In this file the entire trajectory data from the beginning to the end of the simulation is
        #written for all WCB trajectories selected by the algorithm.
	
        #"All trajs begin at first runtime output. If trajs begin later (trajs might be started"
	#"at diffferent rtimes), then they are NaN until data is written. Example:"
	#"	traj 901 starts at rtime=0, then some_var = [1,2,3,4,...,192]"
	#" 	if traj 4000 starts at rtime=4, then some_var = [NaN,NaN,NaN,NaN,5,6,...,192]"

	#"The file WCB_tau_vars_realtime.nc is written."
	#"This uses primarily the function get_var_realtime"
        pth = TRAJ_PATH

	ntrajs = length(get_w_a_all()[1])
        ds_out_name = "WCB_tau_vars_realtime.nc"
        ds_out = Dataset(pth*ds_out_name,"c")
        defDim(ds_out,"time",NUMBER_RTIME_STEPS)
        defDim(ds_out,"ntraj",ntrajs)
        ds_out.attrib["title"] = "WCB variables in real time. To get the DateTime of rtime, take rt_to_dt(rt_value)."

        vardict = Dict()
        ncdict = Dict()
	println("------------STARTING------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
        for i in TRAJ_VAR_LIST[TRAJ_VAR_LIST .!= "rtime"]
                println("writing variable $(i)")
                ncdict[i] = defVar(ds_out,i,Float64,("ntraj","time"))
                ncdict[i][:,:] = get_var_realtime(i)
		GC.gc()
        end
	
	global W6,A6,T6 = get_w_a_all()
	t6_ = defVar(ds_out,"tau_600",Float64,("ntraj",))	
	t6_.attrib["units"] = "hours"
	t6_[:] = T6 ./ 2
	
	r_t = defVar(ds_out,"rtime",Float64,("time",))
        r_t[:] = RT_ARRAY_TIMESTEPS

        global PR = get_var_realtime("p") ./ 100
        global NS = [findfirst((!isnan).(PR[i,:]))-1 for i in 1:length(T6)]
        global AINS = [CartesianIndex(i,A6[i] + NS[i]) for i in 1:length(T6)]
        
        starts = [find_start(k) for k in 1:length(T6)]
        stops = [find_stop(k) for k in 1:length(T6)]
       
        inflow = fill(NaN,ntrajs,NUMBER_RTIME_STEPS)
        for k in 1:length(T6)
            inflow[k,1:starts[k]] .= 0
            inflow[k,starts[k]:stops[k]] .= 1
            if stops[k] < 199
                inflow[k,stops[k]+1:end] .= 2
            end
        end
        inflow_ = defVar(ds_out,"in_out",Float64,("ntraj","time"))
        inflow_[:,:] = inflow

        close(ds_out)
	println("------------FINISHED------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
end

function write_tauWCB_var_file_normed()
        #" !! This function must be executed AFTER write_WCB_var_file_realtime !!"
	#"This function won't run is WCB_tau_vars_realtime.nc doesn't exist already"

	#"This function creates the WCB_tau_vars_normed.nc file." 
	#"It writes variables in a normalized ascent time scale from 0 to 1 in 101 steps."
	#"0 is beginning of ascent. 1 is end of ascent. It takes the data from WCB_vars_realtime but"
	#"stretches/contract the ascent part we want such that all arrays have the same length."
	#"The times are still such that the first time entry is at the beginning of the ascent."
	#"rtime is still written."
	#"There are no more NaNs. If we take the example from the function above:"
	#" 	traj 901 has t6 = 2h. Then some_var = [1,2,3,4,NaN,NaN,...,NaN]"
	#" 			       and some_var_normed = [1,1.2,1.4,...,3.8,4]"

	#"the function "get_var_tWCB_normed" is the one that is primarily used here"
	
        pth = TRAJ_PATH
	
        if length(glob("WCB_tau_vars_realtime.nc",TRAJ_PATH)) == 0
            error("You must create the WCB_tau_vars_realtime.nc file first.")
	end

        if (!@isdefined PR)
            global dsr = Dataset(TRAJ_PATH*"WCB_tau_vars_realtime.nc")
            global PR = dsr["p"][:,:] ./ 100
            global W6,A6,T6 = get_w_a_all() #WCB metadata
            global NS = [findfirst((!isnan).(PR[i,:]))-1 for i in 1:length(T6)]
            global AINS = [CartesianIndex(i,A6[i] + NS[i]) for i in 1:length(T6)]
        end

	starts = [find_start(k) for k in 1:length(T6)]
        stops = [find_stop(k) for k in 1:length(T6)]

	ntrajs = length(T6)
        ds_out_name = "WCB_tau_vars_normed.nc"
        ds_out = Dataset(pth*ds_out_name,"c")
        defDim(ds_out,"time",101)
        defDim(ds_out,"ntraj",ntrajs)
        ds_out.attrib["title"] = "variables during normed WCB ascent time"

        vardict = Dict()
        ncdict = Dict()
	println("------------STARTING------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
	for i in TRAJ_VAR_LIST
		print("$(i), ")
                println("writing variable $(i)")
                ncdict[i] = defVar(ds_out,i,Float64,("ntraj","time"))
                ncdict[i][:,:] = get_var_tWCB_normed(i,starts,stops)
	end
	t6_ = defVar(ds_out,"tau_600",Float64,("ntraj",))	
	t6_.attrib["units"] = "hours"
	t6_[:] = T6 ./ 2

	tWCB = stops .- starts
	tWCB_ = defVar(ds_out,"tau_WCB",Float64,("ntraj",))	
	tWCB_.attrib["units"] = "hours"
	tWCB_[:] = tWCB ./ 2
        
	close(ds_out)
	println("------------FINISHED------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
end

function write_PECR_file()
        #" !! This function must be executed AFTER write_WCB_var_file_normed !!"
	#"This function won't run is WCB_tau_vars_normed.nc doesn't exist already"

	#"This function creates the PECR.nc file." 

        pth = TRAJ_PATH
	
        if length(glob("WCB_tau_vars_normed.nc",TRAJ_PATH)) == 0
            error("You must create the WCB_tau_vars_normed.nc file first.")
	end
        
        include(CODE_PATH*"WCB_functions.jl")

	ntrajs = length(T6)
        ds_out_name = "PECR.nc"
        ds_out = Dataset(pth*ds_out_name,"c")
        defDim(ds_out,"ntraj",ntrajs)
        defDim(ds_out,"time",101)
        ds_out.attrib["title"] = "variables at the end of tau_WCB"

        vardict = Dict()
        ncdict = Dict()
	println("------------STARTING------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
       
        println("qhy")
        qhy_ = defVar(ds_out,"qhy",Float64,("ntraj","time"))
        qhy_[:,:] = calc_qhy_all()

        println("qtot")
        qtot_ = defVar(ds_out,"qtot",Float64,("ntraj","time"))
        qtot_[:,:] = qhy_[:,:] .+ dsn["qv"][:,:]

        println("allfl")
        allfl_ = defVar(ds_out,"allfl",Float64,("ntraj","time"))
        allfl_[:,:] = flux("r")
        for hy in ["i","s","g","h"]
            allfl_[:,:] += flux(hy)
        end

        println("PE")
        PE_ = defVar(ds_out,"PE",Float64,("ntraj","time"))
        PE_[:,:] = calc_PE_time(qhy_[:,:],allfl_[:,:])
        
        println("CR")
        CR_ = defVar(ds_out,"CR",Float64,("ntraj","time"))
        CR_[:,:] = calc_CR_time()

        println("DR")
        DR_ = defVar(ds_out,"DR",Float64,("ntraj","time"))
        DR_[:,:] = calc_DR_time()	
        
        println("DR_mix")
        DR_mix_ = defVar(ds_out,"DR_mix",Float64,("ntraj","time"))
        DR_mix_[:,:] = calc_DR_mix_time(qhy_[:,:],allfl_[:,:])

        println("DR_mphys")
        DR_mphys_ = defVar(ds_out,"DR_mphys",Float64,("ntraj","time"))
        DR_mphys_[:,:] = calc_DR_mphys_time(qhy_[:,:],allfl_[:,:])

        println("RHi")
        RHi_ = defVar(ds_out,"RHi",Float64,("ntraj","time"))
        RHi_[:,:] = RH_i.(dsn["t"][:,:],dsn["qv"][:,:],dsn["p"][:,:])

        println("RHw")
        RHw_ = defVar(ds_out,"RHw",Float64,("ntraj","time"))
        RHw_[:,:] = RH_w.(dsn["t"][:,:],dsn["qv"][:,:],dsn["p"][:,:])
        
        close(ds_out)
	println("------------FINISHED------------")
	println("--------------------------------")
	println(Dates.format(now(), "HH:MM-dd-mm-yyyy"))
	println("--------------------------------")
end
#----------------------------------------------------------------------



#----------------------------------------------------------------------
#"       THE FOLLOWING FUNCTIONS ARE THE ONES USED BY THE FOUR BIG FUNCTIONS ABOVE    "#
#----------------------------------------------------------------------

function write_WCB_ncfile(fl::String)
        #"iterate this function over the file list of the trajectories"
        #"for f in nfl write_WCB_ncfile(f) end"

        #"this write only the small files with the metadata of where to find the WCB trajs"
        #"within the traj files"

	#"The selection is done by the funciton "WCB_sel" that returns index (w6), ascent start (a6)"
	#"and ascent time (t600). See it for more details on selection."

        files = TRAJ_FILES[1:length(WCB_fls)]
        pth = TRAJ_PATH

        println("going through file $(fl)")
        ds = Dataset(fl)
        wcb,asc,tau = WCB_sel(ds)
	
	#"If file contains any WCB trajectories, write netCDF file."
        if length(wcb) > 0
                ds_out_name = fl[end-NAME_IND:end]
                ds_out = Dataset(pth*"WCB_tau_$(ds_out_name)","c")
                defDim(ds_out,"pos",size(wcb)[1])
                ds_out.attrib["title"] = "WCB positions for $(fl[end-TITLE_IND:end])"

                pos = defVar(ds_out,"position",Int64,("pos",), 
			      attrib = OrderedDict("units" => "Row number in file"))
                ascent = defVar(ds_out,"ascent",Int64,("pos",),
                              attrib = OrderedDict("units" => "Array pos. at which WCB ascent begins"))
                tau600 = defVar(ds_out,"tau600",Int64,("pos",),
                              attrib = OrderedDict("units" => "Time steps until completion of ascent"))
		pos[:] = wcb
                ascent[:] = asc
		tau600[:] = tau
                println("found $(size(wcb)) WCBs")
                close(ds_out)
        elseif length(wcb) == 0
                println("found no WCB trajectories, not writing file")
        end
	close(ds)
end

function find_tau_and_asc(p)
	#"This function takes a 1 dim pressure array in hPa and returns ascent time"
	#"and where ascent begins (so asc and tau600)."
        #"p must be in hPa and a 1 dim array"
        
        # !!! THE TIME OF ASCENT IS IN TIME STEPS !!! NOT HOURS OR MINUTES
        
        DT = Inf
        t = []
        for i in 1:(length(p)-1)
		for j in i:length(p)
		       dp = p[i] - p[j]
		       dt = j-i
		       if dp > 600
			       if dt < DT
				       push!(t,(i,dt))
				       DT=dt
			       end
			       break
		       end
	       end
        end
        return t[end]
end


function find_tau_x00(p,x)
	#"This function takes a 1 dim pressure array in hPa and returns ascent time"
        #"p must be in hPa and a 1 dim array"
        
        # !!! THE TIME OF ASCENT IS IN TIME STEPS !!! NOT HOURS OR MINUTES
        
        DT = Inf
        t = []
        for i in 1:(length(p)-1)
		for j in i:length(p)
		       dp = p[i] - p[j]
		       dt = j-i
		       if dp > x
			       if dt < DT
				       push!(t,dt)
				       DT=dt
			       end
			       break
		       end
	       end
        end
        return t[end]
end

function find_WCB(ds)	
	#"This function finds trajectories with a 600hpa ascent in a netCDF dataset of trajectories"

	#"convert rtime into hours"
	rt = ds["rtime"][:] ./ 3600

	#"convert pres into hPa"
        pres = ds["p"][:,:] ./ 100

        wcb_p = Int[]
	asc = Int[]
	tau = Int[]
        
	#"find all trajectories that go up 600hPa in at least 48h"
        try Int(48/RT_STEP)
        catch
            @error "RT_STEP must be such that 48/RT_STEP is an integer!"
        end

        for i in 1:(size(rt)[1]-Int(48/RT_STEP))
                p0 = pres[:,i]
                p48 = pres[:,i+Int(48/RT_STEP)]
		
                dp = p0 .- p48
                wcbs = findall(dp .> 600)
                for x in wcbs
                        push!(wcb_p, x)
                end
        end

	#"sort arrays and make sure no double couting"
        uniq = findfirst.(isequal.(unique(wcb_p)),[wcb_p])
        wcb_sort = sortperm(wcb_p[uniq])

	wcb_out = wcb_p[uniq][wcb_sort]

	for i in 1:length(wcb_out)
		a,t600 = find_tau_and_asc(pres[wcb_out[i],:])
		push!(asc,a)
		push!(tau,t600)
	end

	return (wcb_out,asc,tau)
end

function WCB_sel(ds)
	#"This function takes the WCB ascents from "find_WCB" function and then"
	#"does a second selection process to make sure they are part of the cyclone"
	#"that we want them to be"

        println("---------------------------------------------------")
        println("""going through dataset for $(ds.attrib["time"])""")

        rt = ds["rtime"][:] ./ (3600)

        rt0 = rt[1]
        println("there will be $(size(rt)[1]) loops")

        wcbs,asc,tau = find_WCB(ds)

        lat = ds["lat"][:,:][wcbs,:] .* 180/pi
        lon = ds["lon"][:,:][wcbs,:] .* 180/pi

        cyc_mask = Int.(ones(size(lat)[1]))

	#----------------------------------
	#"The second selection is a visual confirmation for specific time." 
	#"If trajs outside of boundary at this time, mask is set to zero (or false)."
	#----------------------------------
        j_check1 = dt_to_nearest_rt(TIME_CHECK_1) #index to check 1
        j_check2 = dt_to_nearest_rt(TIME_CHECK_2) #index to check 2
	
        temp_cyc_mask2 = Int.(ones(size(lat)[1])) #all true

        #go through every trajectory and check if it is in boxes 1 or 2
        for i in 1:size(lat)[1]
            lo11 = LON_CHECK_1[1] 
            lo12 = LON_CHECK_1[2] 
            la11 = LAT_CHECK_1[1] 
            la12 = LAT_CHECK_1[2] 
            
            lo21 = LON_CHECK_2[1] 
            lo22 = LON_CHECK_2[2] 
            la21 = LAT_CHECK_2[1] 
            la22 = LAT_CHECK_2[2] 

            if !(lo11 < lon[i,j_check1] < lo12) .| !(la11 < lat[i,j_check1] < la12)
			temp_cyc_mask2[i] = 0 #remove all not in box 1
            end
            if !(lo21 < lon[i,j_check2] < lo22) .| !(la21 < lat[i,j_check2] < la22)
			temp_cyc_mask2[i] = 0 #remove all not in box 2
            end
	end

        #combine masks. Now only those in boxes 1 or 2 are left
	cyc_mask = cyc_mask .& temp_cyc_mask2

        println("")
        println("---------------------------------------------------")
        return (wcbs[Bool.(cyc_mask)],asc[Bool.(cyc_mask)],tau[Bool.(cyc_mask)])
end

function WCB_info(fl::String)
        #"This function gives you the WCB metadata for the trajectory file"
	ds = Dataset(fl)
        return (ds["position"][:],ds["ascent"][:],ds["tau600"][:])
	close(ds)
end

function get_w_a_all()
	#"This function returns the WCB metadata for ALL trajectories."
        #w = which trajectory is a WCB trajectory?
        #a = index for the beginning of the t600 ascent
        #t = t600 in time steps

        files = WCB_fls
	w_tot = Vector{Float64}()
	a_tot = Vector{Float64}()
	t6_tot = Vector{Float64}()
	for f in files
		w,a,t6 = WCB_info(f)
		w_tot = vcat(w_tot,w)
		a_tot = vcat(a_tot,a)
		t6_tot = vcat(t6_tot,t6)
	end
	return Int.(w_tot),Int.(a_tot),Int.(t6_tot)
end

if (@isdefined PR)
    global W6, A6, T6 = get_w_a_all()
    #as of which index are the values in the realtime file not NaN?
    global NS = [findfirst((!isnan).(PR[i,:]))-1 for i in 1:length(T6)]
    global AINS = [CartesianIndex(i,A6[i] + NS[i]) for i in 1:length(T6)]
end

function find_start(k)
	#"This function finds where the WCB ascent starts before t600 begins."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour before t600"
        if (!@isdefined PR)
            error("""You have to write the WCB_tau_vars_realtime.nc 
                  file first, then include this script again!""")
        elseif (@isdefined PR)
            pk = PR[k,:]
            a_s = A6
            n_s = NS 
            t_s = T6 
        end

	start_max = find_start_max(k,pk,a_s,n_s,t_s)
	start_min = find_start_min(k,pk,a_s,n_s,t_s)
	
	if start_max > start_min
		startlist = start_min:start_max
	elseif start_max < start_min
		startlist = start_max:start_min
	elseif start_max == start_min
		startlist = start_max
	end
	
	pstartlist = [pk[i] for i in startlist]
	startarg = argmax(pstartlist)
	return startlist[startarg]
	
end

function find_start_min(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	a_start = a_s[k] + n_s[k]
	if a_start == 1
		return 1
	end
	for i in 1:a_start
                p_diff_min = -(TRAJ_TIME_STEP_HOURS * 8 * 2) .* (a_start - i)
		if i == a_start
			return a_start
		elseif pk[a_start] .- pk[i] .< p_diff_min
			return i
		end	
	end
end

function find_start_max(k,pk,a_s,n_s,t_s)
	#"This function finds where the WCB ascent starts before t600 begins."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour before t600"
	if a_s[k] == 1
		return a_s[k] + n_s[k]
	elseif a_s[k] == 2
                if (pk[n_s[k] + a_s[k]] - pk[n_s[k] + a_s[k] - 1]) <= -(TRAJ_TIME_STEP_HOURS * 8)
			return n_s[k] + a_s[k] - 1
		else
			return n_s[k] + a_s[k]
		end
	else
		counter=0
		for i in reverse(n_s[k]+3:n_s[k]+a_s[k])
			counter = i-2
                        if (pk[i] - pk[i-2]) .> -(TRAJ_TIME_STEP_HOURS * 8 * 2)
				if (pk[i] - pk[i-2]) .<= 0
					return i-2
				elseif (pk[i] - pk[i-1]) .<= 0
					return i-1
				else
					return i
				end
			end
		end
		if counter == n_s[k] + 1
			return n_s[k] + 1
		end
	end
end

function find_stop(k)
	#"This function finds where the WCB ascent stops after t600 is completed."
        #"It uses find_stop_max and find_stop_min and looks for where in between"
	#"those two the pressure is minimal"
        if (!@isdefined PR)
            error("""You have to write the WCB_tau_vars_realtime.nc 
                  file first, then include this script again!""")
        elseif (@isdefined PR)
            pk = PR[k,:]
            a_s = A6
            n_s = NS 
            t_s = T6 
        end
	
        stop_max_1 = find_stop_max_1(k,pk,a_s,n_s,t_s)
	stop_max_2 = find_stop_max_2(k,pk,a_s,n_s,t_s)
	stop_max = minimum([stop_max_1,stop_max_2])	
	stop_min = find_stop_min(k,pk,a_s,n_s,t_s)

	if stop_max > stop_min
		stoplist = stop_min:stop_max
	elseif stop_max < stop_min
		stoplist = stop_max:stop_min
	elseif stop_max == stop_min
		stoplist = [stop_max]
	end
	pstoplist = [pk[i] for i in stoplist]
	stoparg = argmin(pstoplist)
	return stoplist[stoparg]

end

function find_stop_max_1(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	t6stop = a_s[k] + t_s[k] + n_s[k]
	k_to_stop = NUMBER_RTIME_STEPS - t6stop
	if k_to_stop == 0
		return NUMBER_RTIME_STEPS
	end
	for i in k_to_stop:-1:0
                p_diff_min = -(TRAJ_TIME_STEP_HOURS * 8) .* i
		if i == 0
			return t6stop
		elseif pk[t6stop + i] .- pk[t6stop] .< p_diff_min
			return t6stop + i
		end	
	end
end

function find_stop_max_2(k,pk,a_s,n_s,t_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after tWCB is completed."
	#"It goes backwards from the end and looks if the mean ascent is"
	#"faster than 8hPa per hour. This function is used by find_stop"
	tWCBstop = find_stop_min(k,pk,a_s,n_s,t_s)
	k_to_stop = NUMBER_RTIME_STEPS - tWCBstop
	if k_to_stop == 0
		return NUMBER_RTIME_STEPS
	end
	for i in k_to_stop:-1:0
                p_diff_min = -(TRAJ_TIME_STEP_HOURS * 8 .+ 0.5) .* i
		if i == 0
			return tWCBstop
		elseif pk[tWCBstop + i] .- pk[tWCBstop] .< p_diff_min
			return tWCBstop + i
		end	
	end
end

function find_stop_min(k,pk,a_s,t_s,n_s)
	#"This function finds the minimum time after which WCB ascent no"
	#"longer is larger than 8hPa/h after t600 is completed."
	#"It looks for where the pressure velocity is on average smaller than"
	#"8hPa per hour in the hour after t600. This function is used by find_stop."

        ik1 = NUMBER_RTIME_STEPS
        ik2 = NUMBER_RTIME_STEPS - 2
        ik3 = NUMBER_RTIME_STEPS - 1
        if a_s[k] + t_s[k] + n_s[k] == ik1
                if (pk[ik1] - pk[ik2]) .<= -(TRAJ_TIME_STEP_HOURS * 8 * 2)
			return ik1
                elseif (pk[ik3] - pk[ik2]) .<= -(TRAJ_TIME_STEP_HOURS * 8)
			return ik3
		else
			return ik2
		end
	elseif a_s[k] + t_s[k] + n_s[k] == ik3
		if (pk[ik1] - pk[ik3]) .<= -(TRAJ_TIME_STEP_HOURS * 8 * 2)
			return ik2
		else
			return ik3
		end
	elseif a_s[k] + t_s[k] + n_s[k] == ik1
		return ik1
	else
		counter = 0
		for i in (a_s[k]+t_s[k]+n_s[k]):ik2
			counter=i+2
			if (pk[i+2] - pk[i]) .> -(TRAJ_TIME_STEP_HOURS * 8 * 2)
				if (pk[i+2] - pk[i]) .<= 0
					return i+2
				elseif (pk[i+1] - pk[i]) .<= 0
					return i+1
				else
					return i
				end
			end
		end
		if counter == ik1
			return counter
		end
	end
end

function get_var_tWCB_normed(vr::String,starts,stops)

	#" !! This function can only be executed once WCB_tau_vars_realtime.nc has been written !!" 
        files = TRAJ_FILES
        fl = TRAJ_PATH*"WCB_tau_vars_realtime.nc"
        t6 = T6
	#"Initialize array"
        v_tot = Vector{Matrix{Float64}}()

	#"Load data from WCB_tau_vars.nc"
	ds = Dataset(fl)
	vr_ar = Array(ds[vr])	

	v_f = fill(NaN,(length(t6),101))
	
	#"Loop over all trajectories. LinearInterpolation is used to remap all variable array"
	#"from WCB_tau_vars_realtime.nc onto the size 101."
	if vr != "rtime"
		for i in 1:length(t6)
			v = vr_ar[i,starts[i]:stops[i]]
			f_int = LinearInterpolation(1:length(v),v,extrapolation_bc=Line())
			stps = (length(v)-1)/100
			upto = length(v)
			v_f[i,1:end] = f_int(1:stps:upto)
		end
	elseif vr == "rtime"
		for i in 1:length(t6)
			v = vr_ar[starts[i]:stops[i]]
			f_int = LinearInterpolation(1:length(v),v,extrapolation_bc=Line())
			stps = (length(v)-1)/100
			upto = length(v)
			v_f[i,1:end] = f_int(1:stps:upto)
		end
	end
	v_tot = vcat(v_tot,[v_f[i,:] for i in 1:size(v_f)[1]])	
	close(ds)
        v_out = reduce(vcat,transpose.(v_tot))
	if vr in MICRO_VAR_LIST
		v_out = v_out .- v_out[:,1]
	end
	return v_out
end

#check if 1d array (including NaNs) is monotone
function is_monotone(ar)
    inc = [ar[i+1] .- ar[i] for i in 1:length(ar)-1]
    if nm.minimum(inc) >= 0
        return true
    else
        return false
    end
end

#how many restarts? this takes an array that is like [nm.sum(ds["evap"][:,i]) for i in 1:length(r_t)]
function where_restarts(ar)
    inc = [ar[i+1] .- ar[i] for i in 1:length(ar)-1]
    return findall(inc .< 0)
end


function get_var_realtime(vr::String)
	#"This function retrieves a variable for all WCB trajectories by going through all"
	#"files and using the WCB metadata to select the WCB trajectories and extracting the"
	#"variables FOR THE ENTIRE SIMULATION. Returns them in simulation realtime and is NOT"
	#"reorganized in time. If a trajectory begins after the first rtime step, then the"
	#"variable is NaN until it starts."

        files = TRAJ_FILES
        pth = TRAJ_PATH

	#"Initialize array."
        v_tot = Vector{Matrix{Float64}}()
		
	#"Loop over all trajectory files."
        for flnm in 1:length(WCB_fls)
		f = files[flnm]
                #println(f)
		#"for the nested sim the first four time steps in the first couple of 
                #files are weird, thus they are skipped."
			
		#"Get data from trajectory file"
                ds = Dataset(f)

		#"Get metadata"
                w,a,t6 = WCB_info(WCB_fls[flnm])

		#"Initialize array"
                v_f = fill(NaN,(length(w),NUMBER_RTIME_STEPS))

		#"Convert rtime to hours"
		r_t = ds["rtime"][:] ./ (3600)
        
                frm1 = NUMBER_RTIME_STEPS - length(r_t) + 1

                micro_var = ds["cond"][:,:]
                for i in 1:size(micro_var)[1]
                    ix = findfirst((!).(isnan).(micro_var[i,:]))
                    micro_var[i,:] = micro_var[i,:] .- micro_var[i,ix]
                end

                micro_summed = [nm.sum(micro_var[:,i]) for i in 1:size(micro_var)[2]]
                #is this array monotone? If not, there is a restart
                monotone = is_monotone(micro_summed)
                if monotone == false
                    restart_inds = where_restarts(micro_summed)
                    restart = true
                elseif monotone == true
                    restart = false
                    restart_inds = false
                end

		#"Here variables are loaded in a loop that is multithreaded."
		#"If variable is NOT microphysical, take the ascent part and write to v_f array."
		#"If it is microphysical, it needs to go through the funciton "micro_singleline"."
		#"That function makes sure the variable is set to zero at beginning of ascent and"
		#"if fixes problems associated with ICON restarts. See the function for more info."
		if vr != "rtime"
                        vr_ar = ds[vr][:,:]
                       	print(size(vr_ar))
			println(f)
			for i in 1:length(w)

		 	    if vr ∉ MICRO_VAR_LIST
                            	v_f[i,frm1:end] = vr_ar[w[i],:]
			    elseif vr in MICRO_VAR_LIST
				v_f[i,frm1:end] = micro_singleline(vr_ar[w[i],:],restart,restart_inds)
			    end	
			end
		elseif vr == "rtime"
			error("rtime isn't necessary in realtime")
		end
		v_tot = vcat(v_tot,[v_f[i,:] for i in 1:size(v_f)[1]])	
                close(ds)
        end
        GC.gc()
        return reduce(vcat,transpose.(v_tot))
end


function micro_singleline(vr_ar,restart,restart_inds)
	#"This function takes a 1-dim rtime and variable array and makes it such that"
	#"the variable resets for every restart."
	#"This is for the microphysical variables which we want to start at zero at"
	#"the beginning of the ascent and accumulate with time."
	#"If the ICON simulation has to be restarted at a certain rtime, then the variable"
	#"will reset to zero in the simultion. To compensate for this, this function brings"
	#"the variable to the value it had before the restart."
        if restart == false
            return vr_ar .- vr_ar[1]
        end
        out = copy(vr_ar)
        for ix in reverse(restart_inds)
                out[ix+1:end] .= out[ix+1:end] .+ out[ix]
        end
        #make first time step 0
        out = out .- out[1]
        return out
end
