using Glob, NCDatasets, Plots
using Plots.PlotMeasures
using Statistics
using BenchmarkTools, Dates, DataStructures
using DelimitedFiles
using NaNMath; nm=NaNMath
using StatsBase
using LaTeXStrings
using ColorSchemeTools

include("NAMELIST.jl")
include(CODE_PATH*"WCB_selection.jl")
include(CODE_PATH*"thermodynamic_functions.jl")
include(CODE_PATH*"statistical_functions.jl")


#---------------------------------------------------------------------------
#       "Skript for general work with ICON WCB trajectory files etc"
#       "Created by Cornelis Schwenk 27.11.2023"
#       "c.schwenk@uni-mainz.de"

#	"This skript should only be included once all of the WCB_vars... files"
#	"have been created. See WCB_selection.jl"
#---------------------------------------------------------------------------


#	"define some files and arrays we can use generally"
#	"paths and so on are defined in WCB_selection.jl"

#----------glob NWP-------------------
NWP_fls = glob("*mod_NWP*.nc",NWP_PATH)

#glob_ds_NWP = Dataset(glob_NWP_fls[1])
#glob_NWP_lon = glob_ds_NWP["lon"][:]
#glob_NWP_lat = glob_ds_NWP["lat"][:]

#----------nest NWP-------------------
nest_NWP_fls_dom1 = glob("*mod_NWP*DOM01*.nc",NWP_PATH)
nest_NWP_fls_dom2 = glob("*mod_NWP*DOM02*.nc",NWP_PATH)
nest_NWP_fls_dom3 = glob("*mod_NWP*DOM03*.nc",NWP_PATH)

nest_dom1_ds_NWP = Dataset(nest_NWP_fls_dom1[1])
nest_NWP_lon_dom1 = nest_dom1_ds_NWP["lon"][:]
nest_NWP_lat_dom1 = nest_dom1_ds_NWP["lat"][:]

#nest_dom2_ds_NWP = Dataset(nest_NWP_fls_dom2[1])
#nest_NWP_lon_dom2 = nest_dom2_ds_NWP["lon"][:]
#nest_NWP_lat_dom2 = nest_dom2_ds_NWP["lat"][:]

nest_dom3_ds_NWP = Dataset(nest_NWP_fls_dom3[1])
nest_NWP_lon_dom3 = nest_dom3_ds_NWP["lon"][:]
nest_NWP_lat_dom3 = nest_dom3_ds_NWP["lat"][:]

#------------WCB-----------------------

#	"Domain borders"
maxlon1 = maximum(nest_NWP_lon_dom1)
minlon1 = minimum(nest_NWP_lon_dom1)
maxlat1 = maximum(nest_NWP_lat_dom1)
minlat1 = minimum(nest_NWP_lat_dom1)

#maxlon2 = maximum(nest_NWP_lon_dom2)
#minlon2 = minimum(nest_NWP_lon_dom2)
#maxlat2 = maximum(nest_NWP_lat_dom2)
#minlat2 = minimum(nest_NWP_lat_dom2)

maxlon3 = maximum(nest_NWP_lon_dom3)
minlon3 = minimum(nest_NWP_lon_dom3)
maxlat3 = maximum(nest_NWP_lat_dom3)
minlat3 = minimum(nest_NWP_lat_dom3)

close(nest_dom3_ds_NWP)
close(nest_dom1_ds_NWP)

global bin_ar = collect(240:5:845)
global Tbin_ar = collect(-50:2:20)
global hymet = ["c","r","i","s","g","h"]


xar_norm = LinRange(0,1,101)
dsd = Dict()
dsdn = Dict()
dsdr = Dict()
t6_d = Dict()

if length(glob("WCB_tau_vars_normed.nc",TRAJ_PATH)) > 0
	global dsn = Dataset(TRAJ_PATH*"WCB_tau_vars_normed.nc")
        global TAU_600 = dsn["tau_600"][:]
        global TAU_WCB = dsn["tau_WCB"][:]
end

#	"Math label used for plotting"
taulabel = L"\tau_{600}\quad [\mathrm{h}]"


function find_max_dp2()
	#"This function finds the fastest 2hr pressure ascent"
	#"during the WCB-ascent"
        p_ar = PR
        ain = AINS	
	dp2_out = zeros(1:size(p_ar)[1])
	for k in 1:length(dp2_out)
		p_k = p_ar[k,ain[k][2]:end]
		dp2_1 = maximum(p_k[1:end-4] .- p_k[5:end])
		dp2_2 = maximum(p_k[2:end-4] .- p_k[6:end])
		dp2_3 = maximum(p_k[3:end-4] .- p_k[7:end])
		dp2_4 = maximum(p_k[4:end-4] .- p_k[8:end])
		dp2_out[k] = maximum([dp2_1,dp2_2,dp2_3,dp2_4])
	end
	return dp2_out
end

function find_max_dp2_op()
	#"This function finds the fastest 2hr pressure DESCENT"
	#"during the WCB-ascent"
        p_ar = PR
        ain = AINS	
	
	dp2_out = zeros(1:size(p_ar)[1])
	for k in 1:length(dp2_out)
		p_k = p_ar[k,ain[k][2]:end]
		dp2_1 = minimum(p_k[1:end-4] .- p_k[5:end])
		dp2_2 = minimum(p_k[2:end-4] .- p_k[6:end])
		dp2_3 = minimum(p_k[3:end-4] .- p_k[7:end])
		dp2_4 = minimum(p_k[4:end-4] .- p_k[8:end])
		dp2_out[k] = minimum([dp2_1,dp2_2,dp2_3,dp2_4])
	end
	return dp2_out
end

function find_dp2(rt::Int)	
        p_ar = PR
        if rt > (NUMBER_RTIME_STEPS - 4)
		return [NaN for i in 1:size(p_ar)[1]]
	end

	return p_ar[:,rt+3] .- p_ar[:,rt]	
end

function find_dp2_eff(p_ar,rt)
        if rt > (NUMBER_RTIME_STEPS - 4)
		return [NaN for i in 1:size(p_ar)[1]]
	end
	return p_ar[:,rt+3] .- p_ar[:,rt]
end

function get_dp2_all()
	dp2all = zeros(size(PR))
        for i in 1:NUMBER_RTIME_STEPS
                dp2all[:,i] = find_dp2_eff(PR,i)
        end
	return dp2all
end

function find_CAP()
	#"This function finds the fastest 2hr pressure ascent"
	#"during the WCB-ascent"
        p_ar = PR
        ain = AINS	
	
	CAP_out = zeros(1:size(p_ar)[1])
	for k in 1:length(CAP_out)
		p_k = p_ar[k,ain[k][2]:end]
		CAP_out[k] = find_max_dpx(p_k,1)
	end
	return CAP_out
end

if !@isdefined dp2
	global DP2 = find_max_dp2()
end

global M_c = TAU_600 .< 5
global M_s = (TAU_600 .> 20) .& (DP2 .< 350)
#global M_n = 18 .>= t6n .>= 6

#---- Cloud parameters-----
global a_geo_cloud = 0.124
global b_geo_cloud = 1/3
global a_geo_rain = 0.124
global b_geo_rain = 1/3
global a_geo_ice = 0.835
global b_geo_ice = 0.390
global a_geo_snow = 2.4
global b_geo_snow = 0.455
global a_geo_graupel = 0.142
global b_geo_graupel = 0.314
global a_geo_hail = 0.1366
global b_geo_hail = 1/3

#------------------------------------------------------------
#		Functions
#------------------------------------------------------------

#	"Conversion of rtime to DateTime and vice versa"
function calc_D_ice(q,qn)
	m = q / qn
	return a_geo_ice * (m ^ b_geo_ice)
end

function calc_radius_ice(q,qn)
	m = q ./ qn
        m[qn .< 10] .= NaN
	return 0.5 .* a_geo_ice .* (m .^ b_geo_ice) 
end

function calc_radius_hyd(q,qn,hyd::String)
        if (hyd == "c") | (hyd == "cloud")
            a_geo = a_geo_cloud
            b_geo = b_geo_cloud
        elseif (hyd == "r") | (hyd == "rain")
            a_geo = a_geo_rain
            b_geo = b_geo_rain
        elseif (hyd == "i") | (hyd == "ice")
            a_geo = a_geo_ice
            b_geo = b_geo_ice
        elseif (hyd == "s") | (hyd == "snow")
            a_geo = a_geo_snow
            b_geo = b_geo_snow
        elseif (hyd == "g") | (hyd == "graupel")
            a_geo = a_geo_graupel
            b_geo = b_geo_graupel
        elseif (hyd == "h") | (hyd == "hail")
            a_geo = a_geo_hail
            b_geo = b_geo_hail
        end
        m = q ./ qn
        m[qn .< 10] .= NaN
        m[q .<= 0] .= NaN
        r = 0.5 .* a_geo .* (m .^ b_geo)
        return r
end

function calc_fx_ice(q,qn)
	x = q ./ qn
	return 10 .* qn .* exp( 60 .^ (1/3) ) ./ x 
end

function calc_vertical_velocity()
	#"returns ascent velocity in m/s derived from altitude change over runtime"
	ds = dsn
	alt = Array(ds["alt"])
	rt = Array(ds["rtime"])
	altd = alt[:,2:end] .- alt[:,1:end-1]
	rtd = rt[:,2:end] .- rt[:,1:end-1]
	wv = altd ./ (rtd .* 60 .* 60)
	return wv
end

function calc_vertical_velocity_realtime()
        #"returns ascent velocity in m/s derived from altitude change over runtime"
        ds = dsr
        alt = Array(ds["alt"])
        altd = alt[:,2:end] .- alt[:,1:end-1]
        rtd = 0.5
        wv = altd ./ (rtd .* 60 .* 60)
        return wv
end


function t6_binning(ar,t6,func)
	#"Bin func(values) into tau_600 bins."
	#"Used to calculate median and mean as function of t6."
	#"Example: CR_mean = t6_binning(CR,tau_600,mean)."
	#"Uses bin_two_ars function from statistical_functions.jl"
	bins = 1:48
	return bin_two_ars(ar,t6,bins,func)
end

function t6_binning_fine(ar,t6,func)
        v = []
        for t in collect(1:0.5:48)
                if length(ar[t-1 .< t6 .<= t]) == 0
                        push!(v,NaN)
                else
                        push!(v,func(ar[t-1 .< t6 .<= t]))
                end
        end
        return v
end

function t6_prcntl(ar,t6,prc)
	v = Float64[]
	for t in collect(1:0.5:48)
		if length(ar[t-1 .< t6 .<= t]) == 0
			push!(v,NaN)
		else
			s_ar = sort(ar[t-1 .< t6 .<= t])
			iprc = Int(round(length(s_ar)*prc/100,digits=0))
			push!(v,s_ar[iprc])
		end
	end
	return v
end

function anyar_prcntl(ar1,ar2,bins,prc)
	#ar1 is binned into bins of ar2
	v = Float64[]
	for i in 1:(length(bins)-1)
		if length(ar1[bins[i] .< ar2 .<= bins[i+1]]) < 100
			push!(v,NaN)
		else
			s_ar = sort(ar1[bins[i] .< ar2 .<= bins[i+1]])
			push!(v,prcntl_1d(s_ar,prc))
		end
	end
	return v
end


function arnorm(ar)
	return (ar .- mean(ar)) ./ std(ar)
end

function taumask(a1,a2)
	#"Create a mask for a1 .< tau_600 .< a2"
        tau = dsn["tau_600"][:] ./ 2
        return (tau .> a1) .&& (tau .<= a2)
end


function find_min_dpx(p,x)
	#"takes 1-dim p_array and time steps x"
	#"finds the minimum ascent in x time steps"
	dpx = p[1:end-x] .- p[x+1:end]
	return minimum(dpx)
end

function find_max_dpx(p,x)
	#"takes 1-dim p_array and time steps x"
	#"finds the maximum ascent in x time steps"
	dpx = p[1:end-x] .- p[x+1:end]
	return maximum(dpx)
end

function flux(hy::String,form="normed",masking=false,k1=1,k2=6)
	#"Calculate precipitation flux of specific species"
	#"Chose normed, regular or realtime to specify which WCB_tau_vars.... to use"
	if form == "normed"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_normed.nc")
	elseif form =="realtime"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_realtime.nc")
	else
		error("form must be 'realtime', or 'normed'")
	end
	ins = Array(ds["q$(hy)_in"])
	out = Array(ds["q$(hy)_out"])

	if masking==true
		ins = in[taumask(k1,k2),:]
		out = out[taumask(k1,k2),:]
	end
	return ins .- out
end

function flux_fast_ind(hy::String,form="normed",ind=101,masking=false,k1=1,k2=6)
	#"Faster version for specific time step"
	if form == "normed"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_normed.nc")
	elseif form =="realtime"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_realtime.nc")
	else
		error("form must be 'realtime', or 'normed'")
	end
	ins = ds["q$(hy)_in"][:,ind]
	out = ds["q$(hy)_out"][:,ind]

	if masking==true
		ins = in[taumask(k1,k2)]
		out = out[taumask(k1,k2)]
	end
	return ins .- out
end

function flux_superquick_ind(hy::String,ds,ind::Int)
	#"Faster version for specific time step, input ds"
	ins = ds["q$(hy)_in"][:,ind]
	out = ds["q$(hy)_out"][:,ind]
	return ins .- out
end

#	"Faster version for just the last time step"
function flux_fast_last(hy::String,form::String,masking=false,k1=1,k2=6)
	if form == "normed"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_normed.nc")
	elseif form =="realtime"
		ds = Dataset(TRAJ_PATH*"WCB_tau_vars_realtime.nc")
	else
		error("form must be 'realtime', or 'normed'")
	end
	ins = ds["q$(hy)_in"][:,end]
	out = ds["q$(hy)_out"][:,end]

	if masking==true
		ins = in[taumask(k1,k2)]
		out = out[taumask(k1,k2)]
	end
	return ins .- out
end

function calc_flux_realtime(hy::String,ts::Int)
	#"This calculates the flux at ts timesteps after end of ascent from the"
	#"realtime file and makes sure at beginning of ascent, the flux is zero."
	ins = get_realtime_var_t6_plus_ts("q$(hy)_in",ts,true)
	outs = get_realtime_var_t6_plus_ts("q$(hy)_out",ts,true)
	Pts = ins .- outs
	return Pts
end

function calc_tot_flux_realtime(ts::Int)
	#"This calculates the flux at ts timesteps after end of ascent from the"
	#"realtime file and makes sure at beginning of ascent, the flux is zero."
	#"This is TOTAL flux"
	ins = sum([get_realtime_var_t6_plus_ts("q$(hy)_in",ts,true) for hy in hymet[2:end]])
	outs = sum([get_realtime_var_t6_plus_ts("q$(hy)_out",ts,true) for hy in hymet[2:end]])
	Pts = ins .- outs
	return Pts
end

function get_t6_plus_ts_indices(ts::Int)
	#"This function gives the CartesianIndices of trajectories in the"
	#"realtime file ts-timestepts after the end of tau_600 ascent. The msk"
	#"is a mask that is false if ts-hours is outside of the sim."
	ds = dsr
        t6_p_ts = A6 .+ NS .+ T6 .+ ts
        upto = length(T6)

	ts_inds = [CartesianIndex(i,t6_p_ts[i]) for i in 1:upto]
	msk = t6_p_ts .<= NUMBER_RTIME_STEPS
	return (ts_inds,msk)
end

function get_tWCB_plus_ts_indices(ts::Int)
	#"This function gives the CartesianIndices of trajectories in the"
	#"realtime file ts-timestepts after the end of WCB ascent. The msk" 
	#"is a mask that is false if ts-hours is outside of the sim."
	ds = dsr
	upto = length(dsn["tau_600"])
	ts_inds = [CartesianIndex(k,find_stop(k) + ts) for k in 1:upto]
	msk = [ts_inds[k][2] for k in 1:upto] .<= NUMBER_RTIME_STEPS
	return (ts_inds,msk)
end

function get_realtime_var_t6_plus_ts(vr::String,ts::Int,micro=false)
	#"This function gets a variable ts-timesteps after the ascent"
	ts_inds,msk=get_t6_plus_ts_indices(ts)
	vrar = Array(dsr[vr])
	vrout = vrar[ts_inds[msk]]
	if micro==true
		#"if micro=true, then we need to subtract the value at the beginning"
		#"of the ascent. As per definition of our microphysical variables"
		vrst = vrar[AINS[msk]]
		vrout = vrout .- vrst
	end
	return vrout
end

#------------------------------------------------
#	"Calculations for water budget and residuals etc."
#	"See publication for details."
#	"The following functions only use the time-normalized data from WCB_tau_vars_normed.nc"
#------------------------------------------------

function calc_P()
        return sum([flux(hy,"normed") for hy in hymet[2:end]])
end

function calc_dqv_traj(ind=101::Int)
        #"calculate change in moisture along ascent"
        ds = dsn
        dqv_traj = -ds["qihh"][:,ind] .- ds["qc_nuc"][:,ind] .- ds["cond"][:,ind] .-
                                ds["evap"][:,ind] .+ ds["r_evap"][:,ind] .+ ds["f_evap"][:,ind] .-
                                ds["qx_dep"][:,ind] .- ds["satad2"][:,ind] .+ ds["qvturb"][:,ind] .+
                                ds["qvconv"][:,ind]
        return dqv_traj
end

function calc_dqv_traj_all(ds)
        #"calculate change in moisture along ascent"
        dqv_traj = -Array(ds["qihh"]) .- Array(ds["qc_nuc"]) .- Array(ds["cond"]) .-
                                Array(ds["evap"]) .+ Array(ds["r_evap"]) .+ Array(ds["f_evap"]) .-
                                Array(ds["qx_dep"]) .- Array(ds["satad2"]) .+ Array(ds["qvturb"]) .+
                                Array(ds["qvconv"])
        return dqv_traj
end


function calc_dqv_traj_realtime(ts::Int)
        #"calculate change in vapor along ascent"
        ds = dsr
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	vrlist1 = ["qihh","qc_nuc","cond","evap","qx_dep","satad2"]
	vrlist2 = ["r_evap","f_evap","qvturb","qvconv"]
	dqhy_traj2 = sum([get_realtime_var_t6_plus_ts(vr,ts,true) for vr in vrlist2])
	dqhy_traj1 = sum([get_realtime_var_t6_plus_ts(vr,ts,true) for vr in vrlist1])
        return -dqhy_traj1 .+ dqhy_traj2
end

function calc_dqhy_traj(ind=101::Int)
        #"calculate change in hydrometeors along ascent"
        ds = dsn
        allfl = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
        dqhy_traj = ds["qihh"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["cond"][:,ind] .+
                        ds["evap"][:,ind] .-
                        ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["qx_dep"][:,ind] .+
                        ds["satad2"][:,ind] .+ ds["qcturb"][:,ind] .+ ds["qcconv"][:,ind] .+
                        ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .+ ds["qrconv"][:,ind] .- allfl;
        return dqhy_traj
end

function calc_dqhy_traj_superquick(ds,ind::Int)
        #"calculate change in hydrometeors along ascent"
        allfl = sum([flux_superquick_ind(hy,ds,ind) for hy in hymet[2:end]])
        dqhy_traj = ds["qihh"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["cond"][:,ind] .+
                        ds["evap"][:,ind] .-
                        ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["qx_dep"][:,ind] .+
                        ds["satad2"][:,ind] .+ ds["qcturb"][:,ind] .+ ds["qcconv"][:,ind] .+
                        ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .+ ds["qrconv"][:,ind] .- allfl;
        return dqhy_traj
end

function calc_dqhy_traj_all(ds,allfl)
        #"calculate change in hydrometeors along ascent"
	dqhy_traj = ds["qihh"][:,:] .+ ds["qc_nuc"][:,:] .+ ds["cond"][:,:] .+ ds["evap"][:,:] .-
                        ds["r_evap"][:,:] .- ds["f_evap"][:,:] .+ ds["qx_dep"][:,:] .+
                        ds["satad2"][:,:] .+ ds["qcturb"][:,:] .+ ds["qcconv"][:,:] .+
                        ds["qsconv"][:,:] .+ ds["qiconv"][:,:] .+ ds["qrconv"][:,:] .- allfl;
        return dqhy_traj
end

function calc_dqhy_traj_realtime(ts::Int,Pts)
        #"calculate change in hydrometeors along ascent"
	#"for efficiency, calculate total precip at ts (Pts) beforehand and insert as arg"
        ds = dsr
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	vrlist1 = ["qihh","qc_nuc","cond","evap","qx_dep","satad2","qcturb","qcconv",
		  "qsconv","qiconv","qrconv"]
	vrlist2 = ["r_evap","f_evap"]
	dqhy_traj2 = sum([get_realtime_var_t6_plus_ts(vr,ts,true) for vr in vrlist2])
	dqhy_traj1 = sum([get_realtime_var_t6_plus_ts(vr,ts,true) for vr in vrlist1])
        return dqhy_traj1 .- dqhy_traj2 .- Pts
end

function calc_res_qv(ind=101::Int)
        #"calculate residual water vapor along ascent"
        ds = dsn
        res_qv = ds["qv"][:,1] .- ds["qv"][:,ind] .+ calc_dqv_traj(ind)
        return res_qv
end

function calc_res_qv_all(ds)
        #"calculate residual water vapor along ascent"
        res_qv = ds["qv"][:,1] .- ds["qv"][:,:] .+ calc_dqv_traj_all(ds)
        return res_qv
end

function calc_res_qv_realtime(ts::Int)
        #"calculate residual vapor along ascent"
        ds = dsr
	ts_inds,msk = get_t6_plus_ts_indices(ts)
        qvts = ds["qv"][:,:][ts_inds[msk]]
	qv1 = ds["qv"][:,:][AINS]
        res_qv = qv1[msk] .- qvts .+ calc_dqv_traj_realtime(ts)
        return res_qv
end

function calc_res_qhy(ind=101::Int)
        #"calculate residual hydrometeors along ascent"
        ds = dsn
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        qhy0 = sum([ds["q$(hy)"][:,1] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj(ind)
        return res_qhy
end

function calc_res_qhy_quick(qhy0,ind=101::Int)
        #"calculate residual hydrometeors along ascent"
        ds = dsn
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj(ind)
        return res_qhy
end

function calc_res_qhy_superquick(ds,qhy0,ind::Int)
        #"calculate residual hydrometeors along ascent"
        hy_tot_ind = sum([ds["q$(hy)"][:,ind] for hy in hymet])
        res_qhy = qhy0 .- hy_tot_ind .+ calc_dqhy_traj_superquick(ds,ind)
        return res_qhy
end

function calc_res_qhy_all(ds,qhy,allfl)
        #"calculate residual hydrometeors along ascent"
        res_qhy = qhy[:,1] .- qhy .+ calc_dqhy_traj_all(ds,allfl)
        return res_qhy
end

function calc_res_qhy_realtime(ts::Int,Pts)
        #"calculate residual hydrometeors along ascent"
        ds = dsr
	ts_inds,msk = get_t6_plus_ts_indices(ts)
        hy_tot_ind = sum([ds["q$(hy)"][:,:][ts_inds[msk]] for hy in hymet])
        hy_tot_1 = sum([ds["q$(hy)"][:,:][AINS[msk]] for hy in hymet])
        res_qhy = hy_tot_1 .- hy_tot_ind .+ calc_dqhy_traj_realtime(ts,Pts)
        return res_qhy
end

function calc_res(ind=101::Int)
        return calc_res_qhy(ind) .+ calc_res_qv(ind)
end

function calc_res_all(ds,qhy,allfl)
        return calc_res_qhy_all(ds,qhy,allfl) .+ calc_res_qv_all(ds)
end

function calc_res_realtime(ts::Int,Pts)
        return calc_res_qhy_realtime(ts,Pts) .+ calc_res_qv_realtime(ts)
end

function calc_HYD(ind=101::Int)
        ds = dsn
        return sum([ds["q$(hy)"][:,1] for hy in hymet]) .+ ds["qcturc"][:,ind] .+
                        ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ ds["qsconv"][:,ind] .+
                        ds["qiconv"][:,ind] .- calc_res_qhy(ind)
end

function calc_HYD_quick(qhy0,ind=101)
	ds = dsn
	return qhy0 .+ ds["qcturc"][:,ind] .+ ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ 
		ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .- 
		calc_res_qhy_quick(qhy0,ind)
end

function calc_HYD_superquick(ds,qhy0,ind::Int)
	return qhy0 .+ ds["qcturc"][:,ind] .+ ds["qcconv"][:,ind] .+ ds["qrconv"][:,ind] .+ 
		ds["qsconv"][:,ind] .+ ds["qiconv"][:,ind] .- 
		calc_res_qhy_superquick(ds,qhy0,ind)
end

function calc_HYD_all(ds,qhy,allfl)
	one = qhy[:,1] .+ ds["qcturc"][:,:] .+ ds["qcconv"][:,:] .+ ds["qrconv"][:,:]
	two = ds["qsconv"][:,:] .+ ds["qiconv"][:,:]
	GC.gc()
	return one .+ two .- calc_res_qhy_all(ds,qhy,allfl)
end

function calc_HYD_realtime(ts::Int,Pts)
	#"for efficiency, Pts is calculated beforehand and entered as argument"
        ds = dsr
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	qcturc = get_realtime_var_t6_plus_ts("qcturc",ts,true)
	qcconv = get_realtime_var_t6_plus_ts("qcconv",ts,true)
	qrconv = get_realtime_var_t6_plus_ts("qrconv",ts,true)
	qsconv = get_realtime_var_t6_plus_ts("qsconv",ts,true)
	qiconv = get_realtime_var_t6_plus_ts("qiconv",ts,true)
	qhy0 = sum([ds["q$(hy)"][:,:][AINS[msk]] for hy in hymet])
	resqhy = calc_res_qhy_realtime(ts,Pts)
	return qhy0 .+ qcturc + qcconv .+ qrconv .+ qsconv .+ qiconv .- resqhy
end

function calc_VAP(ind=101::Int)
        ds = dsn
        return ds["qv"][:,1] .+ ds["qvturc"][:,ind] .+ ds["qvconv"][:,ind] .- calc_res_qv(ind)
end

function calc_VAP_all(ds)
        return ds["qv"][:,1] .+ ds["qvturc"][:,:] .+ ds["qvconv"][:,:] .- calc_res_qv_all(ds)
end

function calc_VAP_realtime(ts::Int)
	ts_inds,msk = get_t6_plus_ts_indices(ts)
        ds = dsr
	qvturc = get_realtime_var_t6_plus_ts("qvturc",ts,true)
	qvconv = get_realtime_var_t6_plus_ts("qvconv",ts,true)
        return ds["qv"][:,:][ains[msk]] .+ qvturc .+ qvconv .- calc_res_qv_realtime(ts)
end

function calc_qtot(ind=101::Int)
        ds = dsn
        return sum([ds["q$(hy)"][:,ind] for hy in ["v","c","r","i","s","g","h"]])
end

function calc_qhy(ind=101::Int)
        ds = dsn
        return sum([ds["q$(hy)"][:,ind] for hy in hymet])
end

function calc_qtot_all()
        ds = dsn
	qtot = ds["qv"][:,:]
	for hy in ["c","r","i","s","g","h"]
		qtot = qtot .+ ds["q$(hy)"][:,:]
	end
	GC.gc()
	return qtot
end

function calc_qhy_all()
        ds = dsn
	qhy = ds["qc"][:,:]
	for hy in ["r","i","s","g","h"]
		qhy = qhy .+ ds["q$(hy)"][:,:]
	end
	GC.gc()
	return qhy
end

function calc_qtot_realtime(ts::Int)
	ts_inds,msk = get_t6_plus_ts_indices(ts)
        ds = dsr
        return sum([ds["q$(hy)"][:,:][ts_inds[msk]] for hy in ["v","c","r","i","s","g","h"]])
end

function calc_convs(ind=101::Int)
        ds = dsn
        return sum([ds["q$(x)conv"][:,ind] for x in ["v","c","r","i","s"]])
end

function calc_convs_all(ds)
        convs = ds["qvconv"][:,:]
	for hy in ["c","r","i","s"]
		convs = convs .+ ds["q$(hy)conv"][:,:]
	end
	return convs
end

function calc_convs_realtime(ts::Int)
	qvconv = get_realtime_var_t6_plus_ts("qvconv",ts,true)
	qcconv = get_realtime_var_t6_plus_ts("qcconv",ts,true)
	qrconv = get_realtime_var_t6_plus_ts("qrconv",ts,true)
	qiconv = get_realtime_var_t6_plus_ts("qiconv",ts,true)
	qsconv = get_realtime_var_t6_plus_ts("qsconv",ts,true)
        return qvconv .+ qcconv .+ qrconv .+ qiconv .+ qsconv
end

function calc_turbs(ind=101::Int)
        ds = dsn
        return sum([ds["q$(hy)turb"][:,ind] for hy in ["c","v"]])
end

function calc_turbs_all(ds)
        return sum([ds["q$(hy)turb"][:,:] for hy in ["c","v"]])
end

function calc_turbs_realtime(ts::Int)
	qvturb = get_realtime_var_t6_plus_ts("qvturb",ts,true)
	qcturb = get_realtime_var_t6_plus_ts("qcturb",ts,true)
        return qvturb .+ qcturb
end

function calc_turcs(ind=101::Int)
        ds = dsn
        return sum([ds["q$(hy)turc"][:,ind] for hy in ["c","v"]])
end

function calc_turcs_realtime(ts::Int)
	qvturb = get_realtime_var_t6_plus_ts("qvturc",ts,true)
	qcturb = get_realtime_var_t6_plus_ts("qcturc",ts,true)
        return qvturc .+ qcturc
end

function calc_qtc(ind=101::Int)
        return calc_convs(ind) .+ calc_turbs(ind)
end

function calc_qtc_all(ds)
        return calc_convs_all(ds) .+ calc_turbs_all(ds)
end

function calc_qtc_realtime(ts::Int)
        return calc_convs_realtime(ts) .+ calc_turbs_realtime(ts)
end

function calc_qtcr(ind=101::Int)
        return calc_res(ind) .- calc_qtc(ind)
end

function calc_qtcr_all(ds,qhy,allfl)
        return calc_res_all(ds,qhy,allfl) .- calc_qtc_all(ds)
end

function calc_qtcr_realtime(ts::Int,Pts)
        return calc_res_realtime(ts,Pts) .- calc_qtc_realtime(ts)
end

function calc_C_hy(ind=101::Int)
        ds = dsn
        C_hy = ds["cond"][:,ind] .+ ds["qc_nuc"][:,ind] .+ ds["qihh"][:,ind] .+ ds["depo"][:,ind]
        return C_hy
end

function calc_C_hy_realtime(ts::Int)
	cond_ = get_realtime_var_t6_plus_ts("cond",ts,true)
	qc_nuc = get_realtime_var_t6_plus_ts("qc_nuc",ts,true)
	qihh = get_realtime_var_t6_plus_ts("qihh",ts,true)
	depo = get_realtime_var_t6_plus_ts("depo",ts,true)
	C_hy = cond_ .+ qc_nuc .+ qihh .+ depo	
        return C_hy
end

function calc_E_v(ind=101::Int)
        ds = dsn
        E_v = ds["evap"][:,ind] .- ds["r_evap"][:,ind] .- ds["f_evap"][:,ind] .+ ds["subl"][:,ind]
        return E_v
end

function calc_E_v_realtime(ts::Int)
	evap = get_realtime_var_t6_plus_ts("evap",ts,true)
	r_evap = get_realtime_var_t6_plus_ts("r_evap",ts,true)
	f_evap = get_realtime_var_t6_plus_ts("f_evap",ts,true)
	subl = get_realtime_var_t6_plus_ts("subl",ts,true)
        E_v = evap .- r_evap .- f_evap .+ subl
	return E_v
end

function calc_eps(ind=101::Int)
        return calc_HYD(ind) ./ calc_VAP(ind)
end

function calc_eps_realtime(ts::Int,Pts)
        return calc_HYD_realtime(ts,Pts) ./ calc_VAP_realtime(ts)
end

function calc_PE(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
	ds = dsn
        tot_prec = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        qhy0 = sum([ds["q$(hy)"][:,1] for hy in hymet])
	HYD = calc_HYD_quick(qhy0,ind)
        PE = tot_prec ./ (C_hy .+ E_v .+ HYD)
	PE[PE .> 1] .= NaN
	PE[PE .< 0] .= NaN
	return PE
end

function calc_PE_time(qhy,allfl)
	ds = dsn
	
	C_hy = ds["cond"][:,:] .+ ds["qc_nuc"][:,:] .+ ds["qihh"][:,:] .+ ds["depo"][:,:]
        E_v = ds["evap"][:,:] .- ds["r_evap"][:,:] .- ds["f_evap"][:,:] .+ ds["subl"][:,:] 
	
	GC.gc()
	HYD = calc_HYD_all(ds,qhy,allfl)

        PE = allfl ./ (C_hy .+ E_v .+ HYD)
        PE[PE .> 1] .= NaN
        PE[PE .< 0] .= NaN

	GC.gc()
        return PE
end

function calc_PE_realtime(ts::Int)
        #"This function calculates PE ts-timesteps (1ts = 0.5h) after finishing the ascent."
	#"This takes a while ~2min"
        C_hy = calc_C_hy_realtime(ts)
	E_v = calc_E_v_realtime(ts)
	GC.gc()
        Pts = calc_tot_flux_realtime(ts)
        GC.gc()
	HYD = calc_HYD_realtime(ts,Pts)
        GC.gc()
        return Pts ./ (C_hy .+ E_v .+ HYD)
end

function calc_PE_Eul(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        P1 = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        PE = P1 ./ (C_hy .+ E_v)
        return PE
end

function calc_PE_notc(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        tot_prec = sum([flux_fast_ind(i,"normed",ind) for i in hymet[2:end]])
        qhy0 = sum([dsn["q$(hy)"][:,1] for hy in hymet])
        return tot_prec ./ (C_hy .+ E_v .+ qhy0)
end

function calc_CR(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        VAP = calc_VAP(ind)
        CR = (C_hy .+ E_v) ./ VAP
	CR[CR .> 1] .= NaN
	CR[CR .< 0] .= NaN
	return CR
end

function calc_CR_time()
	ds = dsn
        C_hy = ds["cond"][:,:] .+ ds["qc_nuc"][:,:] .+ ds["qihh"][:,:] .+ ds["depo"][:,:]
        E_v = ds["evap"][:,:] .- ds["r_evap"][:,:] .- ds["f_evap"][:,:] .+ ds["subl"][:,:] 
        VAP = calc_VAP_all(ds)
        CR = (C_hy .+ E_v) ./ VAP
	CR[CR .> 1] .= NaN
	CR[CR .< 0] .= NaN
	return CR
end

function calc_CR_Eul(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        Qtot0 = calc_qtot(1)
        CR = (C_hy .+ E_v) ./ Qtot0
        return CR
end

function calc_CR_notc(ind=101::Int)
        C_hy = calc_C_hy(ind)
        E_v = calc_E_v(ind)
        qv0 = dsn["qv"][:,1]
        return (C_hy .+ E_v) ./ qv0
end


function calc_DR(ind=101::Int)
        qtot_ind = calc_qtot(ind)
        qtot_1 = calc_qtot(1)
        return (qtot_1 .- qtot_ind) ./ (qtot_1)
end

function calc_DR_time()
        qtot_t = mapreduce(permutedims,vcat,[calc_qtot(i) for i in 1:101])'
        qtot_1 = calc_qtot(1)
        return (qtot_1 .- qtot_t) ./ (qtot_1)
end

function calc_DR_realtime(ts::Int)
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	qtotr = sum([dsr["q$(hy)"][:,:][ts_inds[msk]] for hy in ["v","c","r","i","s","g","h"]])
	qtot0 = sum([dsr["q$(hy)"][:,:][AINS] for hy in ["v","c","r","i","s","g","h"]])
	return (qtot0[msk] .- qtotr) ./ qtot0[msk]
end

function calc_DR_mix(ind=101::Int)
        qtot0 = calc_qtot(1)
        qtcr = calc_qtcr(ind)
        return qtcr ./ qtot0
end

function calc_DR_mix_time(qhy,allfl)
	ds = dsn
        qtot0 = calc_qtot(1)
	qtcr = calc_qtcr_all(ds,qhy,allfl)
        return qtcr ./ qtot0
end

function calc_DR_mix_realtime(ts::Int)
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	qtot0 = sum([dsr["q$(hy)"][:,:][ains] for hy in ["v","c","r","i","s","g","h"]]) 
        Pts = calc_tot_flux_realtime(ts)
        GC.gc()
	qtcr = calc_qtcr_realtime(ts,Pts)
        return qtcr ./ qtot0[msk]
end

function calc_DR_mphys(ind=101::Int)
        P1 = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
        qtot0 = calc_qtot(1)
        qtcr = calc_qtcr(ind)
	DR_mphys =  P1 ./ (qtot0 .- qtcr)
	DR_mphys[DR_mphys .> 1] .= NaN
	DR_mphys[DR_mphys .< 0] .= NaN
	return DR_mphys
end

function calc_DR_mphys_time(qhy,allfl)
	ds = dsn
	#ins = sum([ds["q$(hy)_in"][:,:] for hy in hymet[2:end]])
        #outs = sum([ds["q$(hy)_out"][:,:] for hy in hymet[2:end]])
	#allfl = ins .- outs
	#qhy = sum([dsn["q$(hy)"][:,:] for hy in hymet])
	qtot0 = calc_qtot(1)
        #GC.gc()
	qtcr = calc_qtcr_all(ds,qhy,allfl)
	DR_mphys =  allfl ./ (qtot0 .- qtcr)
	DR_mphys[DR_mphys .> 1] .= NaN
	DR_mphys[DR_mphys .< 0] .= NaN
	return DR_mphys
end

function calc_DR_mphys_realtime(ts::Int)
        Pts = calc_tot_flux_realtime(ts)
	GC.gc()
	ts_inds,msk = get_t6_plus_ts_indices(ts)
	qtot0 = sum([dsr["q$(hy)"][:,:][AINS] for hy in ["v","c","r","i","s","g","h"]])
        qtcr = calc_qtcr_realtime(ts,Pts)
	return Pts ./ (qtot0[msk] .- qtcr)
end

function calc_DR_res(ind=101::Int)
        resqv = calc_res_qv(ind)
        resqhy = calc_res_qhy(ind)
        qtot = calc_qtot(1)
        return (resqv .+ resqhy) ./ qtot
end

function calc_q_removal(ind=101::Int)
        #"this function tracks the changes in qtot budget"
        #"the last time-entry is equal to change in qtot"
        ds = dsn
        allfl = sum([flux_fast_ind(hy,"normed",ind) for hy in hymet[2:end]])
        turbs = sum([ds["q$(hy)turb"][:,ind] for hy in ["c","v"]])
        convs = sum([ds["q$(x)conv"][:,ind] for x in ["v","c","r","i","s"]])
        resqv = calc_res_qv(ind)
        resqhy = calc_res_qhy(ind)
        return allfl .- turbs .- convs .+ resqv .+ resqhy
end

function calc_delta_H()
        ds = dsn
        return ds["qihh"][:,:] .+ ds["qc_nuc"][:,:] .+ ds["cond"][:,:] .+ ds["evap"][:,:] .+ ds["r_evap"][:,:] .+
                ds["f_evap"][:,:] .+ ds["qx_dep"][:,:] .+ ds["satad2"][:,:]
end
