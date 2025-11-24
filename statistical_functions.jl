using NaNMath; nm=NaNMath
using StatsBase, Glob, NCDatasets, Plots

#---------------------------------------------------------------------------
#	"Skript for some useful statistical manipulations of data."
#	"Created by Cornelis Schwenk 19.10.2023"
#	"c.schwenk@uni-mainz.de"
#---------------------------------------------------------------------------

function prcntl(ar,n::Int)
	#------------------------------
        #"This function gives you the n-th percentile"
	#"of a 2-dim array along the first dimension"
	#"This can handle NaNs"
	#------------------------------
       
	out = similar(ar, size(ar, 2)) 
	rows = size(ar,1)
		
        for j in 1:size(ar,2)
		col = view(ar,:,j)
		col_sorted = sort(col[.!isnan.(col)])
		if isempty(col_sorted)
			out[j] = NaN
		else
			i = Int(round(length(col_sorted) * n / 100, digits=0))
			out[j] = i > 0 ? col_sorted[i] : NaN
		end
	end
	return out
end

function prcntl_1d(ar,n::Int)
	#------------------------------
        #"This function gives you the n-th percentile"
	#"of a 1-dim array"
	#"This can handle NaNs"
	#------------------------------ 
	ar = sort(ar)
        a = ar[(!isnan).(ar)]
       	if length(a) == 0
        	return NaN
        else
       		i = Int(round(length(a)*n/100,digits=0))
		if i == 0
                    i = 1
                end
                return a[i]
	end
end

function anyar_prcntl(ar1,ar2,bins,prc)
	#------------------------------
        #"This function gives you the n-th percentile"
	#"of a 1-dim array ar1 binned into bins of ar2."
	#"This can handle NaNs"
	#------------------------------
        v = Float64[]
        for i in 1:(length(bins)-1)
                if length(ar1[bins[i] .< ar2 .<= bins[i+1]]) == 0
                        push!(v,NaN)
                else
                        s_ar = sort(ar1[bins[i] .< ar2 .<= bins[i+1]])
                        push!(v,prcntl_1d(s_ar,prc))
                end
        end
        return v
end

function pbin_array(vr_ar,p_ar,bins)
	#------------------------------
        #"This function bins values from vr_ar into pressure bins of p_ar."
	#"vr_ar[i,j] must also correspond to p_ar[i,j]."
	#"This returns a dictionary with each entry having all vr_ar points"
	#"that fall into that pressure bin"
	#"if you want the mean for example then take the output as P_d and:"
	#	"mean_ar = [mean(P_d[k]) for k in sort(keys(P_d))]"
	#------------------------------
        P_d = Dict()
        for pr in bins
                inds = findall((p .>= pr) .&& (p .<= pr+5))
                P_d[pr] = y[inds]
        end
        return P_d
end

function prcntl_of_pbins(bindict,prc)
	#-----------------------------
	#"This function gives the percentiles like prcntl but instead"
	#"of a dictionary created by pbin_array"
	#"Example:"
	#	"temp_d = pbin_array(temp_data,pressure_data)"
	#	"temp_10th_prcntl = prcntl_of_pbins(temp_d,10)"
	#-----------------------------
        x_ar = sort([k for k in keys(bindict)])
        prc_ar = []
        for x in x_ar
                l = length(bindict[x])/100
                ind = Int(round(l*prc,digits=0))
                push!(prc_ar,sort(bindict[x])[ind])
        end
        return prc_ar
end

function cumsum(x,bs::String)
	#-----------------------------
	#"This function takes the cumulative sum of a 1-dim array, but..."
	#"only of the positive or negative values."
	#"Example: cumsum([-2,-1,1,2,3,4],"big") = 1+2+3+4"
	#"Example: cumsum([-2,-1,1,2,3,4],"small") = -2 + -1"
	#"Its like a selective integration of only the positive or negative values"
	#-----------------------------
	k = 0
        v = []
        push!(v,k)
        y = copy(x)
        if bs == "big"
                y[y .< 0] .= 0
        elseif bs == "small"
                y[y .> 0] .= 0
        else
                println("error, bs must be big or small")
        end
        for i in 1:length(y)
                k = k + y[i]
                push!(v,k)
        end
        return v
end


function bin_two_ars(ar1,ar2,bins,func)
        #"This funciton bins func(ar1) into bins of ar2"
        #"used to calculate median and mean of ar1 as function of ar2"
        #"example: y_mean = t6_binning(ydata,xdata,xbins,mean)"
        v = []
        for t in 1:(length(bins)-1)
                if length(ar1[bins[t] .< ar2 .<= bins[t+1]]) == 0
                        push!(v,NaN)
                else
                        push!(v,func(ar1[bins[t] .< ar2 .<= bins[t+1]]))
                end
        end
        return v
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

function derv_2dar(x) 
        #"take derivative of 2dimensional array"
        return x[:,2:end] .- x[:,1:end-1]
end

function derv_ar(x)
        #"take derivative of 1dimensional array"
        return x[2:end] .- x[1:end-1]
end
function new_var(ar,bs::String)
        #"Return derivative of 2-dimensional array"
        return Matrix{Float64}(mapreduce(permutedims, vcat, map(x->cumsum(x,bs),eachrow(derv_2dar(ar)))))
end

function arnorm(ar)
        return (ar .- mean(ar)) ./ std(ar)
end

