include("/lustre/project/nhr-tpchange/gziyan/WCB_julia_code/WCB_selection.jl")

println("writing meta files")
println("---------------")
write_all_WCB_meta_files()

include("/lustre/project/nhr-tpchange/gziyan/WCB_julia_code/WCB_selection.jl")

print("writing realtime file")
println("---------------")
write_WCB_var_file_realtime()
include("/lustre/project/nhr-tpchange/gziyan/WCB_julia_code/WCB_selection.jl")

print("writing normalized file")
println("---------------")
write_tauWCB_var_file_normed()

include("/lustre/project/nhr-tpchange/gziyan/WCB_julia_code/WCB_functions.jl")
print("writing PECR.nc")
println("---------------")
write_PECR_file()


