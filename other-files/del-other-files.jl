
all_files = readdir()
keep = ["cf-wavefunc.jl","mc-cf.jl","read-CF-data.jl","write-accmat-hdf5.jl"]
for a in all_files
	if a !in keep
		rm(a)
	end
end
