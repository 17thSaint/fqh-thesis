using HDF5, DelimitedFiles

function read_data(num_parts,m,axis)
	cd("MC-Data/MC-Data/")
	if axis == 1
		stuff = readdlm("xpos-P$num_parts-M$m.csv",',',Float64)
	else
		stuff = readdlm("ypos-P$num_parts-M$m.csv",',',Float64)
	end
	cd("../../")
	return stuff
	#println("Done Reading M=",m,", ","P=",num_parts)
end

function read_hist_data(num_parts,m)
	cd("MC-Data/MC-Data/")
	stuff = readdlm("histo-P$num_parts-M$m.csv",',',Float64)
	cd("../../")
	return stuff
	#println("Done Reading M=",m,", ","P=",num_parts)
end

function read_hdf5_data(num_parts,m,data_type)
	file = h5open("$type-mc2000000-notherm.hdf5","r")
	return read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
end


binary_file_ypos = h5open("ypos-mc2000000-notherm.hdf5","w")
create_group(binary_file_ypos,"metadata")
metadata = binary_file_ypos["metadata"]
metadata["mc_steps"] = 2000000
metadata["step_size"] = 0.1
println("Metadata Added")
 
create_group(binary_file_ypos,"all-data")
alldata = binary_file_ypos["all-data"]

for i in 1:3
	create_group(alldata,"m-$i")
	m = alldata["m-$i"]
	for j in 2:8
		create_group(m,"parts-$j")
		parts = m["parts-$j"] 
		parts["data"] = read_data(j,i,0)
		println("Added P=",j,", ","M=",i)
		parts = 0
	end
	m = 0
end

close(binary_file_ypos)




"fin"
