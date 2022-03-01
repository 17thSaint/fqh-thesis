using HDF5,PyPlot

function write_pos_data_hdf5(mc_steps,particles,n,p,data,count,qpart=[0,[0]])
	println("Starting Data Write")
	qpart_count = qpart[1]
	if qpart_count > 0
		binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$count-rad.hdf5","w")
	else 
		binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$count.hdf5","w")
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	for i in 1:qpart_count
		metadata["qpart_position_$i"] = qpart[2][i]
	end
	println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["pos_x"] = real(data[1])
	alldata["pos_y"] = -imag(data[1])
	alldata["wavefunc"] = data[2]
	close(binary_file_pos)
	println("Data Added, File Closed")
end

function read_CF_hdf5(folder,mc_steps,particles,n,p,count,qpart_count=0)
	cd("..")
	cd("$folder")
	if qpart_count >= 1
		file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$count-rad.hdf5","r")
		qpart_data = [read(file["metadata"],"qpart_position_$i") for i in 1:qpart_count]
		full_qpart_data = [qpart_count,qpart_data]
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		full_data = [full_qpart_data,positions,wavefunc_data]
	elseif qpart_count < 1
		file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$count.hdf5","r")
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		full_data = [positions,wavefunc_data]
	end
	
	cd("..")
	cd("Codes")
	return full_data
end


