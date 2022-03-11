#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("LinearAlgebra")
using HDF5,LinearAlgebra

function write_pos_data_hdf5(folder,mc_steps,particles,n,p,data,count,which,qpart=[0,[0]],log_form=false)
	println("Starting Data Write")
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	qpart_count = qpart[1]
	if log_form
		if qpart_count > 0
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$count-rad-log-$which.hdf5","w")
		else 
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$count-log.hdf5","w")
		end
	else
		if qpart_count > 0
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$count-rad-$which.hdf5","w")
		else 
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$count.hdf5","w")
		end
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
        if folder != "NA"
            cd("..")
            cd("Codes")
        end
	println("Data Added, File Closed")
end

function read_CF_hdf5(folder,mc_steps,particles,n,p,which,rad_count=1,qpart_count=0,log_form=false)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if qpart_count >= 1
		if log_form
			file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad-log-$which.hdf5","r")
		else
			file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad.hdf5","r")
		end
		qpart_data = [read(file["metadata"],"qpart_position_$i") for i in 1:qpart_count]
		full_qpart_data = [qpart_count,qpart_data]
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		full_data = [full_qpart_data,positions,wavefunc_data]
	elseif qpart_count < 1
		if log_form
			file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$which-log.hdf5","r")
		else
			file = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$which.hdf5","r")
		end
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		full_data = [positions,wavefunc_data]
	end
        if folder != "NA"
            cd("..")
            cd("Codes")
        end	
	return full_data
end

function find_CF_data(particles,n,p,rad_choice,qpart_count,log_form)
	cd("..")
	cd("cf-data-pract")
	file_list = readdir()
	number = length(file_list)
	all_data = []
	for i in 1:number
		name = file_list[i]
		separated = split(name,"-")
		
		parts_here = parse(Int,separated[6])
		if parts_here == particles
			n_here = parse(Int,separated[8])
			p_here = parse(Int,separated[10])
			if n_here == n && p_here == p
				radcount_here = parse(Int,separated[13])
				if radcount_here == rad_choice
					qpartcount_here = parse(Int,separated[12])
					if qpartcount_here == qpart_count
						which_here = parse(Int,split(separated[end],".")[1])
						mcsteps_here = parse(Int,separated[4])
						data_here = read_CF_hdf5("NA",mcsteps_here,parts_here,n_here,p_here,which_here,radcount_here,qpartcount_here,log_form)
						append!(all_data,[data_here])
					end
				end
			end
		end	
	end
	cd("..")
	cd("Codes")
	return all_data
end


function combine_CF_data(all_data,folder)
	files_count = length(all_data)
	particles = size(all_data[1][2])[1]
	qpart_data = all_data[1][1]
	full_wavefunc_data = []
	full_pos_data = Matrix{ComplexF64}(undef,particles,0)
	for i in 1:files_count
		local_wavefunc_data = all_data[i][3]
		append!(full_wavefunc_data,local_wavefunc_data)
		
		local_pos_data = all_data[i][2]
		full_pos_data = cat(full_pos_data,local_pos_data,dims=2)
		println(full_pos_data)
	end
	
	return qpart_data,full_pos_data,full_wavefunc_data
end









"fin"
