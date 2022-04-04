#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("LinearAlgebra")
using HDF5,LinearAlgebra

function write_pos_data_hdf5(folder,mc_steps,particles,n,p,data,rad_count,which,qpart=[0,[0]],log_form=false)
	println("Starting Data Write")
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	qpartcount = qpart[1]
	if log_form
		if qpartcount > 0
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-log-$which.hdf5","w")
		else 
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$rad_count-log.hdf5","w")
		end
	else
		if qpartcount > 0
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-$which.hdf5","w")
		else 
			binary_file_pos = h5open("CF-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-$rad_count.hdf5","w")
		end
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	for i in 1:qpartcount
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

function write_comb_CF_hdf5(folder,particles,n,p,data,rad_choice,log_form)
	println("Starting Data Write")
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	qpart_count = data[1][1]
	if log_form
		binary_file_pos = h5open("CF-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","w")
	else
		binary_file_pos = h5open("CF-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","w")
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	for i in 1:qpart_count
		metadata["qpart_position_$i"] = data[1][2][i]
	end
	println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["pos_x"] = real(data[2])
	alldata["pos_y"] = -imag(data[2])
	alldata["wavefunc"] = data[3]
	close(binary_file_pos)
        if folder != "NA"
            cd("..")
            cd("Codes")
        end
	println("Data Added, File Closed")
end

function read_comb_CF_hdf5(folder,particles,n,p,rad_choice,qpart_count,log_form=false)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if log_form
	file = h5open("CF-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","r")
	else
		file = h5open("CF-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","r")
	end
	qpart_data = [read(file["metadata"],"qpart_position_$i") for i in 1:qpart_count]
	full_qpart_data = [qpart_count,qpart_data]
	positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
	wavefunc_data = read(file["all-data"],"wavefunc")
	full_data = [full_qpart_data,positions,wavefunc_data]
	
        if folder != "NA"
            cd("..")
            cd("Codes")
	end	
	return full_data
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



function find_CF_data(folder,particles,n,p,rad_choice,qpart_count,log_form)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	file_list = readdir()
	number = length(file_list)
	all_data = []
	files_combined = []
	found = false
	for i in 1:number
		name = file_list[i]
		separated = split(name,"-")
		if length(separated) < 10
			continue
		end
		#parts_here = parse(Int,separated[6])
		parts_here = parse(Int,separated[findall(i->i=="part",separated)[1] + 1])
		if parts_here == particles
			#n_here = parse(Int,separated[8])
			#p_here = parse(Int,separated[10])
			n_here = parse(Int,separated[findall(i->i=="n",separated)[1] + 1])
			p_here = parse(Int,separated[findall(i->i=="p",separated)[1] + 1])
			if n_here == n && p_here == p
				#radcount_here = parse(Int,separated[13])
				if separated[3] == "comb"
					radcount_here = parse(Int,separated[findall(i->i=="rad",separated)[1] + 1])
				else
					radcount_here = parse(Int,separated[findall(i->i=="rad",separated)[1] - 1])
				end
				if radcount_here == rad_choice
					#qpartcount_here = parse(Int,separated[12])
					#println(separated[findall(i->i=="qpart",separated)[1] + 1])
					qpartcount_here = parse(Int,separated[findall(i->i=="qpart",separated)[1] + 1])
					if qpartcount_here == qpart_count
						if separated[3] == "comb"
							data_here = read_comb_CF_hdf5("NA",parts_here,n_here,p_here,radcount_here,qpartcount_here,log_form)
						else
							which_here = parse(Int,split(separated[end],".")[1])
							#mcsteps_here = parse(Int,separated[4])
							mcsteps_here = parse(Int,separated[findall(i->i=="mc",separated)[1] + 1])
							println(name)
							data_here = read_CF_hdf5("NA",mcsteps_here,parts_here,n_here,p_here,which_here,radcount_here,qpartcount_here,log_form)
						end
						append!(all_data,[data_here])
						append!(files_combined,[name])
						
					end
				end
			end
		end	
	end
	if length(files_combined) >= 2
		found = true
		#
		#mkdir("indiv-comb-cf-data")
		for i in 1:length(files_combined)
			file_name = files_combined[i]
			cd("indiv-comb-cf-data")
			repeat = findall(na->na==file_name,readdir())
			cd("..")
			if length(repeat) > 0
				println("Found overlap: $file_name")
			else
				mv("$file_name","indiv-comb-cf-data/$file_name",force=false)
			end
		end
		#
	else
		println("Only Single File Found")
	end
	if folder != "NA"
            cd("..")
            cd("Codes")
        end	
	return all_data,found
end


function combine_CF_data(all_data)
	println("Combining Data")
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
	end
	#full_wavefunc_data = reshape(full_wavefunc_data,length(full_wavefunc_data),1)
	result_wavefunc_data = fill(0.0+im*0.0,(length(full_wavefunc_data)))
	for i in 1:length(result_wavefunc_data)
		result_wavefunc_data[i] = full_wavefunc_data[i]
	end
	
	return qpart_data,full_pos_data,result_wavefunc_data
end

#=
log_form = true
particles = 16
np_vals = [[1,1],[1,2],[2,1]]
for k in 1:1
	n,p = np_vals[k]
	for j in 1:2
		qpart_count = j
		for i in 3:7
			rad_choice = i
			alldats = find_CF_data("cf-data",particles,n,p,rad_choice,qpart_count,log_form)
			if alldats[2]
				println("Found files to Combine: n=$n p=$p part=$particles rad=$rad_choice qparts=$qpart_count")
				combined = combine_CF_data(alldats[1])
				write_comb_CF_hdf5("cf-data",particles,n,p,combined,rad_choice,log_form)
				xdats = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,qpart_count,log_form)
			else
				println("No Combinations for n=$n p=$p part=$particles rad=$rad_choice qparts=$qpart_count")
			end
		end
	end
end
=#





"fin"
