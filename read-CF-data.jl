#import Pkg; Pkg.add("HDF5")
#import Pkg; Pkg.add("LinearAlgebra")
using HDF5,LinearAlgebra

function write_pos_data_hdf5(folder,vers,mc_steps,particles,n,p,data,rad_count,which,qpart=[0,[0]],log_form=false,wb=false)
	#println("Starting Data Write")
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	qpartcount = qpart[1]
	if log_form
		if qpartcount > 0
			if wb
				binary_file_pos = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-log-$which.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-log-$which.hdf5","w")
			end
		else
			if wb 
				binary_file_pos = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-log-$which.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-log-$which.hdf5","w")
			end
		end
	else
		if qpartcount > 0
			if wb
				binary_file_pos = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-$which.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$rad_count-rad-$which.hdf5","w")
			end
		else 
			if wb
				binary_file_pos = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$which.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpartcount-$which.hdf5","w")
			end
		end
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	for i in 1:qpartcount
		metadata["qpart_position_$i"] = qpart[2][i]
	end
	#println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["pos_x"] = real(data[1])
	alldata["pos_y"] = -imag(data[1])
	alldata["wavefunc"] = data[2]
	if wb
		alldata["berry"] = data[3]
	end
	close(binary_file_pos)
        if folder != "NA"
            for s in 1:length(split(folder,"/"))	
        		cd("..")
            end
            cd("Codes")
        end
	println("Data Added, File Closed")
end

function write_berry(folder,file_name,data)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	
	binary_file = h5open(file_name,"cw")
	alldata = binary_file["all-data"]
	alldata["berry"] = data
	close(binary_file)
	
	println("Berry Data Added")
	
	if folder != "NA"
            for s in 1:length(split(folder,"/"))	
        		cd("..")
            end
            cd("Codes")
        end
end

function write_comb_CF_hdf5(folder,vers,particles,n,p,data,rad_choice,log_form,wb=false)
	println("Starting Data Write")
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	qpart_count = data[1][1]
	if log_form
		if qpart_count == 0
			binary_file_pos = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-log.hdf5","w")
		else
			if wb
				binary_file_pos = h5open("$vers-pos-wb-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","w")
			end
		end
	else
		if qpart_count == 0
			binary_file_pos = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count.hdf5","w")
		else
			if wb
				binary_file_pos = h5open("$vers-pos-wb-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","w")
			else
				binary_file_pos = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","w")
			end
		end
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
	if wb
		alldata["berry"] = data[4]
	end
	close(binary_file_pos)
        if folder != "NA"
            for s in 1:length(split(folder,"/"))	
        		cd("..")
            end
            cd("Codes")
        end
	println("Data Added, File Closed")
end

function read_comb_CF_hdf5(folder,vers,particles,n,p,rad_choice,qpart_count,log_form=false,wb=false)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if log_form
		if qpart_count == 0
			file = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-log.hdf5","r")
		else
			if wb
				file = h5open("$vers-pos-wb-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","r")
			else
				file = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice-log.hdf5","r")
			end
		end
	else
		if qpart_count == 0
			file = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count.hdf5","r")
		else
			if wb
				file = h5open("$vers-pos-wb-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","r")
			else
				file = h5open("$vers-pos-comb-part-$particles-n-$n-p-$p-qpart-$qpart_count-rad-$rad_choice.hdf5","r")
			end
		end
	end
	qpart_data = [convert(Complex{Float64},read(file["metadata"],"qpart_position_$i")) for i in 1:qpart_count]
	full_qpart_data = [qpart_count,qpart_data]
	positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
	wavefunc_data = read(file["all-data"],"wavefunc")
	if wb
		berry_data = read(file["all-data"],"berry")
		full_data = [full_qpart_data,positions,wavefunc_data,berry_data]
	else
		full_data = [full_qpart_data,positions,wavefunc_data]	
	end
	close(file)
	
        if folder != "NA"
            for s in 1:length(split(folder,"/"))	
        		cd("..")
            end
            cd("Codes")
	end	
	return full_data
end

function read_CF_hdf5(folder,vers,mc_steps,particles,n,p,which,rad_count=1,qpart_count=0,log_form=false,wb=false)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if qpart_count >= 1
		if log_form
			if wb
				file = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad-log-$which.hdf5","r")
			else
				file = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad-log-$which.hdf5","r")
			end
		else
			if wb
				file = h5open("$vers-pos-wb-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad-$which.hdf5","r")
			else
				file = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$rad_count-rad-$which.hdf5","r")
			end
		end
		qpart_data = [read(file["metadata"],"qpart_position_$i") for i in 1:qpart_count]
		full_qpart_data = [qpart_count,qpart_data]
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		if wb
			berry_data = read(file["all-data"],"berry")
			full_data = [full_qpart_data,positions,wavefunc_data,berry_data]
		else
			full_data = [full_qpart_data,positions,wavefunc_data]
		end
	#
	elseif qpart_count < 1
		if log_form
			file = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-log-$which.hdf5","r")
		else
			file = h5open("$vers-pos-mc-$mc_steps-part-$particles-n-$n-p-$p-qpart-$qpart_count-$which.hdf5","r")
		end
		positions = read(file["all-data"],"pos_x") - im.*read(file["all-data"],"pos_y")
		wavefunc_data = read(file["all-data"],"wavefunc")
		qpart_data = [0,[0.0]]
		full_data = [qpart_data,positions,wavefunc_data]
	end
	#
	close(file)
        if folder != "NA"
        	for s in 1:length(split(folder,"/"))	
        		cd("..")
        	end
            	cd("Codes")
        end	
	return full_data
end

function rename_comb_file(given_file)
	this_file = h5open("$given_file","r")
	total_length = length((this_file["all-data"])["wavefunc"])
	og_split = split(given_file,".")
	new_name = og_split[1] * "-len-$total_length." * og_split[2]
	close(this_file)
	return new_name
end

function rename_pos_file(given_file)
	all_atrs = split(given_file,".")[1]
	sepd_atrs = split(all_atrs,"-")
	if length(sepd_atrs[end]) < 2
		sepd_atrs[end] = "10" * sepd_atrs[end]
	else
		sepd_atrs[end] = "1" * sepd_atrs[end]
	end
	new_name_array = ["" for i in 1:Int(2*length(sepd_atrs)-1)]
	for i in 1:Int(2*length(sepd_atrs)-1)
		if i%2 == 0
			new_name_array[i] = "-"
		else
			new_name_array[i] = sepd_atrs[Int((i+1)/2)] 
		end
	end
	new_name = prod(new_name_array) * ".hdf5"
	return new_name
end

function find_CF_data(folder,vers,particles,n,p,rad_choice,qpart_count,log_form,wb=false)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	low_vers = lowercase(vers)
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
		if separated[1] != vers
			continue
		end
		check_wb = length(findall(b->b=="wb",separated))
		if !wb || check_wb > 0
			#parts_here = parse(Int,separated[6])
			parts_here = parse(Int,separated[findall(i->i=="part",separated)[1] + 1])
			if parts_here == particles
				#n_here = parse(Int,separated[8])
				#p_here = parse(Int,separated[10])
				n_here = parse(Int,separated[findall(i->i=="n",separated)[1] + 1])
				p_here = parse(Int,separated[findall(i->i=="p",separated)[1] + 1])
				if n_here == n && p_here == p
					qpartcount_here = parse(Int,separated[findall(i->i=="qpart",separated)[1] + 1])
					if qpartcount_here == qpart_count
					#radcount_here = parse(Int,separated[13])
						#qpartcount_here = parse(Int,separated[12])
						#println(separated[findall(i->i=="qpart",separated)[1] + 1])
						if qpartcount_here != 0
							if separated[3] == "comb"
								radcount_here = parse(Int,separated[findall(i->i=="rad",separated)[1] + 1])
							else
								radcount_here = parse(Int,separated[findall(i->i=="rad",separated)[1] - 1])
							end
						else
							radcount_here = rad_choice
						end
						if radcount_here == rad_choice
							if separated[3] == "comb"
								data_here = read_comb_CF_hdf5("NA",vers,parts_here,n_here,p_here,radcount_here,qpartcount_here,log_form)
							else
								which_here = parse(Int,split(separated[end],".")[1])
								#mcsteps_here = parse(Int,separated[4])
								mcsteps_here = parse(Int,separated[findall(i->i=="mc",separated)[1] + 1])
								println(name)
								data_here = read_CF_hdf5("NA",vers,mcsteps_here,parts_here,n_here,p_here,which_here,radcount_here,qpartcount_here,log_form,wb)
							end
							append!(all_data,[data_here])
							append!(files_combined,[name])
							
						end
					end
				end
			end
		end	
	end
	if length(files_combined) >= 2
		found = true
		#
		if !isdir("indiv-comb-$low_vers-data")
			println("No Existing Indiv $vers Folder")
			mkdir("indiv-comb-$low_vers-data")
		end
		for i in 1:length(files_combined)
			file_name = files_combined[i]
			cd("indiv-comb-$low_vers-data")
			repeat = findall(na->na==file_name,readdir())
			cd("..")
			if length(repeat) > 0
				println("Found overlap: $file_name")
				if split(file_name,"-")[3] == "comb"
					println("Moving Comb File")
					new_name = rename_comb_file(file_name)
					mv("$file_name","indiv-comb-$low_vers-data/$new_name",force=false)
				else
					new_name = rename_pos_file(file_name)
					println("Moving $file_name to $new_name")
					mv("$file_name","indiv-comb-$low_vers-data/$new_name",force=false)
				end
			else
				if split(file_name,"-")[3] == "comb"
					new_name = rename_comb_file(file_name)
					mv("$file_name","indiv-comb-$low_vers-data/$new_name",force=false)
				else
					mv("$file_name","indiv-comb-$low_vers-data/$file_name",force=false)
				end
			end
		end
		#
	else
		println("Only Single File Found")
	end
	if folder != "NA"
		for s in 1:length(split(folder,"/"))	
        		cd("..")
            	end
        	cd("Codes")
        end	
	return all_data,found
end


function combine_CF_data(all_data)
	println("Combining Data")
	wb = false
	files_count = length(all_data)
	particles = size(all_data[1][2])[1]
	qpart_data = all_data[1][1]
	full_wavefunc_data = []
	full_pos_data = Matrix{ComplexF64}(undef,particles,0)
	if length(all_data[1]) > 3
		wb = true
		full_berry_data = []
	end
	for i in 1:files_count
		local_wavefunc_data = all_data[i][3]
		append!(full_wavefunc_data,local_wavefunc_data)
		
		local_pos_data = all_data[i][2]
		full_pos_data = cat(full_pos_data,local_pos_data,dims=2)
		if wb
			append!(full_berry_data,all_data[i][4])
		end
	end
	#full_wavefunc_data = reshape(full_wavefunc_data,length(full_wavefunc_data),1)
	result_wavefunc_data = fill(0.0+im*0.0,(length(full_wavefunc_data)))
	if wb
		result_berry_data = fill(0.0+im*0.0,(length(full_berry_data)))
	end
	for i in 1:length(result_wavefunc_data)
		result_wavefunc_data[i] = full_wavefunc_data[i]
		if wb
			result_berry_data[i] = full_berry_data[i]
		end
	end
	if wb
		full_data = [qpart_data,full_pos_data,result_wavefunc_data,result_berry_data]
	else
		full_data = [qpart_data,full_pos_data,result_wavefunc_data]
	end
	return full_data
end

function main(args)
	args = string(args)
	if args == "T"
		return true
	elseif args == "F"
		return false
	end
end
make_new = main(ARGS[1])

if make_new
	println("Combining Stuff")
	log_form_comb = true
	berry_comb = false
	particles_comb = 22
	vers_comb = "CF"
	low_vers_comb = lowercase(vers_comb)
	np_vals_comb = [[1,1],[1,2],[2,1]]
	for k in 1:1
		n_comb,p_comb = np_vals_comb[k]
		for j in 2:2
			qpart_count_comb = j
			for i in [1]
				rad_choice_comb = i
				alldats_comb = find_CF_data("$low_vers_comb-data",vers_comb,particles_comb,n_comb,p_comb,rad_choice_comb,qpart_count_comb,log_form_comb,berry_comb)
				if alldats_comb[2]
					println("Found files to Combine: $vers_comb n=$n_comb p=$p_comb part=$particles_comb rad=$rad_choice_comb qparts=$qpart_count_comb WB=$berry_comb")
					combined = combine_CF_data(alldats_comb[1])
					write_comb_CF_hdf5("$low_vers_comb-data",vers_comb,particles_comb,n_comb,p_comb,combined,rad_choice_comb,log_form_comb,berry_comb)
					#xdats = read_comb_CF_hdf5("$low_vers_comb-data",vers_comb,particles_comb,n_comb,p_comb,rad_choice_comb,qpart_count_comb,log_form_comb)
				else
					println("No Combinations for n=$n_comb p=$p_comb part=$particles_comb rad=$rad_choice_comb qparts=$qpart_count_comb WB=$berry_comb")
				end
			end
		end
	end

end




"fin"
