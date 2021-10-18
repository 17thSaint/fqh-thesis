#import Pkg; Pkg.add("HDF5")
using HDF5

function start_rand_config(num_parts,m)
	rm = sqrt(2*num_parts*m)
	config = [[rand(Float64)*rand(-1:2:1)*rm,rand(Float64)*rand(-1:2:1)*rm] for i = 1:num_parts]
	return config
end

function dist_btw(part_1,part_2)
	return sqrt((part_1[1] - part_2[1])^2 + (part_1[2] - part_2[2])^2)
end

function prob_wavefunc(config, m, qhole)
	full = 0
	for j = 1:size(config)[1]
		for i = 1:j-1
			dist = dist_btw(config[i],config[j]) 
			full += -2*m*log( dist )
		end
		if qhole[1] > 0
			full += -2*log( dist_btw(config[j],qhole[2]))
		end
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)
	end
	return full
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [[0.0,0.0] for i = 1:num_parts]
	shift_matrix[chosen][1] += rand(-1:2:1)*rand(Float64)*step_size
	shift_matrix[chosen][2] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,chosen,num_parts,m,step_size,qhole)
	start_ham = prob_wavefunc(config,m,qhole)
	shift_matrix = move_particle(num_parts,chosen,step_size)
	new_ham = prob_wavefunc(config+shift_matrix,m,qhole)
	rand_num = rand(Float64)
	if exp(-(new_ham - start_ham)) > rand_num
		return config+shift_matrix, 1
	else
		return config, 0
	end
	
	return "Acceptance Calculation Error"
end

function main(steps,num_parts,m,step_size,qhole)
	running_config = start_rand_config(num_parts,m)
	samp_freq = 1
	acc_rate = 0.0
	therm_time = Int(0.1*steps)
	collection_time = Int(steps*samp_freq*0.9)
	time_config_x = fill(0.0,(num_parts,Int(0.9*steps)))
	time_config_y = fill(0.0,(num_parts,Int(0.9*steps)))
	index = 1
	for i_therm in 1:therm_time
		if i_therm%(therm_time*0.05) == 0
			println("Thermalizing:"," ",100*i_therm/therm_time,"%")
		end
		for j_therm in 1:num_parts
			movement = acc_rej_move(running_config,j_therm,num_parts,m,step_size,qhole)
			running_config = movement[1]
		end
	end
	println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		if i%(collection_time*0.05) == 0
			println("Running:"," ",100*i/collection_time,"%")
		end
		for j in 1:num_parts
			movement = acc_rej_move(running_config,j,num_parts,m,step_size,qhole)
			acc_rate += movement[2]/(collection_time*num_parts)
			running_config = movement[1]
		end
		
		if i%samp_freq == 0
			time_config_x[:,index] = [running_config[x][1] for x in 1:num_parts]
			time_config_y[:,index] = [running_config[x][2] for x in 1:num_parts]
			index += 1
		end
		
	end
	

	return acc_rate,time_config_x,time_config_y
end

function write_pos_data_hdf5(axis,mc_steps,particles,m,step_size,qhole,data)
	println("Starting Data Write: $axis")
	binary_file_pos = h5open("$axis-pos-mc-$mc_steps-p-$particles-m-$m-qhole.hdf5","w")
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	metadata["mc_steps"] = mc_steps
	metadata["step_size"] = step_size
	metadata["parts"] = particles
	metadata["filling_factor"] = m
	metadata["qhole_position"] = qhole[2]
	metadata["qhole_count"] = qhole[1]
	println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["deets"] = data
	close(binary_file_pos)
	println("Data Added, File Closed: $axis")
end

mcs = 200
particles = 20
step_size = 0.5
quasihole = [1,[1.0,0.0]]
for i in 1:3
	data_here = main(mcs,particles,i,step_size,quasihole)
	write_pos_data_hdf5("x",mcs,particles,i,step_size,quasihole,data_here[2])
	write_pos_data_hdf5("y",mcs,particles,i,step_size,quasihole,data_here[3])
end




#=
println("Making Data File")

saving_data = h5open("acc-rate-data-mk2.hdf5","w")
create_group(saving_data,"metadata")
metadata = saving_data["metadata"]
metadata["mc_steps"] = 2000000
metadata["step_size"] = 0.1
metadata["description"] = "Calculate acceptance rate for a range of M values and particle numbers to figure out proper sampling frequency"
create_group(saving_data,"all-data")
alldata = saving_data["all-data"]

println("File Structure Created; Starting Data Calculation")

max_particles = 20
rates = fill(0.0,(3,max_particles))
particle_numbers = [i for i in 1:max_particles]
for j in 1:3
	create_group(alldata,"m-$j")
	m = alldata["m-$j"]
	for i in 1:max_particles
		rates[j,i] = main(2000000,i,j,0.1)
		println("Done P=",i,", ","M=",j)
	end
	m["data"] = rates[j,:]
end

close(saving_data)
=#




"fin"
