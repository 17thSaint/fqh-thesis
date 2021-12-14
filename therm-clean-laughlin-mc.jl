#import Pkg; Pkg.add("HDF5")
using HDF5,PyPlot

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
		for k in 1:qhole[1]
			full += -2*log( dist_btw(config[j],qhole[k+1]))
		end
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)
	end
	return full
end

function wavefunc_CFn2_N3(config,p)
	imag_config = [config[i][1] - im*config[i][2] for i in 1:3]
	full = log(16*(p^2)) + 2*log(abs(2*imag_config[1]^2 - imag_config[2]^2 - imag_config[3]^2 + 4*imag_config[2]*imag_config[3] - 2*imag_config[1]*(imag_config[2]+imag_config[3])))
	for j in 1:3
		for i in 1:j-1
			dist = dist_btw(config[i],config[j])
			full += (4*p-2)*log(dist)
		end
		full += -0.5 * (config[j][1]^2 + config[j][2]^2) 
	end
	return full
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [[0.0,0.0] for i = 1:num_parts]
	shift_matrix[chosen][1] += rand(-1:2:1)*rand(Float64)*step_size
	shift_matrix[chosen][2] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move_CF(config,p,chosen,num_parts,step_size)
	starting = wavefunc_CFn2_N3(config,p)
	shift_matrix = move_particle(num_parts,chosen,step_size)
	next = wavefunc_CFn2_N3(config+shift_matrix,p)
	rand_num = log(rand(Float64))
	if next - starting > rand_num
		return config+shift_matrix, 1
	else
		return config, 0
	end
	
	return "Acceptance Calculation Error"
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
	samp_freq = 10
	#acc_rate = 0.0
	therm_time = Int(0.1*steps)
	collection_time = Int(steps*0.9)
	time_config_x = fill(0.0,(num_parts,Int(0.9*steps/samp_freq)))
	time_config_y = fill(0.0,(num_parts,Int(0.9*steps/samp_freq)))
	index = 1
	for i_therm in 1:therm_time
		#if i_therm%(therm_time*0.05) == 0
		#	println("Thermalizing:"," ",100*i_therm/therm_time,"%")
		#end
		for j_therm in 1:num_parts
			#movement = acc_rej_move(running_config,j_therm,num_parts,m,step_size,qhole)
			movement = acc_rej_move_CF(running_config,m,j_therm,num_parts,step_size)
			running_config = movement[1]
		end
	end
	println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		#if i%(collection_time*0.05) == 0
		#	println("Running:"," ",100*i/collection_time,"%")
		#end
		#println("Calc New Config",DateTime(now()))
		for j in 1:num_parts
			#movement = acc_rej_move(running_config,j,num_parts,m,step_size,qhole)
			movement = acc_rej_move_CF(running_config,m,j,num_parts,step_size)
			#acc_rate += movement[2]/(collection_time*num_parts)
			running_config = movement[1]
		end
		#println("Found New Config",DateTime(now()))
		#println("Checking to add Data",DateTime(now()))
		if i%samp_freq == 0
			time_config_x[:,index] = [running_config[x][1] for x in 1:num_parts]
			time_config_y[:,index] = [running_config[x][2] for x in 1:num_parts]
			index += 1
		end
		#println("Data Added",DateTime(now()))
	end
	

	return time_config_x,time_config_y#,acc_rate
end

function write_pos_data_hdf5(axis,mc_steps,particles,m,step_size,qhole,data,count,q2loc=1)
	println("Starting Data Write: $axis")
	qhole_count = qhole[1]
	if qhole_count > 1
		#binary_file_pos = h5open("$axis-pos-mc-$mc_steps-p-$particles-m-$m-qhole-$qhole_count-q2loc-$q2loc-$count.hdf5","w")
		binary_file_pos = h5open("$axis-pos-p-$particles-m-$m-qhole-long-$qhole_count-q2loc-$q2loc-$count.hdf5","w")
	else 
		binary_file_pos = h5open("$axis-pos-mc-$mc_steps-p-$particles-m-$m-qhole-$count.hdf5","w")
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	metadata["mc_steps"] = mc_steps
	metadata["step_size"] = step_size
	metadata["parts"] = particles
	metadata["filling_factor"] = m
	for i in 1:qhole_count
		metadata["qhole_position_$i"] = qhole[i+1]
	end
	metadata["qhole_count"] = qhole[1]
	println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["deets"] = data
	close(binary_file_pos)
	println("Data Added, File Closed: $axis")
end



mcs = 100000
particles = 3
step_size = 0.5
#=
rezz = []
Threads.@threads for i in 1:3
	results = main(mcs,particles,i,step_size,0.0)
	append!(rezz,results)
end
=#


# single set 30 seconds, 3 sets in 1 minute
# need to set environment variable export JULIA_NUM_THREADS=n
# or use julia -t n

#=
q_rad_count = 9
i = parse(Int64,ARGS[1])
rm = sqrt(2*particles*i)
qhole2_locations = [[0,0],[-rm*0.5,0],[-1.4*rm,0]]
qhole2_choice = 1
for j in 1:q_rad_count
	quasihole = [2,[j*rm*1.5/q_rad_count,0],qhole2_locations[qhole2_choice]]
	data_here = main(mcs,particles,i,step_size,quasihole)
	write_pos_data_hdf5("x",mcs,particles,i,step_size,quasihole,data_here[1],j,qhole2_choice)
	write_pos_data_hdf5("y",mcs,particles,i,step_size,quasihole,data_here[2],j,qhole2_choice)
end
=#

# 1 minute 15 sec = 1.25 1.25*1000= 1250/60= 20.83


#=  Acceptance rate data collection
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
