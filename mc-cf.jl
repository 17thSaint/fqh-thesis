using HDF5,LaTeXStrings

include("cf-wavefunc.jl")

function write_pos_data_hdf5(mc_steps,particles,n,step_size,qpart,data,count)
	println("Starting Data Write")
	qpart_count = qpart[1]
	if qpart_count > 0
		binary_file_pos = h5open("CF-pos-p-$particles-n-$n-qpart-$qpart_count-$count.hdf5","w")
	else 
		binary_file_pos = h5open("CF-pos-mc-$mc_steps-p-$particles-n-$n-$count.hdf5","w")
	end
	create_group(binary_file_pos,"metadata")
	metadata = binary_file_pos["metadata"]
	for i in 1:qpart_count
		metadata["qpart_position_$i"] = qpart[2][i]
	end
	println("Metadata Added")
	create_group(binary_file_pos,"all-data")
	alldata = binary_file_pos["all-data"]
	alldata["deets_x"] = real(data)
	alldata["deets_y"] = -imag(data)
	close(binary_file_pos)
	println("Data Added, File Closed")
end


function start_rand_config(num_parts,n)
	filling = n/(2*1*n+1)
	rm = sqrt(2*num_parts/filling)
	config = [rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm for i in 1:num_parts]
	return config
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,n,chosen,num_parts,step_size,qpart=[0,[0]],qhole=[0,0])
	start_ham = abs2(get_wavefunc(config,n,qpart,qhole))
	shift_matrix = move_particle(num_parts,chosen,step_size)
	new_ham = abs2(get_wavefunc(config+shift_matrix,n,qpart,qhole))
	rand_num = rand(Float64)
	if new_ham/start_ham >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return config+shift_matrix, 1
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return config, 0
	end
	
	return "Acceptance Calculation Error"
end

function main(n,steps,num_parts,step_size,qpart=[0,[0]],qhole=[0,0])
	running_config = start_rand_config(num_parts,n)
	samp_freq = Int(0.0001*steps)
	acc_count = 0
	therm_time = Int(0.01*steps)
	collection_time = steps-therm_time
	time_config = fill(0.0+im*0.0,(num_parts,Int(collection_time/samp_freq)))
	#energies = [0.0 for i in 1:Int(collection_time/samp_freq)]
	index = 1
	for i_therm in 1:therm_time
		if i_therm%(therm_time*0.05) == 0
			println("Thermalizing:"," ",100*i_therm/therm_time,"%")
		end
		for j_therm in 1:num_parts
			movement = acc_rej_move(running_config,n,j_therm,num_parts,step_size,qpart,qhole)
			running_config = movement[1]
		end
	end
	println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		#println("Calc New Config",DateTime(now()))
		for j in 1:num_parts
			movement = acc_rej_move(running_config,n,j,num_parts,step_size,qpart,qhole)
			acc_count += movement[2]
			running_config = movement[1]
		end
		acc_rate = acc_count/(num_parts*i)
		#println("Found New Config",DateTime(now()))
		#println("Checking to add Data",DateTime(now()))
		if i%samp_freq == 0
			#local_prob = abs2(get_wavefunc(running_config))
			#energies[index] = local_prob
			time_config[:,index] = [running_config[x] for x in 1:num_parts]
			index += 1
		end
		if i%(collection_time*0.05) == 0
			println("Running:"," ",100*i/collection_time,"%, Acc Rate: ",acc_count,"/",num_parts*i)
		end
		#println("Data Added",DateTime(now()))
	end
	

	return time_config#,energies#,acc_rate
end

particles = 4
mc_steps = 1000000
step_size = 0.5
n = 2
filling = round(n/(2*1*n+1),digits=3)
rm = 0.1*sqrt(2*particles/filling)
qpart = [1,[0.0+im*1.0]]#rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm
rezz = main(n,mc_steps,particles,step_size,qpart)
#write_pos_data_hdf5(mc_steps,particles,n,step_size,qpart,rezz,1)
#=
plot(real(transpose(rezz)),imag(transpose(rezz)))
for i in 1:length(qpart[2])
	scatter([real(qpart[2][i])],[imag(qpart[2][i])])
end
=#
xs = collect(Iterators.flatten([real(rezz[i,:]) for i in 1:particles]))
ys = collect(Iterators.flatten([-imag(rezz[i,:]) for i in 1:particles]))
hist2D(xs,ys,bins=100)
title(latexstring("Histogram with Quasihole and \$ \\nu = $filling \$: Mine"))

#=
starting_config = start_rand_config(particles,n)
data_count = 20
xs = [-1.2*sqrt(2*particles*5/2) + i*(2*1.2*sqrt(2*particles*5/2))/data_count for i in 0:data_count]
xs_plot = []
ys_plot = []
probs = []
for i in 1:length(xs)
	local_x = xs[i]
	for j in 1:length(xs)
		#println(i,", ",j)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		prob_local = abs2(get_wavefunc(starting_config,n,qpart))
		append!(probs,[prob_local])
	end
end

scatter3D(xs_plot,ys_plot,probs)
=#




"fin"
