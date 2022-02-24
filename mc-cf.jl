using HDF5,LaTeXStrings

include("cf-wavefunc.jl")

function start_rand_config(num_parts)
	rm = sqrt(2*num_parts*5/2)
	config = [rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm for i in 1:num_parts]
	return config
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,chosen,num_parts,step_size,qpart=[0,[0]],qhole=[0,0])
	start_ham = abs2(get_wavefunc(config,qpart,qhole))
	shift_matrix = move_particle(num_parts,chosen,step_size)
	new_ham = abs2(get_wavefunc(config+shift_matrix,qpart,qhole))
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

function main(steps,num_parts,step_size,qpart=[0,[0]],qhole=[0,0])
	running_config = start_rand_config(num_parts)
	samp_freq = 10
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
			movement = acc_rej_move(running_config,j_therm,num_parts,step_size,qpart,qhole)
			running_config = movement[1]
		end
	end
	println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		#println("Calc New Config",DateTime(now()))
		for j in 1:num_parts
			movement = acc_rej_move(running_config,j,num_parts,step_size,qpart,qhole)
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
mc_steps = 500000
step_size = 0.5
rm = 0.1*sqrt(2*particles*5/2)
qpart = [1,[1.0+im*1.0]]#rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm
rezz = main(mc_steps,particles,step_size,qpart)
#=
plot(real(transpose(rezz)),imag(transpose(rezz)))
for i in 1:length(qpart[2])
	scatter([real(qpart[2][i])],[imag(qpart[2][i])])
end
=#
xs = collect(Iterators.flatten([real(rezz[i,:]) for i in 1:particles]))
ys = collect(Iterators.flatten([-imag(rezz[i,:]) for i in 1:particles]))
hist2D(xs,ys,bins=50)
#title(latexstring("Histogram of Particle Position for \$ \\nu = 2/5 \$ and Quasihole"))

#=
starting_config = start_rand_config(particles)
data_count = 10
xs = [-1.2*sqrt(2*particles*5/2) + i*(2*1.2*sqrt(2*particles*5/2))/data_count for i in 0:data_count]
xs_plot = []
ys_plot = []
probs = []
for i in 1:length(xs)
	local_x = xs[i]
	for j in 1:length(xs)
		println(i,", ",j)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		prob_local = abs2(get_wavefunc(starting_config))
		append!(probs,[prob_local])
	end
end

scatter3D(xs_plot,ys_plot,probs)
=#




"fin"
