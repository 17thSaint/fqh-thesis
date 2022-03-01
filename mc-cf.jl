using LaTeXStrings

include("cf-wavefunc.jl")
include("read-CF-data.jl")

function start_rand_config(num_parts,n,p)
	filling = n/(2*p*n+1)
	rm = sqrt(2*num_parts/filling)
	config = [rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm for i in 1:num_parts]
	return config
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,n,p,chosen,step_size,qpart=[0,[0]],log_form=false)
	num_parts = length(config)
	start_wavefunc = get_wavefunc(config,n,p,qpart)
	start_ham = abs2(start_wavefunc)
	shift_matrix = move_particle(num_parts,chosen,step_size)
	new_wavefunc = get_wavefunc(config+shift_matrix,n,p,qpart)
	new_ham = abs2(new_wavefunc)
	check = new_ham/start_ham
	rand_num = rand(Float64)
	#=
	if log_form
		start_wavefunc = get_wavefunc_fromlog(config,n,p,qpart)
		start_ham = 2*real(start_wavefunc)
		new_wavefunc = get_wavefunc_fromlog(config+shift_matrix,n,p,qpart)
		new_ham = 2*real(new_wavefunc)
		check = new_ham - start_ham
		rand_num = log(rand_num)*2*0.5
	end
	=#
	if check >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return config+shift_matrix, 1, new_wavefunc
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return config, 0, start_wavefunc
	end
	
	return "Acceptance Calculation Error"
end

function main(n,p,steps,num_parts,step_size,qpart=[0,[0]],log_form=false)
	running_config = start_rand_config(num_parts,n,p)
	wavefunc = 0.0+im*0.0
	samp_freq = 100#Int(steps/samp_count)
	acc_count = 0
	therm_time = Int(0.001*steps)
	collection_time = steps-therm_time
	time_config = fill(0.0+im*0.0,(num_parts,Int(collection_time/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(collection_time/samp_freq)))
	index = 1
	for i_therm in 1:therm_time
		if i_therm%(therm_time*0.05) == 0
			println("Thermalizing:"," ",100*i_therm/therm_time,"%")
		end
		for j_therm in 1:num_parts
			movement = acc_rej_move(running_config,n,p,j_therm,step_size,qpart,log_form)
			running_config = movement[1]
		end
	end
	println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		for j in 1:num_parts
			movement = acc_rej_move(running_config,n,p,j,step_size,qpart,log_form)
			acc_count += movement[2]
			running_config = movement[1]
			wavefunc = movement[3]
		end
		acc_rate = acc_count/(num_parts*i)
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:num_parts]
			time_wavefunc[index] = wavefunc
			index += 1
		end
		if i%(collection_time*0.05) == 0
			println("Running:"," ",100*i/collection_time,"%, Acc Rate: ",acc_count,"/",num_parts*i)
		end
	end
	

	return time_config,time_wavefunc
end

particles = 8
mc_steps = 100000
step_size = 0.5
n = 1
p = 2
fill_denom = 2*n*p + 1
filling = round(n/(2*1*n+1),digits=3)
rm = sqrt(2*particles/filling)
for i in 0:9
	println(i)
	x_rad = 0.01*rm + i*(1.29*rm)/10
	qpart = [1,[x_rad+im*0.0]]
	rezz = main(n,p,mc_steps,particles,step_size,qpart)
	write_pos_data_hdf5(mc_steps,particles,n,p,rezz,i+1,qpart)
end


#write_pos_data_hdf5(mc_steps,particles,n,step_size,qpart,rezz,1)
#=
plot(real(transpose(rezz)),imag(transpose(rezz)))
for i in 1:length(qpart[2])
	scatter([real(qpart[2][i])],[imag(qpart[2][i])])
end
xs = collect(Iterators.flatten([real(rezz[i,:]) for i in 1:particles]))
ys = collect(Iterators.flatten([-imag(rezz[i,:]) for i in 1:particles]))
hist2D(xs,ys,bins=100)
title(latexstring("Histogram with Quasihole and \$ \\nu = $filling \$"))

radii = collect(Iterators.flatten([abs2.(rezz[i,:]) for i in 1:particles]))
hist(radii,bins=100)
title_string = ["without Quasiparticle","with Quasihole at Origin"][qpart[1]+1]
title(latexstring("Radial Distribution $title_string, \$ \\nu = $n / $fill_denom \$"))
=#
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
		prob_local = abs2(get_wavefunc(starting_config,n,qpart)[1])
		append!(probs,[prob_local])
	end
end

scatter3D(xs_plot,ys_plot,probs)
=#




"fin"
