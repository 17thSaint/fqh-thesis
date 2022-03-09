using LaTeXStrings,Statistics

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
	shift_matrix = move_particle(num_parts,chosen,step_size)
	rand_num = rand(Float64)	

	#
	if log_form
		start_wavefunc = get_wavefunc_fromlog(config,n,p,qpart)
		start_ham = 2*real(start_wavefunc)
		new_wavefunc = get_wavefunc_fromlog(config+shift_matrix,n,p,qpart)
		new_ham = 2*real(new_wavefunc)
		check = new_ham - start_ham
		rand_num = log(rand_num)*2*0.5
	else
		start_wavefunc = get_wavefunc(config,n,p,qpart)
		start_ham = abs2(start_wavefunc)
		new_wavefunc = get_wavefunc(config+shift_matrix,n,p,qpart)
		new_ham = abs2(new_wavefunc)
		check = new_ham/start_ham
	end
	#
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
	therm_time = 0#Int(0.0001*steps)
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
		if i%(collection_time*0.001) == 0
			println("Running $n $p:"," ",100*i/collection_time,"%, Acc Rate: ",acc_count,"/",num_parts*i)
		end
	end
	

	return time_config,time_wavefunc
end

particles = 16
mc_steps = 1000
step_size = 0.5
log_form = true
np_vals = [[1,1],[1,2],[2,1]]
which_np = 1#parse(Int64,ARGS[1])
n,p = np_vals[which_np]
fill_denom = 2*n*p + 1
filling = n/(2*p*n+1)
rm = sqrt(2*particles/filling)
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
Threads.@threads for k in 1:10
	println("$n/$fill_denom: ",k)
	rad_choice = k
	x_rad = x_rads[rad_choice]
	println(x_rad)
	qpart = [2,[x_rad+im*0.0,0.0+im*0.0]]
	rezz = main(n,p,mc_steps,particles,step_size,qpart,log_form)
	write_pos_data_hdf5("cf-data",mc_steps,particles,n,p,rezz,k,qpart,log_form)
end





"fin"
