#import Pkg; Pkg.add("Statistics")
using Statistics,PyPlot,LaTeXStrings

include("cf-wavefunc.jl")
include("write-accmat-hdf5.jl")
ARGS = false
include("read-CF-data.jl")

function move_particle(num_parts::Int,chosen::Int,step_size::Float64)
	shift_matrix::Vector{ComplexF64} = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(vers::String,config::Vector{ComplexF64},n::Int,p::Int,chosen::Int,step_size::Float64,start_wavefunc::ComplexF64,reject_sets_matrix=Matrix{Vector{Int}}(undef),all_pascal=[],all_deriv_orders=[],qpart=[0,[0]],log_form=false)
	num_parts = length(config)
	shift_matrix = move_particle(num_parts,chosen,step_size)
	rand_num = rand(Float64)	

	#
	if log_form
		if vers == "CF"
			new_wavefunc = get_wavefunc_fromlog(config+shift_matrix,n,p,qpart)
		elseif vers == "RFA"
			new_wavefunc = get_rf_wavefunc(config+shift_matrix,reject_sets_matrix,all_pascal,all_deriv_orders,qpart,log_form)
		end
		start_ham = 2*real(start_wavefunc)
		new_ham = 2*real(new_wavefunc)
		check = new_ham - start_ham
		if isinf(new_ham) && isinf(start_ham)
			println("Both Inf")
			return true, config, 0, start_wavefunc
		elseif isinf(new_ham)
			if new_ham > 0
				println("New Pos Inf -> Accept")
				return true, config+shift_matrix, 1, new_wavefunc
			else
				#println("New Neg Inf -> Reject")
				return true, config, 0, start_wavefunc
			end
		elseif isinf(start_ham)
			if start_ham > 0
				#println("Start Pos Inf -> Reject")
				return true, config, 0, start_wavefunc
			else
				println("Start Neg Inf -> Accept")
				return true, config+shift_matrix, 1, new_wavefunc
			end
		end
		rand_num = log(rand_num)*2*0.5
	else
		if vers == "RFA"
			new_wavefunc = get_rf_wavefunc(config+shift_matrix,reject_sets_matrix,all_pascal,all_deriv_orders,qpart)
		elseif vers == "CF"
			new_wavefunc = get_wavefunc(config+shift_matrix,n,p,qpart)
		end
		start_ham = abs2(start_wavefunc)
		new_ham = abs2(new_wavefunc)
		check = new_ham/start_ham
	end
	#
	if check >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return true, config+shift_matrix, 1, new_wavefunc
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return true, config, 0, start_wavefunc
	end
	
	return "Acceptance Calculation Error"
end

function main(vers,n,p,steps,num_parts,step_size,rad_count,qpart=[0,[0]],log_form=false)
	allowed_sets_matrix = Matrix{Vector{Int}}(undef,(0,0))
	full_pasc_tri = Vector{Vector{Int}}(undef,0)
	full_deriv_ords = Vector{Vector{Any}}(undef,0)
	if vers == "RFA"
		allowed_sets_matrix = get_full_acc_matrix(num_parts)
		full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:num_parts]
		full_derivs = get_deriv_orders_matrix(num_parts)
		println("Made All Presets")
	end
	starting_check = true
	low_vers = lowercase(vers)
	running_config = start_rand_config(num_parts,n,p)
	start_count = 0
	output_file_count = 100
	if steps/output_file_count < 1
		steps_per_file = 1
	else
		steps_per_file = Int(steps/output_file_count)
	end
	while starting_check
		start_count += 1
		if vers == "CF"
		starting_wavefunc = get_wavefunc_fromlog(running_config,n,p,qpart)
		elseif vers == "RFA"
		starting_wavefunc = get_rf_wavefunc(running_config,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
		end
		if isinf(real(starting_wavefunc))
			running_config = start_rand_config(num_parts,n,p)
			#println(running_config[1])
		else
			#println("Started in $start_count steps")
			starting_check = false
		end
	end
	next_wavefunc = get_rf_wavefunc(running_config,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
	filling = n/(2*p*n+1)
	rm = sqrt(2*num_parts/filling)
	lstar = sqrt(2*p*n+1)
	wavefunc = 0.0+im*0.0
	samp_freq = 1#Int(steps/samp_count)
	acc_count = 0
	therm_time = 100#Int(0.0001*steps)
	collection_time = steps
	time_config = fill(0.0+im*0.0,(num_parts,Int(collection_time/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(collection_time/samp_freq)))
	index = 1
	number = 0
	for i_therm in 1:therm_time
		#if i_therm%(therm_time*0.05) == 0
		#	println("Thermalizing:"," ",100*i_therm/therm_time,"%")
		#end
		for j_therm in 1:num_parts
			movement = acc_rej_move(vers,running_config,n,p,j_therm,step_size,next_wavefunc,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
			next_wavefunc = movement[4]
			if movement[1]
				running_config = movement[2]
			else
				println("NaN Wavefunc")
				return movement
			end
		end
	end
	#println("Thermalization Done, Starting Data Collection")
	for i in 1:collection_time
		for j in 1:num_parts
			movement = acc_rej_move(vers,running_config,n,p,j,step_size,next_wavefunc,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
			next_wavefunc = movement[4]
			if movement[1]
				acc_count += movement[3]
				running_config = movement[2]
				wavefunc = movement[4]
			else
				println("NaN Wavefunc")
				return movement
			end
		end
		acc_rate = acc_count/(num_parts*i)
		#
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:num_parts]
			time_wavefunc[index] = wavefunc
			index += 1
		end
		#=
		if i%(collection_time*0.1) == 0
			println("Running $n $p:"," ",100*i/collection_time,"%, Acc Rate: ",acc_count,"/",num_parts*i)
		end
		#
		if i%(samp_freq*steps_per_file) == 0
			number += 1
			data = [time_config[:,index-steps_per_file:index-1],time_wavefunc[index-steps_per_file:index-1]]
			write_pos_data_hdf5("NA",vers,steps,num_parts,n,p,data,rad_count,number,qpart,log_form)
		end
		=#
	end
	acc_rate = acc_count/(num_parts*steps)

	return acc_rate,time_config,time_wavefunc
end


function auto_correlation(energies, delta_t)
    average_energy = mean(energies)
    
    points = Int(floor(length(energies)-delta_t))
    
    energy_fluctuations = energies.-average_energy
    
    autocorrelation_top = [0.0 for i in 1:points]
    autocorrelation_bottom = mean(energy_fluctuations.^2)
    
    for i in 1:points
        autocorrelation_top[i] = energy_fluctuations[i]*energy_fluctuations[i+delta_t]
    end
    
    return (mean(autocorrelation_top)/autocorrelation_bottom)
end

function get_autocorr_length(wavefunc_data,samp_freq)
	energy = 2 .*real.(wavefunc_data)
	full_length = length(wavefunc_data)
	#println("Full Length = $full_length")
	len = 10*full_length
	dts = [i for i in 1:Int(0.1*len)-1]
	autocorr = [0.0 for i in 1:Int(0.1*len)-1]
	for i in 1:Int(0.1*len)-1
		autocorr[i] = auto_correlation(energy,dts[i])
	end
	check_below_tol = autocorr .< [0.01 for i in 1:length(autocorr)]
	#println(autocorr)
	corr_length = samp_freq*dts[findall(check_below_tol)[1]]
	return corr_length,dts,autocorr
end

#
particles = 4
mc_steps = 300
log_form = true
np_vals = [[1,1],[1,2],[2,1]]
which_np = 1#parse(Int64,ARGS[1])
n,p = np_vals[which_np]
fill_denom = 2*n*p + 1
filling = n/(2*p*n+1)
rm = sqrt(2*particles/filling)
step_size = rm/3.0#0.4 + 0.175*rm
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
k = 5#parse(Int64,ARGS[1])
rad_choice = k
x_rad = x_rads[rad_choice]
qpart_choices = [[0,[0.0]],[1,[x_rad+im*0.0]],[2,[x_rad+im*0.0,0.0+im*0.0]]]
qpart_selected = 2#parse(Int64,ARGS[1])
qpart = qpart_choices[qpart_selected]
#
datacount = 10
step_sizes = [(0.25 + i*0.5/datacount)*rm for i in 0:datacount]
corr_lengths = [0.0 for i in 1:datacount+1]
accrates = [0.0 for i in 1:datacount+1]
corr_length_errs = [0.0 for i in 1:datacount+1]
accrate_errs = [0.0 for i in 1:datacount+1]
for i in 1:datacount+1
	println(i/(datacount+1))
	step_size = step_sizes[i]
	howmany = 5
	local_corr_lengths = [0.0 for j in 1:howmany]
	local_accrates = [0.0 for j in 1:howmany]
	for j in 1:howmany
		rezz = main("RFA",n,p,mc_steps,particles,step_size,k,qpart,log_form)
		prob = 2*real(rezz[3])
		local_corr_lengths[j] = get_autocorr_length(prob,1)[1]/mc_steps
		local_accrates[j] = rezz[1]
		#rss = round(step_size,digits=3)
		#plot(prob./maximum(prob),label="$rss")
	end
	#
	corr_lengths[i] = mean(local_corr_lengths)
	corr_length_errs[i] = std(local_corr_lengths)
	accrates[i] = mean(local_accrates)
	accrate_errs[i] = std(local_accrates)
	#
end

figure()
errorbar(step_sizes./rm,corr_lengths,yerr=[corr_length_errs,corr_length_errs],fmt="-o",label="$particles")
title("Correlation Length")
legend()
figure()
errorbar(step_sizes,accrates,yerr=[accrate_errs,accrate_errs],fmt="-o")
title("Acceptance Rate")
#=
flat_data = Iterators.flatten([rezz[2][i,:] for i in 1:particles])
figure()
hist2D(real.(flat_data),-imag.(flat_data),bins=100)
title("$qpart")
=#

#=
data_count = 50
xs = [-1.5*rm + i*(2*1.5*rm)/data_count for i in 1:data_count]
coords_config = start_rand_config(particles,n,p)#10*(rand(particles) + im*rand(particles))
#xs = [real(coords_config[2]) - 0.3*rm + 2*0.3*rm*i/data_count for i in 0:data_count]
#ys = [-imag(coords_config[2]) - 0.3*rm + 2*0.3*rm*i/data_count for i in 0:data_count]
rf_wavefunc = fill(0.0,(data_count+1,data_count+1))
cf_wavefunc = fill(0.0,(data_count+1,data_count+1))
for i in 1:length(xs)
	println(i/length(xs))
	local_x = xs[i]
	for j in 1:length(xs)
		local_y = xs[j]
		coords_config[1] = local_x - im*local_y
			
		cf_wavefunc[i,j] = 2*real(get_wavefunc_fromlog(coords_config,n,p,qpart))
		rf_wavefunc[i,j]  = 2*real(get_rf_wavefunc(coords_config,qpart,log_form))
	end
end
figure()
imshow(cf_wavefunc)#./maximum(cf_wavefunc))
title("CF")
colorbar()
figure()
imshow(rf_wavefunc)#./maximum(rf_wavefunc))
title("RF")
colorbar()
=#






"fin"
