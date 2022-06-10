#import Pkg; Pkg.add("Statistics")
using Statistics

include("cf-wavefunc.jl")
include("write-accmat-hdf5.jl")
#include("berry-cf.jl")
ARGS = "F"
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

function main(vers,n,p,steps,num_parts,step_size,rad_count,qpart=[0,[0]],log_form=false,farm=0,allowed_sets_matrix = Matrix{Any}(undef,(0,0)),full_pasc_tri = Vector{Vector{Int}}(undef,0),full_derivs = Vector{Vector{Any}}(undef,0))
	#allowed_sets_matrix = Matrix{Vector{Int}}(undef,(0,0))
	#full_pasc_tri = Vector{Vector{Int}}(undef,0)
	#full_derivs = Vector{Vector{Any}}(undef,0)
	#if vers == "RFA"
	#	allowed_sets_matrix = get_full_acc_matrix(num_parts)
	#	full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:num_parts]
	#	full_derivs = get_deriv_orders_matrix(num_parts)
		#println("Made All Presets")
	#end
	println("Starting")
	starting_check = true
	low_vers = lowercase(vers)
	running_config = start_rand_config(num_parts,n,p)
	start_count = 0
	output_file_count = 10
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
	if vers == "RFA"
		next_wavefunc = get_rf_wavefunc(running_config,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
	else
		next_wavefunc = get_wavefunc_fromlog(running_config,n,p,qpart)
	end
	if vers == "CF"
		filling = n/(2*p*n+1)
		denom = 2*p*n+1
		rm = sqrt(2*num_parts/filling)
		lstar = sqrt(2*p*n+1)
	elseif vers == "RFA"
		filling = n/(2*p*n-1)
		denom = 2*p*n-1
		rm = sqrt(2*num_parts/filling)
		lstar = sqrt(2*p*n-1)
	end
	wavefunc = 0.0+im*0.0
	samp_freq = 1#Int(steps/samp_count)
	acc_count = 0
	therm_time = 0#50#Int(0.0001*steps)
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
	println("Thermalization Done, Starting Data Collection")
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
			#println("Moved Particle $j in Step $i")
		end
		acc_rate = acc_count/(num_parts*i)
		#
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:num_parts]
			time_wavefunc[index] = wavefunc
			index += 1
			#println("Added Data for Sampling Frequency")
		end
		#
		if i%(collection_time*0.1) == 0
			println("Running $n/$denom:"," ",100*i/collection_time,"%, Acc Rate: ",acc_count,"/",num_parts*i)
		end
		#
		if i%(samp_freq*steps_per_file) == 0
			number += 1
			data = [time_config[:,index-steps_per_file:index-1],time_wavefunc[index-steps_per_file:index-1]]
			folderhere = lowercase(vers)
			focused_rad_count = 30 + rad_count
			write_pos_data_hdf5("NA",vers,steps,num_parts,n,p,data,focused_rad_count,number+farm*output_file_count,qpart,log_form)
			#write_pos_data_hdf5("NA",vers,steps,num_parts,n,p,data,rad_count,number,qpart,log_form)
		end
		#
	end

	acc_rate = acc_count/(num_parts*steps)

	return acc_rate,time_config,time_wavefunc
end

particles = 4
flux_type = "RFA"
which_np = 2
log_form = true
#
mc_steps = 1000

np_vals = [[1,1],[1,2],[2,1]]
n,p = np_vals[which_np]
if flux_type == "RFA"
	fill_denom = 2*n*p - 1
	filling = n/(2*p*n-1)
	allowed_sets_matrix = get_full_acc_matrix(particles)
	full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
	full_derivs = get_deriv_orders_matrix(particles)
elseif flux_type == "CF"
	fill_denom = 2*n*p + 1
	filling = n/(2*p*n+1)
end
rm = sqrt(2*particles/filling)
step_size = rm/3.0#0.4 + 0.175*rm
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
focused_x_rads = [x_rads[3] + j*(x_rads[4]-x_rads[3])/10 for j in 1:10]

for k in 1:10
rad_choice = k
x_rad = focused_x_rads[rad_choice]
qpart_choices = [[0,[0.0]],[1,[x_rad+0.0*im]],[2,[x_rad+0.0*rm,0.0+im*0.0]]]
qpart_selected = 2#parse(Int64,ARGS[2])
qpart = qpart_choices[qpart_selected]
rezz = main(flux_type,n,p,mc_steps,particles,step_size,rad_choice,qpart,log_form,0,allowed_sets_matrix,full_pasc_tri,full_derivs)
sleep(2.0)
end



#compl = collect(Iterators.flatten([rezz[2][i,:] for i in 1:particles]))
#figure()
#hist2D(real.(compl),-imag.(compl),bins=100)



"fin"
