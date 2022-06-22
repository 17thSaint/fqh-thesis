#import Pkg; Pkg.add("Statistics")
using Statistics

include("cf-wavefunc.jl")
include("write-accmat-hdf5.jl")
include("berry-cf.jl")
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

function main(vers,n,p,steps,num_parts,step_size,rad_count,qpart=[0,[0]],log_form=false,farm=0,wb=false,allowed_sets_matrix = Matrix{Any}(undef,(0,0)),full_pasc_tri = Vector{Vector{Int}}(undef,0),full_derivs = Vector{Vector{Any}}(undef,0))
	#allowed_sets_matrix = Matrix{Vector{Int}}(undef,(0,0))
	#full_pasc_tri = Vector{Vector{Int}}(undef,0)
	#full_derivs = Vector{Vector{Any}}(undef,0)
	#if vers == "RFA"
	#	allowed_sets_matrix = get_full_acc_matrix(num_parts)
	#	full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:num_parts]
	#	full_derivs = get_deriv_orders_matrix(num_parts)
		#println("Made All Presets")
	#end
	#println("Starting: Getting Config")
	starting_check = true
	starting_wavefunc = 0.0+im*0.0
	low_vers = lowercase(vers)
	running_config = start_rand_config(num_parts,n,p)
	start_count = 0
	output_file_count = 10
	if steps/output_file_count < 1
		steps_per_file = 1
	else
		steps_per_file = Int(steps/output_file_count)
	end
	#println("Getting First Wavefunc Value")
	while starting_check
		start_count += 1
		if vers == "CF"
		starting_wavefunc = get_wavefunc_fromlog(running_config,n,p,qpart)
		elseif vers == "RFA"
		starting_wavefunc = get_rf_wavefunc(running_config,allowed_sets_matrix,full_pasc_tri,full_derivs,qpart,log_form)
		end
		if isinf(real(starting_wavefunc))
			running_config = start_rand_config(num_parts,n,p)
		else
			println("Started in $start_count steps")
			starting_check = false
		end
	end
	#println("Found Starting Wavefunc")
	next_wavefunc = starting_wavefunc + 0.1 - 0.1
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
	if qpart[1] == 1
		qpart1_og_location = qpart[2][1] + 0.1 - 0.1
		qpart_og = [1,[qpart1_og_location]]
	elseif qpart[1] == 2
		qpart1_og_location = qpart[2][1] + 0.1 - 0.1
		qpart2_og_location = qpart[2][2] + 0.1 - 0.1
		qpart_og = [2,[qpart1_og_location,qpart2_og_location]]
	end
	wavefunc = 0.0+im*0.0
	samp_freq = 1#Int(steps/samp_count)
	acc_count = 0
	therm_time = 0#50#Int(0.0001*steps)
	collection_time = steps
	time_config = fill(0.0+im*0.0,(num_parts,Int(collection_time/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(collection_time/samp_freq)))
	if wb
		time_berry = fill(0.0,(Int(collection_time/samp_freq)))
	end
	index = 1
	number = 0
	#println("Starting Thermalization")
	#time_start = now()
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
			#println("MC Sample $i, Particle $j")
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
			time_config[:,index] = [Complex(running_config[x]) for x in 1:num_parts]
			local_config_time = time_config[:,index]
			time_wavefunc[index] = wavefunc
			if wb
				time_berry[index] = get_expval(vers,n,p,[qpart_og,local_config_time,[wavefunc]],1,-0.001,1,log_form)[1]
			end
			index += 1
			#println("Added Data for Sampling Frequency")
		end
		#
		if i%(collection_time*0.01) == 0
			println("Running $n/$denom:"," ",100*i/collection_time,"%")
		end
		#
		if i%(samp_freq*steps_per_file) == 0
			number += 1
			if wb
				data = [time_config[:,index-steps_per_file:index-1],time_wavefunc[index-steps_per_file:index-1],time_berry[index-steps_per_file:index-1]]
			else
				data = [time_config[:,index-steps_per_file:index-1],time_wavefunc[index-steps_per_file:index-1]]
			end
			folderhere = lowercase(vers)
			focused_rad_count = 0 + rad_count
			write_pos_data_hdf5("$folderhere-data",vers,steps,num_parts,n,p,data,focused_rad_count,number+farm*output_file_count,qpart,log_form,wb)
			#write_pos_data_hdf5("NA",vers,steps,num_parts,n,p,data,rad_count,number,qpart,log_form)
		end
		#
	end
	#time_end = now()
	#total_time = (time_end - time_start).value
	acc_rate = acc_count/(num_parts*steps)
	
	if wb
		return acc_rate,time_config,time_wavefunc,time_berry,total_time
	else
		return acc_rate,time_config,time_wavefunc,total_time
	end
end

particles = 30
#cluster_types = [ ["CF",1],["CF",2],["CF",3],["RFA",1],["RFA",2],["RFA",3] ]
flux_type = "CF"#cluster_types[parse(Int64,ARGS[2])][1]
which_np = 1
log_form = true
with_berry = true
#
mc_steps = 10000

np_vals = [[1,1],[1,2],[2,1]]
n,p = np_vals[which_np]
if flux_type == "RFA"
	fill_denom = 2*n*p - 1
	filling = n/(2*p*n-1)
	allowed_sets_matrix = get_full_acc_matrix(particles)
	full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
	full_derivs = get_deriv_orders_matrix(particles)
elseif flux_type == "CF"
	allowed_sets_matrix = Matrix{Any}(undef,(0,0))
	full_pasc_tri = Vector{Vector{Int}}(undef,0)
	full_derivs = Vector{Vector{Any}}(undef,0)
	fill_denom = 2*n*p + 1
	filling = n/(2*p*n+1)
end
rm = sqrt(2*particles/filling)
step_size = rm/3.0#0.4 + 0.175*rm
big_x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
starting_rad = 3
#x_rads = [big_x_rads[starting_rad] + j*(big_x_rads[starting_rad+1]-big_x_rads[starting_rad])/10 for j in 1:10]
rad_choice = parse(Int64,ARGS[2])
x_rad = big_x_rads[rad_choice]
x_rad_2 = 0.0#-x_rads[rad_choice]
qpart_choices = [[0,[0.0]],[1,[x_rad+0.0*im]],[2,[x_rad+0.0*rm,x_rad_2+im*0.0]]]
qpart_selected = parse(Int64,ARGS[3])
qpart = qpart_choices[qpart_selected]
#times = [0.0 for i in 1:10]
#for i in 1:10
rezz = main(flux_type,n,p,mc_steps,particles,step_size,rad_choice,qpart,log_form,0,with_berry,allowed_sets_matrix,full_pasc_tri,full_derivs)
#println("Berry Phase = ",mean(rezz[4])," +/- ",std(rezz[4]))
#println("$filling, Calced Filling = ",2*mean(rezz[4])/x_rad^2," +/- ",2*std(rezz[4])/x_rad^2)
#println("Time = ",(rezz[5]/1000)/3600)
#times[i] = rezz[5]/10000
#println(times[i])
#end
#println(mean(times)," +/- ",std(times))
#println("Time for 50,000 Samples: ",50000*mean(times)/3600," +/- ",50000*std(times)/3600)


#compl = collect(Iterators.flatten([rezz[2][i,:] for i in 1:particles]))
#figure()
#hist2D(real.(compl),-imag.(compl),bins=50)



"fin"
