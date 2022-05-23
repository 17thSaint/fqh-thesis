#import Pkg; Pkg.add("LinearAlgebra")
using LinearAlgebra

include("write-accmat-hdf5.jl")

function get_energies_alltimes(xs,ys,m,qhole,num_parts,dt)
	time_counts = Int(size(xs)[2]/dt)
	energies = [0.0 for i in 1:time_counts]
	for i in 1:time_counts
		#if i % (0.05*time_counts) == 0
		#	println("Getting Energies: ",(i-1)*100/time_counts,"%")
		#end
		time = Int(1+(i-1)*dt)
		input_config = 
		local_energy = prob_wavefunc(input_config,m,qhole,num_parts)
		energies[i] = local_energy
	end
	return energies
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
	len = 1000
	dts = [1+(i-1)*1 for i in 1:Int(0.1*len)-1]
	autocorr = [0.0 for i in 1:Int(0.1*len)-1]
	for i in 1:Int(0.1*len)-1
		autocorr[i] = auto_correlation(energy,dts[i])
	end
	check_below_tol = autocorr .< [0.01 for i in 1:length(autocorr)]
	#println(autocorr)
	corr_length = samp_freq*dts[findall(check_below_tol)[1]]
	return corr_length,dts,autocorr
end

function start_rand_config(num_parts::Int,n::Int,p::Int)
	filling = n/(2*p*n+1)
	rm = sqrt(2*num_parts/filling)
	config = [rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm for i in 1:num_parts]
	return config
end

function get_wf_elem(num_parts,element,n,qpart=[0,[0]])
	qpart_shift = qpart[1]
	if n == 2 && element[1] > floor(num_parts/2) + qpart_shift
		m = 1
		pow = element[1] - floor(num_parts/2) - 1 - qpart_shift
		result = [m,pow]
	else
		m = 0
		pow = element[1] - 1 - qpart_shift
		result = [m,pow]	
	end
	return result
end

function get_qpart_wf_exp(config,qpart,which_qpart,part,n,exp_check=true)
	lstar2 = 2*1*n + 1
	position = config[part]
	if exp_check
		result = exp(conj(qpart[2][which_qpart])*position/(2*lstar2) - abs2(qpart[2][which_qpart])/(4*lstar2))
	else
		result = conj(qpart[2][which_qpart])*position/(2*lstar2) - abs2(qpart[2][which_qpart])/(4*lstar2)
	end
	return result
end


function get_Jis(config::Vector{ComplexF64},part::Int,acc_parts=[false,[]])
	ji = 1.0
	num_parts = length(config)
	if !acc_parts[1]
		acc_parts[2] = deleteat!([i for i in 1:num_parts],[i for i in 1:num_parts] .== part)
	end
	for i in acc_parts[2]
		dist_btw = config[part]-config[i]
		ji *= dist_btw
	end
	return Complex(ji)
end


function get_Jiprime(config,part,p)
	jiprime = 0
	num_parts = length(config)
	for j in 1:num_parts
		if j == part
			continue
		end
		pos_selected = config[part]
		jiprime_local = p*(pos_selected - config[j])^(p-1)
		for k in 1:num_parts
			if k == part || k == j
				continue
			end
			dist_btw = pos_selected - config[k]
			jiprime_local *= dist_btw^p
		end
		jiprime += jiprime_local
	end
	#=
	if real(qhole[1]) > 0.1
		jiprime *= config[part] - qhole[2]
		jiprime += get_Jis(config,part,qhole)
	end
	=#
	return jiprime
end

function get_Ji2prime(config,part,p)
	ji2prime = 0
	num_parts = length(config)
	pos_selected = config[part]
	for j in 1:num_parts
		if j == part
			continue
		end
		ji2prime_local = 0
		for k in 1:num_parts
			if k == part || k == j
				continue
			end
			ji2prime_local_local = p*(pos_selected - config[k])^(p-1)
			for m in 1:num_parts
				if m == part || m == j || m == k
					continue
				end
				dist_btw = pos_selected - config[m]
				ji2prime_local_local *= dist_btw^p
			end
			ji2prime_local += ji2prime_local_local
		end
		ji2prime += ji2prime_local*p*(pos_selected - config[j])^(p-1)
	end
	
	if p > 1
		second_part = 0.0+im*0.0
		for j in 1:num_parts
			if j == part
				continue
			end
			ji2prime_second_local = p*(p-1)*(pos_selected - config[j])^(p-2)
			for k in 1:num_parts
				if k == part || k == j
					continue
				end
				dist_btw = pos_selected - config[k]
				next = dist_btw^p
				ji2prime_second_local *= next
			end
			second_part += ji2prime_second_local
		end
		ji2prime += second_part
	end
	
	return ji2prime
end

function get_elem_projection(config,part,row,n,p,qpart=[0,[0]])
	num_parts = length(config)
	matrix_element = [row,part]
	qpart_shift = qpart[1]
	Ji, Jiprime, Ji2prime = get_Jis(config,part)^p,get_Jiprime(config,part,p),get_Ji2prime(config,part,p)
	if row >= qpart_shift + 1
		bar, l = get_wf_elem(num_parts,matrix_element,n)
		if bar > 0
			result = 2*(l*(config[part]^(l-1))*Ji + (config[part]^l)*Jiprime)
		else
			result = (config[part]^l)*Ji
		end
	elseif row < qpart_shift + 1
		lstar2 = 2*1*n + 1
		coeff = [(1/lstar2) - 1,1/lstar2^2 - 2\lstar2 + 1]
		coeff_jp = [2,-4*conj(qpart[2][row])]
		coeff_jpp = [0,4]
		exp_part = get_qpart_wf_exp(config,qpart,row,part,n)
		result = coeff[n]*(conj(qpart[2][row])^n)*exp_part*Ji + coeff_jp[n]*Jiprime*exp_part + coeff_jpp[n]*exp_part*Ji2prime
	end
	return result
end

function get_wavefunc(config,n,p,qpart=[0,[0]])
	num_parts = length(config)
	matrix_full = fill(0.0+0.0*im,(num_parts,num_parts))
	wavefunc = 1.0
	for i in 1:num_parts
		if i <= num_parts
			wavefunc *= exp(-abs2(config[i])/4)
		end
		for j in 1:num_parts
			dats = get_elem_projection(config,j,i,n,p,qpart)
			matrix_full[i,j] = dats
		end
	end
	#println(matrix_full)
	wavefunc *= det(matrix_full)
	return wavefunc#,matrix_full
end

function get_logJi(config::Vector{ComplexF64},part::Int64,acc_parts=[false,[]],shift=0.0)
	logji::ComplexF64 = 0.0
	num_parts::Int = length(config)
	if !acc_parts[1]
		acc_parts[2] = deleteat!([i for i in 1:num_parts],[i for i in 1:num_parts] .== part)
	end
	for i in acc_parts[2]
		dist_btw::ComplexF64 = config[part]+shift-config[i]
		logji += log(Complex(dist_btw))
	end
	return Complex(logji)
end

function get_logJastrowfull(config)
	num_parts = length(config)
	logJastrow = 0.0
	for i in 1:num_parts
		logJastrow += get_logJi(config,i)
	end
	return logJastrow
end

function get_log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
end

function get_logJiprime(config,part,p)
	logjiprime = 0.0+im*0.0
	num_parts = length(config)
	start = 0
	for i in 1:num_parts
		if i == part
			continue
		end
		next_part = log(Complex(p)) + (p-1)*log(Complex(config[part] - config[i]))
		for j in 1:num_parts
			if j == i || j == part
				continue
			end
			dist_btw = config[part] - config[j]
			next_part += p*log(Complex(dist_btw))
		end
		if start == 0
			logjiprime += next_part
		else
			logjiprime = get_log_add(next_part,logjiprime)
		end
		start += 1
	end
	return logjiprime
end	

function get_logJi2prime(config,part,p)
	num_parts = length(config)
	logjiprime = get_logJiprime(config,part,p)
	logji2prime = 0.0+im*0.0
	start = 0
	for i in 1:num_parts
		if i == part
			continue
		end
		start2 = 0
		next_part = 0.0
		for j in 1:num_parts
			if j == i || j == part
				continue
			end
			next_next_part = log(p) + (p-1)*log(Complex(config[part] - config[j]))
			for m in 1:num_parts
				if m == i || m == part || m == j
					continue
				end
				dist_btw = config[part] - config[m]
				next_next_part += p*log(Complex(dist_btw))
			end
			if start2 == 0
				next_part += next_next_part
			else
				local_next_part = next_part + 0.0
				next_part = get_log_add(local_next_part,next_next_part)
			end
			start2 += 1
		end
		if start == 0
			logji2prime += next_part + log(p) + (p-1)*log(Complex(config[part] - config[i]))
			start += 1
		else
			logji2prime = get_log_add(next_part + log(p) + (p-1)*log(Complex(config[part] - config[i])),logji2prime)
		end
	end
	
	if p > 1
		start_second = 0
		second_part = 0.0+im*0.0
		for i in 1:num_parts
			if i == part
				continue
			end
			second_next_part = log(Complex(p)) + log(Complex(p-1)) + (p-2)*log(Complex(config[part] - config[i]))
			for j in 1:num_parts
				if j == i || j == part
					continue
				end
				dist_btw = config[part] - config[j]
				next = p*log(Complex(dist_btw))
				second_next_part += next
			end
			if start_second == 0
				second_part += second_next_part
			else
				local_second = second_part + 0.0
				second_part = get_log_add(second_next_part,local_second)
			end
			start_second += 1
		end
		logji2prime = get_log_add(logji2prime,second_part)
	end
	
	return logji2prime
end	


function get_log_elem_proj(config,part,row,n,p,qpart=[0,[0]])
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element,n,qpart)
	qpart_shift = qpart[1]
	lstar2 = 2*1*n + 1
	coeff = [(1/lstar2) - 1,1/lstar2^2 - 2\lstar2 + 1]
	if row >= qpart_shift + 1
		if bar > 0
			
			if l == 0
				logJiprime = get_logJiprime(config,part,p)
				result = log(2) + logJiprime
			else
				logJi = p*get_logJi(config,part)
				logJiprime = get_logJiprime(config,part,p)
				a = log(2) + log(l) + (l-1)*log(config[part]) + logJi
				b = log(2) + l*log(config[part]) + logJiprime
				result = get_log_add(a,b)
			end
		else
			if l == 0
				logJi = p*get_logJi(config,part)
				result = logJi
			else
				logJi = p*get_logJi(config,part)
				result = logJi + l*log(config[part])
			end
		end
		
	else
		logetabar = log(conj(qpart[2][row]))
		log_exppart = get_qpart_wf_exp(config,qpart,row,part,n,false)
		if n == 1
			logJi = p*get_logJi(config,part)
			logJiprime = get_logJiprime(config,part,p)
			part1 = logetabar + log(Complex(coeff[n])) + logJi
			part2 = log(2) + logJiprime
			result = log_exppart + get_log_add(part1,part2)
		else
			logJi = p*get_logJi(config,part)
			logJiprime = get_logJiprime(config,part,p)
			logJi2prime = get_logJi2prime(config,part,p)
			part1 = 2*logetabar + logJi + log(Complex(coeff[n]))
			part2 = log(4) + logetabar + logJiprime
			part3 = log(4) + logJi2prime
			result = log_exppart + get_log_add(part1,get_log_add(part2,part3))
		end
	end
	return result
end

function get_diag_log_det(matrix)
	num_parts = size(matrix)[1]
	active_matrix = matrix + fill(0.0,(num_parts,num_parts))
	outsides = [active_matrix[i,i] for i in 1:num_parts]
	for i in 1:num_parts
		active_matrix[:,i] .-= outsides[i]
	end
	logdet = sum(outsides) + log(Complex(det(exp.(active_matrix))))
	return logdet,outsides,active_matrix
end

function get_log_det(matrix::Matrix{ComplexF64},reg_input=false)
	num_parts::Int64 = size(matrix)[1]
	maxes::Vector{ComplexF64} = [0.0+0.0*im for i in 1:num_parts]
	rejected_indices::Vector{Int64} = []
	
	changed::Matrix{ComplexF64} = matrix + fill(0.0,(num_parts,num_parts))
	for i in 1:num_parts
		row::Vector{Float64} = real(changed[i,:])
		validation::Bool = true
		j::Int64 = 0
		index::Int64 = 0
		while validation
			val::ComplexF64 = sort(row)[end - j]
			guess_index::Int64 = findfirst(real.(changed[i,:]) .== val )
			if any(rejected_indices .== index)
				j += 1
			else
				index = guess_index
				validation = false
			end
		end
		maxes[i] = changed[i,index]
		changed[i,:] .-= maxes[i]
		append!(rejected_indices,index)
	end
	reduced_logdet::ComplexF64 = sum(maxes) + log(Complex(det(exp.(changed))))

	return reduced_logdet,maxes,changed
end


function get_wavefunc_fromlog(config,n,p,qpart=[0,[0]])
	num_parts = length(config)
	log_matrix = fill(0.0+im*0.0,(num_parts,num_parts))
	log_matrix_sepJi = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts
			data_here = get_log_elem_proj(config,j,i,n,p,qpart)
			
			if isnan(data_here)
				println("NaN: ",j,", ",i)
			#	break
			end
			log_matrix[i,j] = data_here[1]
		end
	end
	
	result = get_log_det(log_matrix)[1]
	
	for i in 1:num_parts
		result += -abs2(config[i])/4
	end
	
	return result#,log_matrix
end

function dist_btw_Laugh(part_1,part_2)
	return sqrt((part_1[1] - part_2[1])^2 + (part_1[2] - part_2[2])^2)
end

function prob_wavefunc_laughlin(complex_config, m)
	full = 0
	num_parts = length(complex_config)
	config = [[real(complex_config[i]),imag(complex_config[i])] for i in 1:num_parts]
	for j = 1:num_parts
		for i in 1:num_parts
			if i == j
				continue
			end
			dist = dist_btw_Laugh(config[i],config[j]) 
			full += -m*log( dist )
		end
		
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)
	end
	return full
end

function nested_loop(loop_level::Int64,allowed_vals_dict::Dict{String,Vector{Any}},all_ji_acc_sets::Vector)
	#sum_parts = [0 for i in 1:order-1]
	order::Int = length(keys(allowed_vals_dict))
	parts_count::Int = length(allowed_vals_dict["s1"]) + 1
	if loop_level == order
		for i in 1:length(allowed_vals_dict["s$order"])
			ji_allowed_vals = deleteat!(allowed_vals_dict["s$order"].+(1-1),i)
			
			#=
			sum_parts = [0 for k in 1:order]
			for k in 1:order-1
				next = k+1
				sum_parts[k] = deleteat!(allowed_vals_dict["s$k"].+(1-1),findall(x->x in allowed_vals_dict["s$next"],allowed_vals_dict["s$k"]))[1]
			end
			sum_parts[end] = deleteat!(allowed_vals_dict["s$order"].+(1-1),findall(x->x in ji_allowed_vals,allowed_vals_dict["s$order"]))[1]
			println(i,", Sum Parts: $sum_parts","Jis: ",ji_allowed_vals)
			=#
			
			append!(all_ji_acc_sets,[ji_allowed_vals])
		end
		
		return 
	
	end
	next_level::Int = loop_level + 1
	for i in 1:length(allowed_vals_dict["s$loop_level"])
		sum_part = allowed_vals_dict["s$loop_level"][i]
		allowed_vals_dict["s$next_level"] = deleteat!(allowed_vals_dict["s$loop_level"].+(1-1),i)
		nested_loop(loop_level + 1,allowed_vals_dict,all_ji_acc_sets)
	end
	
end

function get_all_acc_sets(order::Int64,part::Int64,parts_count::Int64)
	starting_allowed_vals_dict::Dict{String,Vector{Any}} = Dict([("s$i",[]) for i in 1:order])
	starting_allowed_vals_dict["s1"] = deleteat!([i for i in 1:parts_count],[i for i in 1:parts_count] .== part)
	
	all_ji_acc_sets::Vector{Vector{Int64}} = []
	nested_loop(1,starting_allowed_vals_dict,all_ji_acc_sets)
	return all_ji_acc_sets
end

function compress_acc_set(parts_count::Int64,acc_set::Vector{Vector{Int64}})
	new_acc_set = [ [1,[1 for j in 1:length(acc_set[1])]] for i in 1:binomial(parts_count-1,length(acc_set[1]))]
	index = 1
	for i in unique(acc_set)
		found = findall(x->x==i,acc_set)
		new_acc_set[index] = [length(found),i]
		index += 1
	end
	return new_acc_set
end

function get_allowed_sets_matrix(num_parts::Int64)
	allowed_sets_matrix = Matrix{Any}(undef,num_parts,num_parts-1)
	for which_part in 1:num_parts
		for which_order in 1:num_parts-1
			allowed_sets_matrix[which_part,which_order] = compress_acc_set(num_parts,get_all_acc_sets(which_order,which_part,num_parts))
		end
	end
	
	return allowed_sets_matrix
end

function get_nested_logadd(loop_level::Int64,all_vals::Vector{ComplexF64},result::ComplexF64)
	if loop_level == 1
		#println(result,", ",loop_level)
		return Complex(result)
	else
		result = get_log_add(result,all_vals[loop_level-1])
		get_nested_logadd(loop_level - 1,all_vals,result)
	end
end

function split_nested_logadd(all_vals::Vector{ComplexF64})
	full_length::Int64 = length(all_vals)
	max_length::Int64 = 49999
	remainder_length::Int64 = full_length%max_length
	count_max_length::Int64 = Int(floor(full_length/max_length))
	delineated::Vector{Vector{ComplexF64}} = [ all_vals[Int(max_length*(i-1)+1):Int(max_length*i)] for i in 1:count_max_length]
	append!(delineated,[all_vals[full_length-remainder_length+1:full_length]])
	
	each_subnest_logadded::Vector{ComplexF64} = [get_nested_logadd(length(delineated[i]), delineated[i], delineated[i][end]) for i in 1:count_max_length+1]
	
	result_allnests = get_nested_logadd(count_max_length+1, each_subnest_logadded, each_subnest_logadded[end])
	return Complex(result_allnests)
end

function get_nth_deriv_Ji(config::Vector{ComplexF64},part::Int64,acc_sets_column::Vector{Any},log_form=false)
	parts_count::Int64 = length(config)
	result::ComplexF64 = 0.0+im*0.0
	all_acc_sets::Vector{Any} = acc_sets_column[part]
	if !log_form
		all_jis::Vector{ComplexF64} = [all_acc_sets[i][1]*get_Jis(config,part,[true,all_acc_sets[i][2]]) for i in 1:length(all_acc_sets)]
		result = sum(all_jis)
	else
		all_jis = [get_logJi(config,part,[true,all_acc_sets[i][2]]) + log(all_acc_sets[i][1]) for i in 1:length(all_acc_sets)]
		#println(length(all_jis))
		if any(isinf.(all_jis))
			result = -Inf
		else
			if length(all_jis) > 50000
				result = split_nested_logadd(all_jis)
			else
				result = get_nested_logadd(length(all_jis),all_jis,all_jis[end]+1-1)
			end
		end
	end
	
	return result
end

function get_nth_deriv_Ji_old(config::Vector{ComplexF64},part::Int64,acc_sets_column::Vector{Vector{Vector{Int64}}},log_form=false)
	parts_count::Int64 = length(config)
	result::ComplexF64 = 0.0+im*0.0
	all_acc_sets::Vector{Vector{Int64}} = acc_sets_column[part]
	if !log_form
		all_jis::Vector{ComplexF64} = [get_Jis(config,part,[true,all_acc_sets[i]]) for i in 1:length(all_acc_sets)]
		result = sum(all_jis)
	else
		all_jis = [get_logJi(config,part,[true,all_acc_sets[i]]) for i in 1:length(all_acc_sets)]
		#println(length(all_jis))
		if any(isinf.(all_jis))
			result = -Inf
		else
			if length(all_jis) > 50000
				result = split_nested_logadd(all_jis)
			else
				result = get_nested_logadd(length(all_jis),all_jis,all_jis[end]+1-1)
			end
		end
	end
	
	return result
end

function get_nth_deriv_Ji_columncall(config::Vector{ComplexF64},part::Int64,all_acc_sets::Vector{Vector{Int64}},log_form=false)
	parts_count::Int64 = length(config)
	result::ComplexF64 = 0.0+im*0.0
	if !log_form
		all_jis::Vector{ComplexF64} = [get_Jis(config,part,[true,all_acc_sets[i]]) for i in 1:length(all_acc_sets)]
		result = sum(all_jis)
	else
		all_jis = [get_logJi(config,part,[true,all_acc_sets[i]]) for i in 1:length(all_acc_sets)]
		#println(length(all_jis))
		if any(isinf.(all_jis))
			result = -Inf
		else
			if length(all_jis) > 50000
				result = split_nested_logadd(all_jis)
			else
				result = get_nested_logadd(length(all_jis),all_jis,all_jis[end]+1-1)
			end
		end
	end
	
	return result
end

function get_pascals_triangle(n::Int64)
	n == 0 && return [],[0]
	n == 1 && return [[1]],[1]
	t::Vector{Vector{Int64}} = get_pascals_triangle(n-1)[1]
	push!(t, [t[end];0] + [0;t[end]])
	row::Vector{Int64} = t[end]
	folded::Vector{Int64} = [row[i] + row[end-i+1] for i in 1:Int(floor(n/2))]
	n%2 != 0.0 && append!(folded,[ row[Int(ceil(n/2))] ])
	return t,folded
end

function get_rf_elem_proj(config::Vector{ComplexF64},part::Int64,row::Int64,acc_sets_matrix::Matrix{Any},pascals_row::Vector{Int64},deriv_orders::Vector{Vector{Int64}},qpart=[0,[0]],log_form=false)
	lstar::Float64 = sqrt(2*1*1-1)
	qpart_shift::Int64 = qpart[1]
	if !log_form
		if row >= qpart_shift + 1
			jis::Vector{ComplexF64} = [get_Jis(config,part)]
			append!(jis,[get_nth_deriv_Ji(config,part,acc_sets_matrix[:,i]) for i in 1:row-1])
			#string_result = join([string(pascals_row[i],"J(",deriv_orders[i][1]-1,"')J(",deriv_orders[i][2]-1,"')") for i in 1:length(pascals_row)],"+")
			indiv_terms::Vector{ComplexF64} = [pascals_row[i]*jis[deriv_orders[i][1]]*jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]
			result::ComplexF64 = sum(indiv_terms)
			if row > 1
				result *= 2
			end
		else
			shift_part::ComplexF64 = conj(qpart[2][row])/(lstar^2)
			ji_shifted::ComplexF64 = get_Jis(config.-shift_part,part)
			front_term::ComplexF64 = config[part] + shift_part - conj(qpart[2][row])
			result = front_term*exp(-shift_part/4)*ji_shifted^2
			#string_result = latexstring("\$ exp()(z_$part - \\bar{\\eta_{$row}})J_{S$part}^2 \$")
		end
	else
		if row >= qpart_shift + 1
			jis = [get_logJi(config,part)]
			append!(jis,[get_nth_deriv_Ji(config,part,acc_sets_matrix[:,i],log_form) for i in 1:row-1])
			indiv_terms = [log(pascals_row[i]) + jis[deriv_orders[i][1]] + jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]  # get rid of for loop here
			if any(isinf.(indiv_terms))
				result = -Inf
			else
				result = get_nested_logadd(length(indiv_terms),indiv_terms,indiv_terms[end]+1-1)
				if row > 1
					result += log(2)
				end
			end
			if isnan(result)
				println("Nan")
				return indiv_terms
			end
		else
			# this is RFA version
			shift_part = conj(qpart[2][row])/(lstar^2)
			ji_shifted = get_logJi(config,part,[false,[]],shift_part)
			front_term = log(Complex(config[part] + shift_part - conj(qpart[2][row])))
			result = front_term - shift_part*qpart[2][row]/4 + 2*ji_shifted
			#
			#= 
			exp_part = (conj(qpart[2][row])*config[part] - abs2(qpart[2][row])/2)/(2*lstar^2)
			front_term = config[part] - conj(qpart[2][row])
			result = log(Complex(front_term))# + exp_part
			=#
		end
	end
	return result
end

function get_rf_elem_proj_old(config::Vector{ComplexF64},part::Int64,row::Int64,acc_sets_matrix::Matrix{Vector{Vector{Int64}}},pascals_row::Vector{Int64},deriv_orders::Vector{Vector{Int64}},qpart=[0,[0]],log_form=false)
	lstar::Float64 = sqrt(2*1*1-1)
	qpart_shift::Int64 = qpart[1]
	if !log_form
		if row >= qpart_shift + 1
			jis::Vector{ComplexF64} = [get_Jis(config,part)]
			append!(jis,[get_nth_deriv_Ji_old(config,part,acc_sets_matrix[:,i]) for i in 1:row-1])
			#string_result = join([string(pascals_row[i],"J(",deriv_orders[i][1]-1,"')J(",deriv_orders[i][2]-1,"')") for i in 1:length(pascals_row)],"+")
			indiv_terms::Vector{ComplexF64} = [pascals_row[i]*jis[deriv_orders[i][1]]*jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]
			result::ComplexF64 = sum(indiv_terms)
			if row > 1
				result *= 2
			end
		else
			shift_part::ComplexF64 = conj(qpart[2][row])/(lstar^2)
			ji_shifted::ComplexF64 = get_Jis(config.-shift_part,part)
			front_term::ComplexF64 = config[part] + shift_part - conj(qpart[2][row])
			result = front_term*exp(-shift_part/4)*ji_shifted^2
			#string_result = latexstring("\$ exp()(z_$part - \\bar{\\eta_{$row}})J_{S$part}^2 \$")
		end
	else
		if row >= qpart_shift + 1
			jis = [get_logJi(config,part)]
			append!(jis,[get_nth_deriv_Ji_old(config,part,acc_sets_matrix[:,i],log_form) for i in 1:row-1])
			indiv_terms = [log(pascals_row[i]) + jis[deriv_orders[i][1]] + jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]  # get rid of for loop here
			if any(isinf.(indiv_terms))
				result = -Inf
			else
				result = get_nested_logadd(length(indiv_terms),indiv_terms,indiv_terms[end]+1-1)
				if row > 1
					result += log(2)
				end
			end
			if isnan(result)
				println("Nan")
				return indiv_terms
			end
		else
			# this is RFA version
			shift_part = conj(qpart[2][row])/(lstar^2)
			ji_shifted = get_logJi(config,part,[false,[]],shift_part)
			front_term = log(Complex(config[part] + shift_part - conj(qpart[2][row])))
			result = front_term - shift_part*qpart[2][row]/4 + 2*ji_shifted
			#
			#= 
			exp_part = (conj(qpart[2][row])*config[part] - abs2(qpart[2][row])/2)/(2*lstar^2)
			front_term = config[part] - conj(qpart[2][row])
			result = log(Complex(front_term))# + exp_part
			=#
		end
	end
	return result
end

function get_rf_elem_proj_columncall(config::Vector{ComplexF64},part::Int64,row::Int64,pascals_row::Vector{Int64},deriv_orders::Vector{Vector{Int64}},qpart=[0,[0]],log_form=false)
	all_parts = [i for i in 1:length(config)]
	lstar::Float64 = sqrt(2*1*1+1)
	qpart_shift::Int64 = qpart[1]
	if !log_form
		if row >= qpart_shift + 1
			allsets = read_acc_matrix_data("acc-matrix-data",length(config),part,[i for i in 1:row-1])
			jis::Vector{ComplexF64} = [get_Jis(config,part)]
			append!(jis,[get_nth_deriv_Ji_columncall(config,part,allsets[i]) for i in 1:row-1])
			#string_result = join([string(tri_coeffs[i],"J(",deriv_orders[i][1]-1,"')J(",deriv_orders[i][2]-1,"')") for i in 1:length(tri_coeffs)],"+")
			indiv_terms::Vector{ComplexF64} = [pascals_row[i]*jis[deriv_orders[i][1]]*jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]
			result::ComplexF64 = 2*sum(indiv_terms)
		else
			shift_part::ComplexF64 = conj(qpart[2][row])/(lstar^2)
			ji_shifted::ComplexF64 = get_Jis(config.-shift_part,part)
			front_term::ComplexF64 = config[part] + shift_part - conj(qpart[2][row])
			result = front_term*exp(-shift_part/4)*ji_shifted^2
			#string_result = latexstring("\$ exp()(z_$part - \\bar{\\eta_{$row}})J_{S$part}^2 \$")
		end
	else
		if row >= qpart_shift + 1
			allsets = read_acc_matrix_data("acc-matrix-data",length(config),part,[i for i in 1:row-1])
			jis = [get_logJi(config,part)]
			append!(jis,[get_nth_deriv_Ji_columncall(config,part,allsets[i],log_form) for i in 1:row-1])
			#println(part," ,",row,", ",pascals_row)
			indiv_terms = [log(pascals_row[i]) + jis[deriv_orders[i][1]] + jis[deriv_orders[i][2]] for i in 1:length(pascals_row)]  # get rid of for loop here
			if any(isinf.(indiv_terms))
				result = -Inf
			else
				result = log(2) + get_nested_logadd(length(indiv_terms),indiv_terms,indiv_terms[end]+1-1)
			end
			if isnan(result)
				println("Nan")
				return indiv_terms
			end
		else
			shift_part = conj(qpart[2][row])/(lstar^2)
			ji_shifted = get_logJi(config.-shift_part,part)
			front_term = log(Complex(config[part] + shift_part - conj(qpart[2][row])))
			result = front_term - shift_part/4 + 2*ji_shifted
		end
	end
	return result
end

function get_deriv_orders_matrix(num_parts::Int64)
	all_deriv_ords = Vector{Vector{Vector{Int64}}}(undef,num_parts)
	for j_rows in 1:num_parts
		len::Int64 = Int(ceil(j_rows/2))
		all_deriv_ords[j_rows] = [[j_rows-k+1,k] for k in 1:len]
	end
	return all_deriv_ords
end

function get_rf_wavefunc_old(config::Vector{ComplexF64},acc_sets_matrix::Matrix{Vector{Vector{Int64}}},all_pascal::Vector{Vector{Int64}},all_deriv_orders::Vector{Vector{Vector{Int64}}},qpart=[0,[0]],log_form=false)
	num_parts::Int8 = length(config)
	wavefunc::ComplexF64 = 1.0
	if log_form
		wavefunc = 0.0
	end
	full_matrix::Matrix{ComplexF64} = fill(0.0+im*0.0,(num_parts,num_parts))
	#string_matrix = fill("",(num_parts,num_parts))
	for i in 1:num_parts
		#
		pascal_row::Vector{Int64} = all_pascal[i]
		deriv_orders::Vector{Vector{Int64}} = all_deriv_orders[i]
		if i <= num_parts
			#
			if !log_form
				wavefunc *= exp(-abs2(config[i])/4)
			else
				wavefunc -= abs2(config[i])/4
			end
			#
		end
		#
		for j in 1:num_parts
			full_matrix[i,j] = get_rf_elem_proj_old(config,j,i,acc_sets_matrix,pascal_row,deriv_orders,qpart,log_form)
			#string_matrix[i,j] = get_rf_elem_proj(config,j,i,acc_sets_matrix,pascal_row,deriv_orders,qpart,log_form)
		end
	end
	#
	if !log_form
		mat_det::ComplexF64 = det(full_matrix)
		wavefunc *= mat_det
	else
		mat_det = get_log_det(full_matrix)[1]
		wavefunc += mat_det
	end
	#
	return wavefunc#, full_matrix
end

function get_rf_wavefunc(config::Vector{ComplexF64},acc_sets_matrix::Matrix{Any},all_pascal::Vector{Vector{Int64}},all_deriv_orders::Vector{Vector{Vector{Int64}}},qpart=[0,[0]],log_form=false)
	num_parts::Int8 = length(config)
	wavefunc::ComplexF64 = 1.0
	if log_form
		wavefunc = 0.0
	end
	full_matrix::Matrix{ComplexF64} = fill(0.0+im*0.0,(num_parts,num_parts))
	#string_matrix = fill("",(num_parts,num_parts))
	for i in 1:num_parts
		#
		pascal_row::Vector{Int64} = all_pascal[i]
		deriv_orders::Vector{Vector{Int64}} = all_deriv_orders[i]
		if i <= num_parts
			#
			if !log_form
				wavefunc *= exp(-abs2(config[i])/4)
			else
				wavefunc -= abs2(config[i])/4
			end
			#
		end
		#
		for j in 1:num_parts
			full_matrix[i,j] = get_rf_elem_proj(config,j,i,acc_sets_matrix,pascal_row,deriv_orders,qpart,log_form)
			#string_matrix[i,j] = get_rf_elem_proj(config,j,i,acc_sets_matrix,pascal_row,deriv_orders,qpart,log_form)
		end
	end
	#
	if !log_form
		mat_det::ComplexF64 = det(full_matrix)
		wavefunc *= mat_det
	else
		mat_det = get_log_det(full_matrix)[1]
		wavefunc += mat_det
	end
	#
	return wavefunc#, full_matrix
end



function get_rf_wavefunc_columncall(config::Vector{ComplexF64},all_pascal::Vector{Vector{Int64}},all_deriv_orders::Vector{Vector{Vector{Int64}}},qpart=[0,[0]],log_form=false)
	num_parts::Int8 = length(config)
	wavefunc::ComplexF64 = 1.0
	if log_form
		wavefunc = 0.0
	end
	full_matrix::Matrix{ComplexF64} = fill(0.0+im*0.0,(num_parts,num_parts))
	#string_matrix = fill("",(num_parts,num_parts))
	for i in 1:num_parts
		#
		pascal_row::Vector{Int64} = all_pascal[i]
		deriv_orders::Vector{Vector{Int64}} = all_deriv_orders[i]
		if i <= num_parts
			#
			if !log_form
				wavefunc *= exp(-abs2(config[i])/4)
			else
				wavefunc -= abs2(config[i])/4
			end
			#
		end
		#
		for j in 1:num_parts
			full_matrix[i,j] = get_rf_elem_proj_columncall(config,j,i,pascal_row,deriv_orders,qpart,log_form)
			#string_matrix[i,j] = get_rf_elem_proj(config,j,i,acc_sets_matrix,pascal_row,qpart)
		end
	end
	#
	if !log_form
		mat_det::ComplexF64 = det(full_matrix)
		wavefunc *= mat_det
	else
		mat_det = get_log_det(full_matrix)[1]
		wavefunc += mat_det
	end
	#
	return wavefunc
end






"fin"
