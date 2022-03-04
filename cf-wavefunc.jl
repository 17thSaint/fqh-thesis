using LinearAlgebra,PyPlot

function get_wf_elem(num_parts,element,n,qpart=[0,[0]])
	qpart_shift = qpart[1]
	if n == 2 && element[1] > num_parts/2 + qpart_shift
		m = 1
		pow = element[1] - num_parts/2 - 1 - qpart_shift
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

function get_Jis(config,part,p)
	ji = 1.0
	num_parts = length(config)
	for i in 1:num_parts
		if i == part
			continue
		end
		dist_btw = config[part]-config[i]
		ji *= dist_btw^p
	end
	#=
	if real(qhole[1]) > 0.1
		ji *= config[part] - qhole[2]
	end
	=#
	return ji
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
				#println("Reg $k $part: ",next)
			end
			second_part += ji2prime_second_local
			println("Reg 2: ",log(second_part))
		end
		ji2prime += second_part
	end
	
	return ji2prime
end

function get_elem_projection(config,part,row,n,p,qpart=[0,[0]])
	num_parts = length(config)
	matrix_element = [row,part]
	qpart_shift = qpart[1]
	Ji, Jiprime, Ji2prime = get_Jis(config,part,p),get_Jiprime(config,part,p),get_Ji2prime(config,part,p)
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
	wavefunc *= det(matrix_full)
	return wavefunc
end

function get_logJi(config,part,p)
	logji = 0.0
	num_parts = length(config)
	for i in 1:num_parts
		if i == part
			continue
		end
		dist_btw = config[part]-config[i]
		logji += p*log(Complex(dist_btw))
	end
	return logji
end

function get_log_add(a,b)
	if real(a) > real(b)
		ordered = [b,a]
	else
		ordered = [a,b]
	end
	result = ordered[2] + log(1 + exp(ordered[1] - ordered[2]))
	return result
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
		next_part = log(p) + (p-1)*log(Complex(config[part] - config[i]))
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
				next_part = next_next_part
			else
				next_part = get_log_add(next_part,next_next_part)
			end
			start2 += 1
		end
		if start == 0
			logji2prime += next_part
			start += 1
		else
			logji2prime = get_log_add(next_part,logji2prime)
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
				#println("Log $j $part: ",exp(next))
			end
			if start_second == 0
				second_part += second_next_part
			else
				local_second = second_part + 0.0
				second_part = get_log_add(second_next_part,local_second)
			end
			start_second += 1
			println("Log 2: ",second_part)
		end
		logji2prime = get_log_add(logji2prime,second_part)
	end
	
	return logji2prime
end	


function get_log_elem_proj(config,part,row,n,p,qpart=[0,[0]])
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element,n,qpart)
	logJi, logJiprime, logJi2prime = get_logJi(config,part,p), get_logJiprime(config,part,p), get_logJi2prime(config,part,p)
	qpart_shift = qpart[1]
	lstar2 = 2*1*n + 1
	coeff = [(1/lstar2) - 1,1/lstar2^2 - 2\lstar2 + 1]
	if row >= qpart_shift + 1
		if bar > 0
			if l == 0
				result = log(2) + logJiprime
			else
				a = log(2) + log(l) + (l-1)*log(config[part]) + logJi
				b = log(2) + l*log(config[part]) + logJiprime
				result = get_log_add(a,b)
			end
		else
			result = l*log(config[part]) + logJi
		end
	else
		logetabar = log(conj(qpart[2][row]))
		log_exppart = get_qpart_wf_exp(config,qpart,row,part,n,false)
		if n == 1
			part1 = logetabar + logJi + log(Complex(coeff[n]))
			part2 = log(2) + logJiprime
			result = log_exppart + get_log_add(part1,part2)
		else
			part1 = 2*logetabar + logJi + log(Complex(coeff[n]))
			part2 = log(4) + logetabar + logJiprime
			part3 = log(4) + logJi2prime
			result = log_exppart + get_log_add(part1,get_log_add(part2,part3))
		end
	end
	return result
end

function get_log_det(matrix,reg_input=false)
	num_parts = size(matrix)[1]
	maxes = [0.0+0.0*im for i in 1:num_parts]
	rejected_indices = []
	starting = matrix
	#expected_reg_det = det(exp.(matrix))
	if reg_input
		#expected_reg_det = det(matrix)
		starting = log.(matrix)
	end
	changed = starting + fill(0.0,(num_parts,num_parts))
	for i in 1:num_parts
		allowed_indices = [m for m in 1:num_parts]
		overlap = [findall(q->q == m,allowed_indices) for m in rejected_indices]
		true_overlap = sort([overlap[ov][1] for ov in 1:length(overlap)])
		deleteat!(allowed_indices,true_overlap)

		row = [real(changed[i,j]) for j in allowed_indices]
		val, local_index = findmax(row)
		dats = findall(q->q == val,real.(changed[i,:]))
		index = dats[1]
		maxes[i] = changed[i,index]
		changed[i,:] = [changed[i,j] - maxes[i] for j in 1:num_parts]
		append!(rejected_indices,index)
	end
	final = changed
	reduced_logdet = sum(maxes) + log(det(exp.(final)))
	#reduced_regdet = exp(reduced_logdet)
	return reduced_logdet#,expected_reg_det,reduced_regdet,matrix,starting,final,maxes
end


function get_wavefunc_fromlog(config,n,p,qpart=[0,[0]])
	num_parts = length(config)
	log_matrix = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts
			data_here = get_log_elem_proj(config,j,i,n,p,qpart)
			log_matrix[i,j] += data_here
		end
	end
	result = get_log_det(log_matrix)
	for i in 1:num_parts
		result += -abs2(config[i])/4
	end
	return result
end




"fin"
