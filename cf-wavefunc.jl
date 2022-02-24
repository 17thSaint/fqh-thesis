using LinearAlgebra,PyPlot

function get_wf_elem(num_parts,element,qpart=[0,[0]],n=2)
	qpart_shift = qpart[1]
	if element[1] > num_parts/2 + qpart_shift
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

function get_qpart_wf_exp(config,qpart,which_qpart,particle)
	lstar2 = 2*1*2 + 1
	result = exp(conj(qpart[2][which_qpart])*config[particle]/(2*lstar2) - abs2(qpart[2][which_qpart])/(4*lstar2))
	return result
end

function get_Jis(config,part,qhole=[0,0])
	ji = 1.0
	num_parts = length(config)
	for p in 1:num_parts
		if p == part
			continue
		end
		dist_btw = config[part]-config[p]
		ji *= dist_btw
	end
	if real(qhole[1]) > 0.1
		ji *= config[part] - qhole[2]
	end
	return ji
end

function get_Jiprime(config,part,qhole=[0,0])
	jiprime = 0
	num_parts = length(config)
	for j in 1:num_parts
		if j == part
			continue
		end
		jiprime_local = 1
		for k in 1:num_parts
			if k == part
				continue
			end
			if k == j
				continue
			end
			dist_btw = config[part] - config[k]
			jiprime_local *= dist_btw
		end
		jiprime += 2*jiprime_local
	end
	if real(qhole[1]) > 0.1
		jiprime *= config[part] - qhole[2]
		jiprime += 2*get_Jis(config,part,qhole)
	end
	return jiprime
end

function get_elem_projection(config,part,row,qpart=[0,[0]],qhole=[0,0])
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element)
	Ji, Jiprime = get_Jis(config,part,qhole),get_Jiprime(config,part,qhole)
	if qpart[1] > 0 && row <= qpart[1]
		exp_part = get_qpart_wf_exp(config,qpart,row,part)
		result = (16/25)*conj(qpart[2][row])*exp_part
	elseif bar > 0
		result = 2*(l*(config[part]^(l-1))*Ji + (config[part]^l)*Jiprime)
	else
		result = (config[part]^l)*Ji
	end
	return result
end

function get_wavefunc(config,qpart=[0,[0]],qhole=[0,0])
	num_parts = length(config)
	matrix_full = fill(0.0+0.0*im,(num_parts,num_parts))
	wavefunc = 1.0
	for i in 1:num_parts
		wavefunc *= exp(-abs2(config[i])/4)
		for j in 1:num_parts
			matrix_full[i,j] = get_elem_projection(config,j,i,qpart,qhole)
		end
	end
	wavefunc *= det(matrix_full)
	return wavefunc
end

function get_logJi(config,part)
	logji = 0.0
	num_parts = length(config)
	for p in 1:num_parts
		if p == part
			continue
		end
		dist_btw = config[part]-config[p]
		logji += log(Complex(dist_btw))
	end
	return logji
end

function get_logJiprime(config,part)
	logjiprime = 0.0+im*0.0
	num_parts = length(config)
	start = 0
	for i in 1:num_parts
		if i == part
			continue
		end
		next_part = 0.0+im*0.0
		for j in 1:num_parts
			if j == i
				continue
			end
			if j == part
				continue
			end
			dist_btw = config[part] - config[j]
			next_part += log(Complex(dist_btw))
		end
		if start == 0
			logjiprime += next_part
		else
			logjiprime += log(1 + exp(next_part - logjiprime))
		end
		start += 1
	end
	logjiprime += log(2)
	return logjiprime
end	

function get_log_elem_proj(config,part,row)
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element)
	logJi, logJiprime = get_logJi(config,part), get_logJiprime(config,part)
	if bar > 0
		if l == 0
			result = log(2) + logJiprime
		else
			a = log(2) + log(l) + (l-1)*log(config[part]) + logJi
			b = log(2) + l*log(config[part]) + logJiprime
			result = a + log(1 + exp(b-a))
		end
	else
		result = l*log(config[part]) + logJi
	end
	return result
end

function get_wavefunc_fromlog(config)
	num_parts = length(config)
	log_matrix = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts	
			log_matrix[i,j] += get_log_elem_proj(config,j,i)
		end
	end
	reg_matrix = exp.(log_matrix)
	result = det(reg_matrix)
	for i in 1:num_parts
		result *= exp(-abs2(config[i])/4)
	end
	return result
end




"fin"
