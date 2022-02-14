using LinearAlgebra,PyPlot

function get_wf_elem(num_parts,element,n=2)
	if element[1] > num_parts/2
		m = 1
		pow = element[1] - num_parts/2 - 1
	else
		m = 0
		pow = element[1] - 1
	end
	return m,pow
end

function get_Jis(config,part)
	ji = 1.0
	num_parts = length(config)
	for p in 1:num_parts
		if p == part
			continue
		end
		dist_btw = config[part]-config[p]
		ji *= dist_btw
	end
	return ji
end

function get_Jiprime(config,part)
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
	return jiprime
end

function get_elem_projection(config,part,row)
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element)
	Ji, Jiprime = get_Jis(config,part),get_Jiprime(config,part)
	if bar > 0
		result = 2*(l*(config[part]^(l-1))*Ji + (config[part]^l)*Jiprime)
	else
		result = (config[part]^l)*Ji
	end
	return result
end

function get_wavefunc(config)
	num_parts = length(config)
	matrix_full = fill(0.0+0.0*im,(num_parts,num_parts))
	wavefunc = 1.0
	for i in 1:num_parts
		#wavefunc *= exp(-abs2(config[i])/4)
		for j in 1:num_parts
			matrix_full[i,j] = get_elem_projection(config,j,i)
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
		logji += log(dist_btw)
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
			next_part += log(dist_btw)
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

function get_det_fromlog(config)
	num_parts = length(config)
	log_matrix = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts	
			log_matrix[i,j] += get_log_elem_proj(config,j,i)
		end
	end
	reg_matrix = exp.(log_matrix)
	result = det(reg_matrix)
	return result
end

#=
function get_logJi(config,part)
	logJi = 0.0+0.0*im
	num_parts = length(config)
	for i in part + 1:num_parts
		logJi += 2*log(config[part] - config[i])
	end
	return logJi
end
=#

function get_log_elem_proj(config,num_parts,part,row)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element)
	logJi = get_logJi(config,num_parts,part)
	if bar > 0
		J_ratio = sum([2/(config[part]-config[j]) for j in part + 1:num_parts])
		result = l*log(config[part]) + logJi + log((l/config[part]) + J_ratio) + log(2)
	else
		result = l*log(config[part]) + logJi
	end
	return result
end

function get_log_wavefunc(config,num_parts)
	matrix_full = fill(0.0+0.0*im,(num_parts,num_parts))
	log_wavefunc = -0.25*sum(abs2.(config))
	for i in 1:num_parts
		for j in 1:num_parts
			matrix_full[i,j] = get_log_elem_proj(config,num_parts,j,i)
		end
	end
	sub_matrix = matrix_full[2:4,1:3]
	#log_wavefunc += det(sub_matrix)
	return log_wavefunc,matrix_full
end

#=
setup = Complex.([1*i for i in 1:4])
particles = length(setup)
reg_val = fill(0.0+0.0*im,(4,4))
log_val = fill(0.0+0.0*im,(4,4))
for i in 1:4
	for j in 1:4
		reg_val[i,j] = log(get_elem_projection(setup,particles,i,j))
		log_val[i,j] = get_log_elem_proj(setup,particles,i,j)
	end
end
#println(reg_val," ",log_val)
=#



"fin"
