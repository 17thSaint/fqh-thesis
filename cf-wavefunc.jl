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
	println("Config: ",config)
	num_parts = length(config)
	#jiprime_elem = [1.0+0.0*im for i in part + 1:num_parts]
	for p in 1:num_parts
		if p == part
			continue
		end
		dist_btw = config[part]-config[p]
		println("Dist Btw $part and $p: ",dist_btw)
		ji *= dist_btw
		println("J $part: ",ji)
		#=
		for i in part + 1:num_parts
			if i == p
				jiprime_elem[i-part] *= 2*dist_btw
				continue
			end
			jiprime_elem[i-part] *= dist_btw^2
		end
		=#
	end
	#jiprime = sum(jiprime_elem)
	return ji#,jiprime
end

function get_Jiprime(config,part)
	jiprime = 0
	num_parts = length(config)
	for j in part + 1:num_parts
		prod_part = 2*abs(config[part] - config[j])
		for k in part + 1:num_parts
			if k == j
				continue
			end
			prod_part *= abs2(config[part] - config[k])
		end
		jiprime += prod_part
	end
	return jiprime
end

function get_elem_projection(config,num_parts,part,row)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element)
	Ji, Jiprime = get_Jis(config,num_parts,part)
	if bar > 0
		result = 2*(l*(config[part]^(l-1))*Ji + (config[part]^l)*Jiprime)
	else
		result = (config[part]^l)*Ji
	end
	return result
end

function get_wavefunc(config,num_parts)
	matrix_full = fill(0.0+0.0*im,(num_parts,num_parts))
	wavefunc = 1.0
	for i in 1:num_parts
		wavefunc *= exp(-abs2(config[i])/4)
		for j in 1:num_parts
			matrix_full[i,j] = get_elem_projection(config,num_parts,j,i)
		end
	end
	wavefunc *= det(matrix_full)
	return wavefunc
end

function get_logJi(config,part)
	logJi = 0.0+0.0*im
	num_parts = length(config)
	for i in part + 1:num_parts
		logJi += 2*log(config[part] - config[i])
	end
	return logJi
end

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
