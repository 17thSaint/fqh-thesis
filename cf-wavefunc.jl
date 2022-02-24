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

function get_qpart_wf_exp(config,qpart,which_qpart,chosen_qpart,n)
	lstar2 = 2*1*n + 1
	if chosen_qpart[1] == "particle"
		position = config[chosen_qpart[2]]
	else
		position = qpart[2][chosen_qpart[2]]
	end
	result = exp(conj(qpart[2][which_qpart])*position/(2*lstar2) - abs2(qpart[2][which_qpart])/(4*lstar2))
	return result
end

function get_Jis(config,part,qpart=[0,[0]],qhole=[0,0])
	ji = 1.0
	num_parts = length(config)
	for p in 1:num_parts
		if p == part
			continue
		end
		dist_btw = config[part]-config[p]
		ji *= dist_btw
	end
	#
	for q in 1:qpart[1]
		dist_btw = config[part]-qpart[2][q]
		ji *= dist_btw^2
	end
	#
	if real(qhole[1]) > 0.1
		ji *= config[part] - qhole[2]
	end
	return ji
end

function get_Jiprime(config,part,qpart=[0,[0]],qhole=[0,0])
	jiprime = 0
	num_parts = length(config)
	for j in 1:num_parts + qpart[1]
		if j == part
			continue
		end
		if j > num_parts
			pos_selected = qpart[2][part-num_parts]
		else
			pos_selected = config[part]
		end
		jiprime_local = 1
		for k in 1:num_parts + qpart[1]
			if k == part
				continue
			end
			if k == j
				continue
			end
			if k > num_parts
				dist_btw = pos_selected - qpart[2][k-num_parts]
			else
				dist_btw = pos_selected - config[k]
			end
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

function get_elem_projection(config,part,row,n,qpart=[0,[0]],qhole=[0,0])
	num_parts = length(config)
	matrix_element = [row,part]
	qpart_shift = qpart[1]
	if row >= qpart_shift + 1
		bar, l = get_wf_elem(num_parts,matrix_element,n)
		if part <= num_parts # bottom left all particles
			Ji, Jiprime = get_Jis(config,part,qhole),get_Jiprime(config,part,qhole)
			section = "BL"
			if bar > 0
				result = 2*(l*(config[part]^(l-1))*Ji + (config[part]^l)*Jiprime)
			else
				result = (config[part]^l)*Ji
			end
		else # bottom right polynomial is quasiparticle
			section = "BR"
			if bar > 0
				result = 0.0#qpart[2][part-num_parts]*Jiprime
			else
				result = (qpart[2][part-num_parts]^l)#*Ji
			end
		end
	elseif row < qpart_shift + 1
		lstar2 = 2*1*n + 1
		coeff = [(1/lstar2) - 1,1/lstar2^2 - 2\lstar2 + 1]
		if part <= num_parts # top left quasiparticle coherent with particles
			section = "TL"
			Ji, Jiprime = get_Jis(config,part,qhole),get_Jiprime(config,part,qhole)
			exp_part = get_qpart_wf_exp(config,qpart,row,["particle",part],n)
			result = coeff[n]*(conj(qpart[2][row])^n)*exp_part*Ji + 2*Jiprime*exp_part
		else # top right quasiparticle coherent with quasiparticles
			section = "TR"
			exp_part = get_qpart_wf_exp(config,qpart,row,["quasi",part-num_parts],n)
			result = coeff[n]*(conj(qpart[2][row])^n)*exp_part
		end
	end
	return result#,section
end

function get_wavefunc(config,n,qpart=[0,[0]],qhole=[0,0])
	num_parts = length(config)
	matrix_full = fill(0.0+0.0*im,(num_parts+qpart[1],num_parts+qpart[1]))
	#matrix_display = fill("",(num_parts+qpart[1],num_parts+qpart[1]))
	wavefunc = 1.0
	for i in 1:num_parts + qpart[1]
		if i <= num_parts
			wavefunc *= exp(-abs2(config[i])/4)
		end
		for j in 1:num_parts + qpart[1]
			dats = get_elem_projection(config,j,i,n,qpart,qhole)
			matrix_full[i,j] = dats
			#matrix_display[i,j] = dats[2]
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

function get_log_elem_proj(config,part,row,n)
	num_parts = length(config)
	matrix_element = [row,part]
	bar, l = get_wf_elem(num_parts,matrix_element,n)
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

function get_wavefunc_fromlog(config,n)
	num_parts = length(config)
	log_matrix = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts	
			log_matrix[i,j] += get_log_elem_proj(config,j,i,n)
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
