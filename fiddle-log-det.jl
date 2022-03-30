using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")

function get_logdet_diff(log_matrix,exp_shift)
	parts = length(exp_shift)
	multip_log_matrix = fill(0.0+im*0.0,(parts,parts))
	#exp_shift = [rand(Float64)+im*rand(Float64) for j in 1:parts]
	for j in 1:parts
	multip_log_matrix[:,j] = log_matrix[:,j] .+ exp_shift[j]
	end
	
	logdet_multip = get_log_det(multip_log_matrix)
	#println(multip_log_matrix,", ",log_matrix)
	logdet_reg = get_log_det(log_matrix)
	expected_shift = sum(exp_shift)
	logdet_expected = logdet_reg + expected_shift
	
	reg_diff = exp(logdet_multip) - exp(logdet_expected)
	log_diff = logdet_multip - logdet_expected
	#println(round(imag(diff)/pi,digits=0),", ",imag(diff)/pi)
	#if !isapprox(round(imag(diff)/pi,digits=0),imag(diff)/pi,atol=10^(-1))
	#	println(imag(diff)/pi)
	#end
	return multip_log_matrix
end

function get_regdet_diff(log_matrix,exp_shift)
	parts = length(exp_shift)
	multip_reg_matrix = fill(0.0+im*0.0,(parts,parts))
	reg_matrix = exp.(log_matrix)
	reg_shift = exp.(exp_shift)
	for j in 1:parts
	multip_reg_matrix[:,j] = reg_matrix[:,j] .* reg_shift[j]
	end
	det_multip = det(multip_reg_matrix)
	
	det_reg = det(reg_matrix)
	expected_shift = prod(reg_shift)
	det_expected = det_reg * expected_shift
	log_det_multip = log(det_multip)
	log_det_expected = log(det_reg) + log(expected_shift)#log(det_expected)
	
	log_diff = log_det_multip - log_det_expected
	reg_diff = det_multip - det_expected
	
	return log_diff
end


function get_rand_matrix(num_parts)
	mat = fill(0.0+im*0.0,(num_parts,num_parts))
	for i in 1:num_parts
		for j in 1:num_parts
			mat[i,j] = rand(Float64)*rand((-1,1)) + im*rand(Float64)*rand((-1,1))
		end
	end
	return mat
end

function get_rand_shift(num_parts)
	shift = [rand(Float64)*rand((-1,1)) + im*rand(Float64)*rand((-1,1)) for i in 1:num_parts]
	return shift
end

function check_matrix_equality(matrix_one,matrix_two)
	same_real,same_imag = true,true
	counts = size(matrix_one)[1]
	real_perc_diff = abs.((real.(matrix_one) - real.(matrix_two))./real.(matrix_one))
	if mean(real_perc_diff) > 0.01
		same_real = false
	end
	
	imag_diff_mod = (imag.(matrix_one) - imag.(matrix_two))./(2*pi)
	rounded_diffs = round.(imag_diff_mod,digits=0)
	remainder = abs.(imag_diff_mod - rounded_diffs)
	if mean(remainder) > 10^(-1)
		same_imag = false
	end
	return same_real,same_imag
end

#
data_count = 50

cf_ver = [fill(0.0,(data_count,data_count)) for i in 1:4]
laugh_ver = [fill(0.0,(data_count,data_count)) for i in 1:4]

xs_plot = [[],[],[],[]]
ys_plot = [[],[],[],[]]

parts_vals = [6,8,10,12]
for k in 1:4


particles = parts_vals[k]
mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = 1,1
fill_denom = 2*n*p + 1
rm = sqrt(2*particles*fill_denom/n)
#=
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
rad_choice = 5
x_rad = x_rads[rad_choice]
=#
qpart = [0,[0.0+im*0.0]]
starting_config = start_rand_config(particles,n,p)
xs = [-0.5*rm + i*(2*0.5*rm)/data_count for i in 1:data_count]
for i in 1:length(xs)
	local_x = xs[i]
	println(i/length(xs))
	for j in 1:length(xs)
		local_y = xs[j]
		radius = round(abs(local_x - im*local_y),digits=2)
		starting_config[1] = local_x - im*local_y
		
		vals_CF = get_wavefunc_fromlog(starting_config,n,p,qpart)
		cf_ver[k][i,j] = 2*real(vals_CF[2])
		laugh_ver[k][i,j] = -prob_wavefunc_laughlin(starting_config, 2*p*n+1)
		
		#=
		full_matrix = vals_CF
		diag_rez,diag_vals,diag_redmat = get_diag_log_det(full_matrix)
		maxes_rez,max_vals,maxes_redmat = get_log_det(full_matrix)
		actual_rez = log(Complex(det(exp.(full_matrix))))
		
		fromdiag_ogmat = fill(0.0+im*0.0,(particles,particles))
		frommaxes_ogmat = fill(0.0+im*0.0,(particles,particles))
		for l in 1:particles
			fromdiag_ogmat[:,l] = diag_redmat[:,l] .+ diag_vals[l]
			frommaxes_ogmat[:,l] = maxes_redmat[:,l] .+ maxes_redmat[l]
		end
		eq = check_matrix_equality(fromdiag_ogmat,frommaxes_ogmat)
		
		if eq[1]
			println("Same Real $radius")
		end
		if eq[2]
			println("Same Imag $radius")
		end
		
		diff_rezz_logdet[k][i,j] = real(actual_rez-diag_rez)
		diff_maxdiag_logdet[k][i,j] = real(maxes_rez-diag_rez)
		#println("Radius $radius: ",diff_rezz_logdet[k][i,j])
		#sum(max_vals),", ",log(Complex(det(exp.(maxes_redmat)))))
		=#
	end
end
end

#
for i in 1:4
if true
figure()
imshow(cf_ver[i])
title("CF $i")
colorbar()
end
if true
figure()
imshow(laugh_ver[i])
title("Laughlin $i")
colorbar()
end
end

