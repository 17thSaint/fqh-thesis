using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")
include("mc-cf.jl")

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
	return log_diff
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

function get_density(particles,n,p,hdf5_data,data_count)
	rm = sqrt(2*particles*(2*p*n+1)/n)
	xs = [i*1.3*rm/data_count + 0.5*1.3*rm/data_count for i in 0:data_count-1]
	qpart_data,pos_data,wavefunc_data = hdf5_data
	change_counts = [0.0 for i in 1:data_count]
	total_possible_counts = particles*length(wavefunc_data)
	#errors = [0.0 for i in 1:data_count]
	previous_count = 0
	for k in 1:data_count
	    edge = xs[k]
	    #
	    local_count = 0
	    for i in 1:length(wavefunc_data)
	    	xpos = real.(pos_data[1:end,i])
		ypos = -imag.(pos_data[1:end,i])
		for j in 1:length(xpos)
		    if xpos[j] <= edge && xpos[j] >= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
		        local_count += 1
		    end
		end
	    end
	    
	    change_counts[k] = (local_count - previous_count)/total_possible_counts
	    previous_count = local_count
	end
	return xs,change_counts
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
	same = true
	counts = size(matrix_one)[1]
	real_perc_diff = abs.((real.(matrix_one) - real.(matrix_two))./real.(matrix_one))
	if mean(real_perc_diff) > 0.01
		same = false
	end
	
	imag_diff_mod = (imag.(matrix_one) - imag.(matrix_two))./(2*pi)
	rounded_diffs = round.(imag_diff_mod,digits=0)
	remainder = abs.(imag_diff_mod - rounded_diffs)
	if mean(remainder) > 10^(-1)
		same = false
	end
	return same
end

#
data_count = 20
#exp_prob = []
#reg_prob = []
logjast = fill(0.0+im*0.0,(data_count,data_count))
logj1 = fill(0.0+im*0.0,(data_count,data_count))
logslater = fill(0.0+im*0.0,(data_count,data_count))
#part1 = []
exp_prob_CF_sep = fill(0.0,(data_count,data_count))
exp_prob_CF_reg = fill(0.0,(data_count,data_count))
exp_prob_Laugh = fill(0.0,(data_count,data_count))

exact_diffs = fill(0.0+im*0.0,(data_count,data_count))
#randmat_diffs = [fill(0.0+im*0.0,(data_count,data_count)) for i in 1:length(coef_vals)]
randshift_diffs_log = [fill(0.0+im*0.0,(data_count,data_count)) for i in 1:4]
randshift_diffs_reg = [fill(0.0+im*0.0,(data_count,data_count)) for i in 1:4]
randshift_diffs_btw = [fill(0.0+im*0.0,(data_count,data_count)) for i in 1:4]
#givematrix = [fill(0.0+im*0.0,(particles,particles))]
trial_vals = []
xs_plot = [[],[],[],[]]
ys_plot = [[],[],[],[]]

parts_vals = [6,8,10,12]
for k in 1:4


particles = parts_vals[k]#12
mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = 1,0#np_vals[which]
fill_denom = 2*n*p + 1
rm = sqrt(2*particles*fill_denom/n)
#=
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
rad_choice = 5
x_rad = x_rads[rad_choice]
=#
qpart = [0,[3.0+im*3.0]]
#qp2_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,2,true)
#qp1_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,1,true)
#

starting_config = start_rand_config(particles,n,p)#1 .*qp1_data[2][:,3800]
#parts_xs = real.(starting_config[2:end])
#parts_ys = -imag.(starting_config[2:end])
#scatter3D(parts_xs,parts_ys,[1.0 for i in 1:particles-1],c="r")
#xs = [-0.01*rm + i*(2*0.01*rm)/data_count for i in 1:data_count]
xs = [-0.01*rm + i*(2*0.01*rm)/data_count for i in 1:data_count]
for i in 1:length(xs)
	local_x = xs[i]
	#println(i/length(xs))
	for j in 1:length(xs)
		local_y = xs[j]
		#append!(xs_plot,[local_x])
		#append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		radius = abs(starting_config[1])
		#=
		logjast[i,j] = get_logJastrowfull(starting_config)
		logj1[i,j] = get_logJi(starting_config,1)
		logslater[i,j] = log(starting_config[1])
		=#
		vals_CF = get_wavefunc_fromlog(starting_config,n,p,qpart)
		#append!(trial_vals,[vals_CF,starting_config])
		
		sep_matrix = vals_CF[4]
		exp_shift = [p*get_logJi(starting_config,k) for k in 1:particles]
		
		#exact_diff_local = get_logdet_diff(sep_matrix,exp_shift)
		#exact_diffs[i,j] = exact_diff_local
		#=
		multip = coef_vals[k]
		randmat = multip.*get_rand_matrix(particles)
		randmat_diff_local = get_logdet_diff(randmat,exp_shift)
		randmat_diffs[k][i,j] = randmat_diff_local
		=#
		
		randshift = get_rand_shift(particles)
		randshift_logdiff_local = get_logdet_diff(sep_matrix,randshift)
		#randshift_regdiff_local = get_regdet_diff(sep_matrix,randshift)
		#fromlog = randshift_logdiff_local[1]
		#fromreg = log.(randshift_regdiff_local[1])
		#matrix_equality = check_matrix_equality(fromlog,fromreg)
		#if !matrix_equality
		#	println("Not Equal")
		#end
		
		#fromlog_det = get_diag_log_det(fromlog)
		#fromlog_det_reg = log(det(exp.(fromlog)))
		#fromreg_det = get_log_det(fromreg)
		#fromreg_det_reg = log(det(exp.(fromreg)))
		#percent_diff_log = abs((real(fromlog_det) - real(fromlog_det_reg))/real(fromlog_det))
		#percent_diff_reg = abs((real(fromreg_det) - real(fromreg_det_reg))/real(fromreg_det))
		#=
		if percent_diff > 0.01
			rounded_perdiff = round(percent_diff,digits=3)
			println("Different $particles $rounded_perdiff: ")#,fromlog,", ",fromreg)
		end
		=#
		
		
		#randshift_diffs_reg[k][i,j] = percent_diff_reg
		randshift_diffs_log[k][i,j] = randshift_logdiff_local
		#=
		reg_diff = real(randshift_logdiff_local[2] - randshift_regdiff_local[2])
		log_diff = real(randshift_logdiff_local[1] - randshift_regdiff_local[1])
		if real(reg_diff) > 10^(-5)
			println("Reg $local_x $local_y: ",randshift_logdiff_local[2],", ", randshift_regdiff_local[2])
		end
		if real(log_diff) > 10^(-5)
			println("Log $local_x $local_y: ",randshift_logdiff_local[1],", ", randshift_regdiff_local[1])
		end
		if isnan(log_diff)
			println("NaN Log $local_x $local_y: ",randshift_logdiff_local[1],", ", randshift_regdiff_local[1])
		end
		if isnan(reg_diff)
			println("NaN Reg $local_x $local_y: ",randshift_logdiff_local[2],", ", randshift_regdiff_local[2])
		end
		=#
		
		
		#=
		exp_prob_CF_sep_local = 2*real(vals_CF[2])
		exp_prob_CF_reg_local = 2*real(vals_CF[1])
		exp_prob_CF_sep[i,j] = exp_prob_CF_sep_local
		exp_prob_CF_reg[i,j] = exp_prob_CF_reg_local
		
		vals_Laugh = prob_wavefunc_laughlin(starting_config,2*p+1)
		exp_prob_Laugh_local = -vals_Laugh
		exp_prob_Laugh[i,j] = exp_prob_Laugh_local
		=#
	end
end
end

#
if false
figure()
imshow(real.(exact_diffs))#./maximum(exp_prob_CF))
title("Exact Diffs")
colorbar()
end
if false
figure()
#imshow(real.(randmat_diffs))
plot(coef_vals,errors,"-p")
title("Rand Matrix Diff Errors")
#colorbar()
end
for i in 1:4
part_count = parts_vals[i]
if true
figure()
imshow(real.(randshift_diffs_log[i]))
title("Rand Shift Diffs Log")
colorbar()
end
end
if false
figure()
imshow(abs.(exp_prob_CF_sep-exp_prob_CF_reg))#./maximum(exp_prob_CF))
title("Sep CF Wavefunc")
colorbar()
end
if false
figure()
imshow(exp_prob_CF_reg)#./maximum(exp_prob_CF))
title("Reg CF Wavefunc")
colorbar()
end
if false
figure()
imshow(exp_prob_Laugh)#./maximum(exp_prob_Laugh))
title("Laughlin Wavefunc")
colorbar()
end
if false
figure()
imshow(real.(logjast))
title("Jastrow")
colorbar()
end
if false
figure()
imshow(real.(logj1))
title("J_one")
colorbar()
end
if false
figure()
imshow(real.(logslater))
title("Slater")
colorbar()
end
if false
figure()
imshow(exp_prob_CF./real.(logjast))
title("Quotient")
colorbar()
end



#scatter3D(xs_plot,ys_plot,reg_prob./maximum(reg_prob),label="Reg")
#scatter3D(xs_plot,ys_plot,exp_prob./maximum(exp_prob),label="Exp")
#sleep(2.0)
#scatter3D(xs_plot,ys_plot,part1+logjis,label="One")

#legend()
#xlabel("X")
#ylabel("Y")
#






#array_data = Iterators.flatten([qp1_data[2][i,:] for i in 1:1])#Iterators.flatten([only_no_overlaps_data[2][i,:] for i in 2:particles])
#hist2D(real.(array_data),-imag.(array_data),bins=100)
#
#dens_data_1q = get_density(particles,n,p,qp1_data,200)
#dens_data_2q = get_density(particles,n,p,qp2_data,200)
#qpart_location = qp1_data[1][2][1]
#
#plot(dens_data_1q[1]./rm,dens_data_1q[2],label="1-QP")#"$n/$fill_denom")
#plot(dens_data_2q[1]./rm,dens_data_2q[2],label="2-QP")#"$n/$fill_denom")
#scatter([real(qpart_location)]./rm,[dens_data_1q[2][1]],label="QP")
#legend()
#end
#

"fin"
