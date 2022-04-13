using Statistics,PyPlot,LsqFit,LaTeXStrings
ARGS = false
include("read-CF-data.jl")
include("cf-wavefunc.jl")
include("energy-time-correlation.jl")

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

function get_local_count(x_edge,pos_data,iter_count,rm)
	local_count = 0
	for i in 1:iter_count
	    	xpos = real.(pos_data[1:end,i])
		ypos = -imag.(pos_data[1:end,i])
		for j in 1:length(xpos)
		    if xpos[j] <= x_edge && xpos[j] >= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
		        local_count += 1
		    end
		end
	end
	return local_count
end

function get_density(particles,n,p,hdf5_data,data_count,samp_freq,localized=false)
	rm = sqrt(2*particles*(2*p*n+1)/n)
	denom = 2*p*n+1
	qpart_data,pos_data_og,wavefunc_data = hdf5_data
	if localized
		xs = [real(qpart_data[2][1]) - 0.3*rm + 0.6*rm*i/data_count for i in 0:data_count-1]
	else
		xs = [i*1.3*rm/data_count + 0.5*1.3*rm/data_count for i in 0:data_count-1]
	end
	println(length(wavefunc_data),", ",samp_freq,", ",qpart_data[1])
	number_slices = Int(floor(length(wavefunc_data)/samp_freq))
	pos_data = fill(0.0+im*0.0,(particles,number_slices))
	pos_data[:,1] = pos_data_og[:,1]
	if samp_freq != 1
		for i in 1:number_slices
			pos_data[:,i] = pos_data_og[:,i*samp_freq]  
		end
	else
		pos_data = pos_data_og
	end
	qp_count = qpart_data[1]
	change_counts = [0.0 for i in 1:data_count]
	total_possible_counts = particles*length(wavefunc_data)
	#errors = [0.0 for i in 1:data_count]
	if localized
		zeroth_x = real(qpart_data[2][1]) - rm*(0.3 + 0.6/data_count)
		previous_count = get_local_count(zeroth_x,pos_data,number_slices,rm)
	else
		previous_count = 0
	end
	for k in 1:data_count
	    if k%(data_count*0.05) == 0
		    println("Running QP $qp_count $n/$denom: ",100*k/data_count,"%")
	    end	
	    edge = xs[k]
	    local_count = get_local_count(edge,pos_data,number_slices,rm)
	    
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

# need to fix curve fitting on rads 6,7
xs_plot = [[],[],[],[]]
ys_plot = [[],[],[],[]]

particles = 10
log_form = true
flux_type = "RFA"
low_fluxtype = "rfa"
#mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = np_vals[which]
fill_denom = 2*n*p + 1
rm = sqrt(2*particles*fill_denom/n)
#
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
rad_choice = 5
#

#
#for rad_choice in 5:5
qp2_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,2,log_form)
qp1_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,1,log_form)
figure()
array_data_1q = Iterators.flatten([qp1_data[2][i,:] for i in 1:particles])
hist2D(real.(array_data_1q),-imag.(array_data_1q),bins=100)
title("1Q-$rad_choice")
figure()
array_data_2q = Iterators.flatten([qp2_data[2][i,:] for i in 1:particles])
hist2D(real.(array_data_2q),-imag.(array_data_2q),bins=100)
title("2Q-$rad_choice")
#end
#
#=
rad_vals = [5]
exp_pvals = [0.0 for i in 1:length(rad_vals)]
errors = [0.0 for i in 1:length(rad_vals)]
densities_1q = []
densities_2q = []
fitcurves_1q = []
fitcurves_2q = []
for i in 1:length(rad_vals)
#i = 5
rad_choice = rad_vals[i]
x_rad = x_rads[rad_choice]
qp2_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,2,log_form)
qp1_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,1,log_form)

xs_count = 200
autocorr_length_qp1 = 1
autocorr_length_qp2 = 1
dens_data_1q = get_density(particles,n,p,qp1_data,xs_count,autocorr_length_qp1,false)
dens_data_2q = get_density(particles,n,p,qp2_data,xs_count,autocorr_length_qp2,false)
append!(densities_1q,[dens_data_1q])
append!(densities_2q,[dens_data_2q])
qpart_location = qp1_data[1][2][1]
#=
rad_range_1q = findall(x->isapprox(x,real(qpart_location)/rm,atol=0.1),dens_data_1q[1]./rm)
rad_range_2q = findall(x->isapprox(x,real(qpart_location)/rm+0.05,atol=0.1),dens_data_2q[1]./rm)
vert_shift_1q = minimum(dens_data_1q[2][rad_range_1q[1]:rad_range_1q[end]])
vert_shift_2q = minimum(dens_data_2q[2][rad_range_2q[1]:rad_range_2q[end]])
scaling_guess = maximum(dens_data_1q[2])

model_1q(qrad,p) = p[1].*(qrad .- p[2]).^2 .+ vert_shift_1q
model_2q(qrad,p) = p[1].*(qrad .- p[2]).^2 .+ vert_shift_2q
p0 = [scaling_guess,x_rad]
fit_1q = curve_fit(model_1q,dens_data_1q[1][rad_range_1q[1]:rad_range_1q[end]],dens_data_1q[2][rad_range_1q[1]:rad_range_1q[end]],p0)
fit_2q = curve_fit(model_2q,dens_data_2q[1][rad_range_2q[1]:rad_range_2q[end]],dens_data_2q[2][rad_range_2q[1]:rad_range_2q[end]],p0)

append!(fitcurves_1q,[model_1q(dens_data_1q[1][rad_range_1q[1]:rad_range_1q[end]],fit_1q.param)])
append!(fitcurves_2q,[model_2q(dens_data_2q[1][rad_range_2q[1]:rad_range_2q[end]],fit_2q.param)])
=#
if true
#figure()
plot(densities_1q[end][1]./rm,densities_1q[end][2],label="1Q")
plot(densities_2q[end][1]./rm,densities_2q[end][2],label="2Q")
#plot(densities_1q[end][1][rad_range_1q[1]:rad_range_1q[end]]./rm,fitcurves_1q[end],"-c")
#plot(densities_2q[end][1][rad_range_2q[1]:rad_range_2q[end]]./rm,fitcurves_2q[end],"-m")
scatter([real(qpart_location)]./rm,[densities_1q[end][2][1]])#,label="QP")
title(latexstring("Density for \$ \\nu=$n/$fill_denom\$"))
legend()
end
#=
exp_pval = (fit_2q.param[2]^2 - fit_1q.param[2]^2)/4
error = sqrt(2*(stderror(fit_1q)[2]*fit_1q.param[2])^2 + 2*(stderror(fit_2q)[2]*fit_2q.param[2])^2)
#println(exp_pval," +/- ",error)
exp_pvals[i] = exp_pval
errors[i] = error
=#

end

=#

#=
figure()
plot(x_rads[rad_vals[1]:rad_vals[end]]./rm,errors,"-p")
title("Errors vs QP Radii")

figure()
errorbar(x_rads[rad_vals[1]:rad_vals[end]]./rm,exp_pvals,yerr=[errors,errors],fmt="-o")
title("P-Value vs QP Radii")
=#
#=
plot(dens_data_1q[1]./rm,dens_data_1q[2],label="1-QP")#"$n/$fill_denom")
plot(dens_data_2q[1]./rm,dens_data_2q[2],label="2-QP")#"$n/$fill_denom")
scatter([real(qpart_location)]./rm,[dens_data_1q[2][1]],label="QP")
#
plot(dens_data_1q[1][rad_range_1q[1]:rad_range_1q[end]]./rm,model_1q(dens_data_1q[1][rad_range_1q[1]:rad_range_1q[end]],fit_1q.param),label="Fit 1Q")
plot(dens_data_2q[1][rad_range_2q[1]:rad_range_2q[end]]./rm,model_2q(dens_data_2q[1][rad_range_2q[1]:rad_range_2q[end]],fit_2q.param),label="Fit 2Q")
#
legend()
=#


#=
xs = [-0.5*rm + i*(2*0.5*rm)/data_count for i in 1:data_count]
for i in 1:length(xs)
	local_x = xs[i]
	println(i/length(xs))
	for j in 1:length(xs)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		
		
	end
end
=#

#scatter3D(xs_plot,ys_plot,reg_prob./maximum(reg_prob),label="Reg")
#scatter3D(xs_plot,ys_plot,exp_prob./maximum(exp_prob),label="Exp")
#sleep(2.0)
#scatter3D(xs_plot,ys_plot,part1+logjis,label="One")

#legend()
#xlabel("X")
#ylabel("Y")
#
#

"fin"
