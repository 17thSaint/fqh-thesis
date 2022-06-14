using Statistics,PyPlot,LsqFit,LaTeXStrings

include("cf-wavefunc.jl")
ARGS = "F"
include("read-CF-data.jl")

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
			if x_edge > 0.0
				if xpos[j] <= x_edge && xpos[j] >= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
					local_count += 1
			    	end
			#=
			elseif x_edge < 0.0
				if xpos[j] >= x_edge && xpos[j] <= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
					local_count += 1
			    	end
			end
			=#
			end
		end
	end
	return local_count
end

function get_density(flux_type,particles,n,p,hdf5_data,data_count,samp_freq,localized=false)
	if flux_type == "RFA"
		rm = sqrt(2*particles*(2*p*n-1)/n)
		denom = 2*p*n-1
	elseif flux_type == "CF"
		rm = sqrt(2*particles*(2*p*n+1)/n)
		denom = 2*p*n+1
	end
	qpart_data,pos_data_og,wavefunc_data = hdf5_data
	if localized
		xs = [real(qpart_data[2][1]) - 0.3*rm + 0.6*rm*i/data_count for i in 0:data_count-1]
	else
		xs = [i*1.3*rm/data_count + 0.5*1.3*rm/data_count for i in 0:data_count-1]
		#xs = [-0.5*rm + i*1.8*rm/data_count for i in 0:data_count-1]
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
		#zeroth_x = xs[1] - 1.8*rm/data_count
		#previous_count = get_local_count(zeroth_x,pos_data,number_slices,rm)
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

function circle_plot(rad,center,col="r")
	#=
	point_count = 10
	xs = [(-rad + 2*rad*i/point_count)+center[1] for i in 0:point_count]
	println(xs)
	ys = sqrt.(rad^2 .- xs.^2)
	ys_top = ys .+ center[2]
	ys_bot = -ys .+ center[2]
	plot(xs,ys_top,"-r")
	plot(xs,ys_bot,"-r")
	=#
	c1 = plt.Circle((center[1],center[2]),rad,fill=false,color=col)
	plt.gcf().gca().add_artist(c1)
	return
end

particles = 20
log_form = true
for f in [["CF",1,0.7]]
flux_type = f[1]
low_fluxtype = lowercase(flux_type)
#mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = f[2]
n,p = np_vals[which]
if flux_type == "RFA"
fill_denom = 2*n*p - 1
elseif flux_type == "CF"
fill_denom = 2*n*p + 1
end
rm = sqrt(2*particles*fill_denom/n)
#
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
for rad_choice in [4]
#rad_choice = 4
#

#
#for rad_choice in 5:5
#qp0_data = read_comb_CF_hdf5("$low_fluxtype-data/mirror-qpart",flux_type,particles,n,p,rad_choice,0,log_form) /mirror-qpart
qp2_data = read_comb_CF_hdf5("$low_fluxtype-data/mirror-qpart",flux_type,particles,n,p,rad_choice,2,log_form)
qp1_data = read_comb_CF_hdf5("$low_fluxtype-data/mirror-qpart",flux_type,particles,n,p,rad_choice,1,log_form)
#qp0_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,0,log_form)
comp_1q = collect(Iterators.flatten([qp1_data[2][i,:] for i in 1:particles]))
comp_2q = collect(Iterators.flatten([qp2_data[2][i,:] for i in 1:particles]))
q_rad_1 = round(real(qp2_data[1][2][1]),digits=3)
q_rad_2 = round(real(qp2_data[1][2][2]),digits=3)

#comp_0q = collect(Iterators.flatten([qp0_data[2][i,:] for i in 1:particles]))
#=
top = 0.5*rm
selected_dats = []
for i in 1:length(qp0_data[3])
	for j in 1:particles
		rad = abs(qp0_data[2][j,i])
		if rad <= top
			append!(selected_dats,[qp0_data[2][j,i]])
		end
	end
end
=#
#comp_0q = collect(Iterators.flatten([qp0_data[2][i,:] for i in 1:particles]))
#hist(abs.(comp_1q),bins=100)
#hist(abs.(comp_2q),bins=100)
#figure()
#hist2D(real.(selected_dats),-imag.(selected_dats),bins=100)
#title("No QP")
#end
which_map = "nipy_spectral"
if true
#figure()
#hist2D(real.(comp_0q),-imag.(comp_0q),bins=50)
#title("$flux_type 0Q")
#colorbar()
if true
figure()
hist_dats_1q = hist2D(real.(comp_1q),-imag.(comp_1q),bins=40,cmap=which_map)
#zs = collect(Iterators.flatten([hist_dats_1q[1][i,:] for i in 1:40]))
#xs = collect(Iterators.flatten([[j for i in 1:40] for j in 1:40]))
#ys = collect(Iterators.flatten([[40-i for i in 1:40] for j in 1:40]))
#x = 1:1:40
#y = 1:1:40
#figure()
#s1 = surface(x,y,(x,y)->hist_dats_1q[1][x,y])
#display(s1)
#
if flux_type == "RFA"
plot([real(q_rad_1)+f[3]],[0.0],"pr")
circle_plot(2.0,[real(q_rad_1)+f[3],0.0])
elseif flux_type == "CF"
plot([real(q_rad_1)],[0.0],"pr")
circle_plot(2.0,[real(q_rad_1),0.0])
end

title(latexstring("$flux_type Single QP: \$ \\nu=1/3 \$ N=10"))
colorbar()
#
end

if true
figure()
hist_dats_2q = hist2D(real.(comp_2q),-imag.(comp_2q),bins=60,cmap=which_map)

if flux_type == "RFA"
circle_plot(2.0,[real(q_rad_1)+f[3],0.0])
plot([real(q_rad_1)+f[3]],[0.0],"pr")

plot([real(q_rad_1)],[0.0],"pk")
plot([real(q_rad_2)],[0.0],"pk")
circle_plot(2.0,[real(q_rad_1),0.0],"k")
circle_plot(2.0,[real(q_rad_2),0.0],"k")
elseif flux_type == "CF"
circle_plot(2.0,[real(q_rad_1),0.0])
plot([real(q_rad_1)],[0.0],"pr")

plot([real(q_rad_1)+f[3]],[0.0],"pk")
plot([real(q_rad_2)-f[3]],[0.0],"pk")
circle_plot(2.0,[real(q_rad_1)+f[3],0.0],"k")
circle_plot(2.0,[real(q_rad_2)-f[3],0.0],"k")
end

title(latexstring("$flux_type 2 QP: \$ \\nu=1/3 \$ N=10"))
colorbar()

end
#=
reg_rm = circle_plot(rm)
half_rm = circle_plot(sqrt(2*particles*2))
plot(reg_rm[1],reg_rm[2],"-r",label="1/3")
plot(reg_rm[1],reg_rm[3],"-r")
plot(half_rm[1],half_rm[2],"-c",label="1/2")
plot(half_rm[1],half_rm[3],"-c")
=#
#legend()
#figure()
#hist2D(real.(comp_0q),-imag.(comp_0q),bins=200)
#title("0Q")
end
#=
figure()
array_data_1q = Iterators.flatten([qp1_data[2][i,:] for i in 1:particles])
hist2D(real.(array_data_1q),-imag.(array_data_1q),bins=100)
title("1Q-$rad_choice")
figure()
array_data_2q = Iterators.flatten([qp2_data[2][i,:] for i in 1:particles])
hist2D(real.(array_data_2q),-imag.(array_data_2q),bins=100)
title("2Q-$rad_choice")
#end
=#
#

#x_rad = x_rads[rad_choice]
#qp2_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,2,log_form)
#qp1_data = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n,p,rad_choice,1,log_form)

if false
xs_count = 50
autocorr_length_qp1 = 1
autocorr_length_qp2 = 1
dens_data_1q = get_density(flux_type,particles,n,p,qp1_data,xs_count,autocorr_length_qp1,false)
dens_data_2q = get_density(flux_type,particles,n,p,qp2_data,xs_count,autocorr_length_qp2,false)
#=
qp1_data_left = [qp1_data[1],qp1_data[2].*-1,qp1_data[3]]
qp2_data_left = [qp2_data[1],qp2_data[2].*-1,qp2_data[3]]
dens_data_1q_left = get_density(flux_type,particles,n,p,qp1_data_left,xs_count,autocorr_length_qp1,false)
dens_data_2q_left = get_density(flux_type,particles,n,p,qp2_data_left,xs_count,autocorr_length_qp2,false)
=#
#=
qp1_data[2] .*= -1.0
qp2_data[2] .*= -1.0
dens_data_1q_left = get_density(flux_type,particles,n,p,qp1_data,xs_count,autocorr_length_qp1,false)
dens_data_2q_left = get_density(flux_type,particles,n,p,qp2_data,xs_count,autocorr_length_qp2,false)
=#
#append!(densities_1q,[dens_data_1q])
#append!(densities_2q,[dens_data_2q])
qpart_location_1 = qp2_data[1][2][1]
qpart_location_2 = qp2_data[1][2][2]
end
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

if false
figure()
plot(dens_data_1q[1]./rm,dens_data_1q[2],"-b",label="1Q")
plot(dens_data_2q[1]./rm,dens_data_2q[2],"-r",label="2Q")
#plot(dens_data_1q_left[1]./(-rm),dens_data_1q_left[2],"-b")
#plot(dens_data_2q_left[1]./(-rm),dens_data_2q_left[2],"-r")
#xlim((0.0,1.0))
#ylim((0.00015,0.00075))
#plot(densities_1q[end][1][rad_range_1q[1]:rad_range_1q[end]]./rm,fitcurves_1q[end],"-c")
#plot(densities_2q[end][1][rad_range_2q[1]:rad_range_2q[end]]./rm,fitcurves_2q[end],"-m")
scatter([real(qpart_location_1)]./rm,[dens_data_2q[2][Int(floor(end/2))]],color="g")#[densities_1q[end][2][1]])#,label="QP")
scatter([real(qpart_location_2)]./rm,[dens_data_2q[2][Int(floor(end/2))]],color="g")
#title(latexstring("Density for \$ \\nu=$n/$fill_denom, N=$particles, Rad=$rad_choice\$"))
title("$flux_type Density $rad_choice")
legend(loc="upper left")
end
end
end
#=
exp_pval = (fit_2q.param[2]^2 - fit_1q.param[2]^2)/4
error = sqrt(2*(stderror(fit_1q)[2]*fit_1q.param[2])^2 + 2*(stderror(fit_2q)[2]*fit_2q.param[2])^2)
#println(exp_pval," +/- ",error)
exp_pvals[i] = exp_pval
errors[i] = error
=#

#end

#

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
