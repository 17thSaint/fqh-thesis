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
	diff = logdet_multip - (logdet_reg + expected_shift)
	#println(round(imag(diff)/pi,digits=0),", ",imag(diff)/pi)
	#if !isapprox(round(imag(diff)/pi,digits=0),imag(diff)/pi,atol=10^(-1))
	#	println(imag(diff)/pi)
	#end
	return diff
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

#
particles = 12
mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = 1,1#np_vals[which]
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
xs_plot = []
ys_plot = []
data_count = 30
starting_config = start_rand_config(particles,n,p)#1 .*qp1_data[2][:,3800]
parts_xs = real.(starting_config[2:end])
parts_ys = -imag.(starting_config[2:end])
#scatter3D(parts_xs,parts_ys,[1.0 for i in 1:particles-1],c="r")
#xs = [-0.01*rm + i*(2*0.01*rm)/data_count for i in 1:data_count]
xs = [-0.01*rm + i*(2*0.01*rm)/data_count for i in 1:data_count]
#exp_prob = []
#reg_prob = []
logjast = fill(0.0+im*0.0,(data_count,data_count))
logj1 = fill(0.0+im*0.0,(data_count,data_count))
logslater = fill(0.0+im*0.0,(data_count,data_count))
#part1 = []
exp_prob_CF_sep = fill(0.0,(data_count,data_count))
exp_prob_CF_reg = fill(0.0,(data_count,data_count))
exp_prob_Laugh = fill(0.0,(data_count,data_count))

diffs = fill(0.0+im*0.0,(data_count,data_count))
#givematrix = [fill(0.0+im*0.0,(particles,particles))]
trial_vals = []
for i in 1:length(xs)
	local_x = xs[i]
	println(i/length(xs))
	for j in 1:length(xs)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		radius = abs(starting_config[1])
		#=
		logjast[i,j] = get_logJastrowfull(starting_config)
		logj1[i,j] = get_logJi(starting_config,1)
		logslater[i,j] = log(starting_config[1])
		=#
		vals_CF = get_wavefunc_fromlog(starting_config,n,p,qpart)
		append!(trial_vals,[vals_CF,starting_config])
		
		sep_matrix = vals_CF[4]
		exp_shift = [p*get_logJi(starting_config,k) for k in 1:particles]
		diff_local = get_logdet_diff(sep_matrix,exp_shift)
		diffs[i,j] = diff_local
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
#
if true
figure()
imshow(real.(diffs))#./maximum(exp_prob_CF))
title("Real Diffs")
colorbar()
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
