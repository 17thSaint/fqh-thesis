using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")
include("energy-time-correlation.jl")

function get_expval(n,p,hdf5_data,chosen_qpart,dtheta,samp_freq,log_form=false)
	qpart_data,pos_data_og,wavefunc_data = hdf5_data
	rm = sqrt(2*size(pos_data_og)[1]*(2*p*n+1)/n)
	len = length(wavefunc_data)
	single_qpart_location = qpart_data[2][chosen_qpart]
	qpart_shifted = single_qpart_location*exp(im*dtheta)
	qpart_data[2][chosen_qpart] = qpart_shifted
	calced_vals = []#[0.0 for i in 1:len]
	count_origin = 0
	count_two = 0
	
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

	for i in 1:number_slices # averaged over a bunch of time slices separated by sampling frequency
			# sampling frequency found from time autocorrelation: ~100
		if i%(number_slices*0.1) == 0
		    println("Running: ",100*i/number_slices,"%")
	    	end
		local_config = pos_data[:,i]
		#dist_to_qpart_origin = abs.(local_config)# .- single_qpart_location)
		#dist_to_qpart_two = abs.(local_config .- single_qpart_location)
		#dist_to_qpart_shifted = abs.(local_config .- qpart_shifted) 
		#too_close_origin = findall(x->x<=0.1,dist_to_qpart_origin)
		#too_close_two = findall(x->x<=0.001*rm,dist_to_qpart_two)
		#count_origin += length(too_close_origin)
		#count_two += length(too_close_two)
		#if length(too_close_two) > 0 #|| length(too_close_origin) > 0
		#	println("Too Close: $dist_to_qpart_two")
		#	continue
		#end
		if log_form
			#=
			new_wavefunc = get_wavefunc_fromlog(local_config,n,p,qpart_data)
			ratio_wavefunc = new_wavefunc - wavefunc_data[i]
			#println(exp(real(ratio_wavefunc)),", ",sin(imag(ratio_wavefunc)),", ",i)
			calced_vals[i] = exp(real(ratio_wavefunc))*sin(imag(ratio_wavefunc))/dtheta
			#println(calced_vals[i])
			=#
			new_wavefunc = get_wavefunc_fromlog(local_config,n,p,qpart_data)
			ratio_wavefunc = exp(new_wavefunc - wavefunc_data[i])
			#calced_vals[i] = imag(ratio_wavefunc)/dtheta
			append!(calced_vals,[imag(ratio_wavefunc)/dtheta])
		else
			new_wavefunc = get_wavefunc(local_config,n,p,qpart_data)
			ratio_wavefunc = new_wavefunc / wavefunc_data[i]
			#calced_vals[i] = imag(ratio_wavefunc)/dtheta
			append!(calced_vals,[imag(ratio_wavefunc)/dtheta])
		end
		if isnan(calced_vals[end]) || isinf(calced_vals[end])
			println("Bad Berry: Stopping")
			println(new_wavefunc,", ",wavefunc_data[i],", ",i,", ",det_part)
			throw(DomainError(calced_vals[end],"Irrational"))
		end
		#
	end
	num_avgd = length(calced_vals)
	println("Averaged over $num_avgd Slices: Removed OG=$count_origin, Two=$count_two")
	if num_avgd > 0
		final_val = mean(calced_vals)
		val_error = std(calced_vals)
	else
		final_val = Inf
		val_error = 0.0
	end
	return final_val,val_error#,count_two,count_origin
end


particles = 16
#mc_steps = 1000
top = 5
#
rads_1q = fill(0.0,(top,3))
berries_1q = fill(0.0,(top,3))
errors_1q = fill(0.0,(top,3))
rads_2q = fill(0.0,(top,3))
berries_2q = fill(0.0,(top,3))
errors_2q = fill(0.0,(top,3))
np_vals = [[1,1],[1,2],[2,1]]
origin_counts = fill(0.0,(top,3))
two_counts = fill(0.0,(top,3))
log_form = true
for j in 2:2
        n,p = np_vals[j]
        for i in 3:5
		#
		radii_data_1q = read_comb_CF_hdf5("cf-data",particles,n,p,i,1,true)
		corr_length_1q = 1#get_autocorr_length(radii_data_1q[3],1)[1]
		println("Corr Length: $corr_length_1q")
		rads_1q[i,j] = real(radii_data_1q[1][2][1])
		berry_calc_1q = get_expval(n,p,radii_data_1q,1,-0.001,corr_length_1q,log_form)
		berries_1q[i,j] = berry_calc_1q[1]
		errors_1q[i,j] = berry_calc_1q[2]
		#two_counts[i,j] = berry_calc_1q[3]
		#
		radii_data_2q = read_comb_CF_hdf5("cf-data",particles,n,p,i,2,true)
		#radii_data_2q = read_CF_hdf5("cf-data",100000,particles,n,p,1,i,2,log_form)
		corr_length_2q = 1#get_autocorr_length(radii_data_2q[3],1)[1]
		println("Corr Length: $corr_length_2q")
		rads_2q[i,j] = real(radii_data_2q[1][2][1])
		berry_calc_2q = get_expval(n,p,radii_data_2q,1,-0.001,corr_length_2q,log_form)
		berries_2q[i,j] = berry_calc_2q[1]
		errors_2q[i,j] = berry_calc_2q[2]
		#origin_counts[i,j] = berry_calc_2q[3]
		#two_counts[i,j] = berry_calc_2q[4]
		#
        end
end
#

# 16 parts 2 qparts 1000 mcs 1/3 rad 7 ver 1, has particle stuck near origin qpart causing bad berry calc
# rad 2 ver 2 has first particle stuck too, same rad 5 and rad 9

which = 2
n,p = np_vals[which]
th_coeff_vals = [2,1,1]
th_coeff = th_coeff_vals[which]
fill_denom = 2*p*n + 1
rm = sqrt(2*particles*fill_denom/n)

#plot(rads_2q[:,which],origin_counts[:,which],"-p",label="OG")
#plot(rads_2q[:,which],two_counts[:,which],"-p",label="Two")
#legend()

#= direct compare berry phase
theory_here_1q = rads_1q[:,which].*rads_1q[:,which]./(th_coeff*(2*p*n+1))
#
#theory_here_2q_m = rads_2q[:,which].*rads_2q[:,which]./(th_coeff*(2*p*n+1)) - [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,which])]
#theory_here_2q_p = rads_2q[:,which].*rads_2q[:,which]./(th_coeff*(2*p*n+1)) + [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,which])]
plot(rads_1q[:,which]./rm,theory_here_1q,"-k",label="TH")# $n / $fill_denom")
#scatter(rads_2q[:,which],theory_here_2q_m,c="r",label="~TH 2Q M")# $n / $fill_denom")
#scatter(rads_2q[:,which],theory_here_2q_p,c="g",label="~TH 2Q P")# $n / $fill_denom")
errorbar(rads_1q[:,which]./rm,-berries_1q[:,which],yerr=[errors_1q[:,which],errors_1q[:,which]],fmt="-o",c="r",label="1Q")# $n / $fill_denom")
errorbar(rads_2q[:,which]./rm,-berries_2q[:,which],yerr=[errors_2q[:,which],errors_2q[:,which]],fmt="-o",label="2Q")# $n / $fill_denom")
xlabel("Quasiparticle Radius")
ylabel("Berry Phase")
title(latexstring("Berry Phase for \$ \\nu = $n/$fill_denom \$, N=$particles"))
legend()
=#
j = 1
# compare shift in berry phase with 2QP
#diff = -berries_1q[:,j]+berries_2q[:,j]
diff = -berries_1q[:,j]+berries_2q[:,j]
#diff_errors = sqrt.(errors_1q[:,j].^2 + errors_2q[:,j].^2)
diff_errors = sqrt.(errors_2q[:,j].^2)
theory_diff = [2*p/(th_coeff*(2*p*n + 1)),-2*p/(th_coeff*(2*p*n + 1))]
errorbar(rads_2q[:,j],diff,yerr=[diff_errors,diff_errors],fmt="-o")
plot(rads_2q[:,j],[theory_diff[1] for i in 1:length(rads_2q[:,j])],label="TH P")
plot(rads_2q[:,j],[theory_diff[2] for i in 1:length(rads_2q[:,j])],label="TH M")
legend()
xlabel("Quasiparticle Radius")
title(latexstring("Shift from Enclosed CFQP at \$ \\nu=$n/$fill_denom \$"))
#







"fin"
