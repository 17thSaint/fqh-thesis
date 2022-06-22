using Statistics
ARGS = "F"
include("read-CF-data.jl")
#include("cf-wavefunc.jl")
include("energy-time-correlation.jl")

function get_expval(vers,n,p,hdf5_data,chosen_qpart,dtheta,samp_freq,log_form=false,acc_matrix = Matrix{Any}(undef,(0,0)),all_pascal = Vector{Vector{Int}}(undef,0),all_derivs = Vector{Vector{Any}}(undef,0))
	qpart_data,pos_data_og,wavefunc_data = hdf5_data
	qpart_og_loc = qpart_data[2][chosen_qpart] + 0.0
	if qpart_data[1] == 1
		qpart_og = [1,[qpart_og_loc]]
	else
		qpart_og = [2,[0.0+im*0.0,qpart_og_loc]]
	end
	if vers == "CF"
		rm = sqrt(2*size(pos_data_og)[1]*(2*p*n+1)/n)
	elseif vers == "RFA"
		rm = sqrt(2*size(pos_data_og)[1]*(2*p*n-1)/n)
	end
	len = length(wavefunc_data)#size(pos_data_og)[2]
	single_qpart_location = qpart_data[2][chosen_qpart]
	qpart_shifted = single_qpart_location*exp(im*dtheta)
	#println(qpart_data[2][chosen_qpart])
	qpart_data[2][chosen_qpart] = qpart_shifted
	#count_origin = 0
	#count_two = 0
	particles = size(pos_data_og)[1]
	#=
	if vers == "RFA"
		acc_matrix = get_full_acc_matrix(particles)
		all_pascal = [get_pascals_triangle(i)[2] for i in 1:particles]
		all_derivs = get_deriv_orders_matrix(particles)
		println("Presets done")
	end
	=#
	number_slices = Int(floor(len/samp_freq))
	calced_vals = [0.0 for i in 1:number_slices]
	pos_data = fill(0.0+im*0.0,(particles,number_slices))
	pos_data[:,1] = pos_data_og[:,1]
	if samp_freq != 1
		for i in 1:number_slices
			pos_data[:,i] = pos_data_og[:,i*samp_freq]  
		end
	else
		pos_data = pos_data_og
	end
	#=
	num_thrds = 10
	slices_per = Int(number_slices/num_thrds)
	slice_sets = [[(k-1)*slices_per+i for i in 1:slices_per] for k in 1:num_thrds]
	Threads.@threads for j in 1:num_thrds
	=#
	for i in 1:number_slices #slice_sets[j]  averaged over a bunch of time slices separated by sampling frequency
			# sampling frequency found from time autocorrelation: ~100
		#if j==1 && i%(slices_per*0.1) == 0
		#    println("Running: ",100*i/slices_per,"%")
	    	#end
		local_config = Complex.(pos_data[:,i])
		#dist_to_qpart_origin = abs.(local_config)# .- single_qpart_location)
		#dist_to_qpart_two = abs.(local_config .- single_qpart_location)
		#dist_to_qpart_shifted = abs.(local_config .- qpart_shifted) 
		#too_close_origin = findall(x->x<=0.001*rm,dist_to_qpart_origin)
		#too_close_two = findall(x->x<=0.001*rm,dist_to_qpart_two)
		#count_origin += length(too_close_origin)
		#count_two += length(too_close_two)
		#if length(too_close_two) > 0 || length(too_close_origin) > 0
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
			if vers == "CF"
				new_wavefunc = get_wavefunc_fromlog(local_config,n,p,qpart_data)
			elseif vers == "RFA"
				new_wavefunc = get_rf_wavefunc(local_config,acc_matrix,all_pascal,all_derivs,qpart_data,log_form)
			end
			#og_wavefunc = get_rf_wavefunc(local_config,acc_matrix,all_pascal,all_derivs,qpart_og,log_form)
			ratio_wavefunc = exp(new_wavefunc - wavefunc_data[i])
			calced_vals[i] = imag(ratio_wavefunc)/dtheta

		else
			if vers == "CF"
				new_wavefunc = get_wavefunc(local_config,n,p,qpart_data)
			elseif vers == "RFA"
				new_wavefunc =  get_rf_wavefunc(local_config,acc_matrix,all_pascal,all_derivs,qpart_data,log_form)
			end
			ratio_wavefunc = new_wavefunc / wavefunc_data[i]
			calced_vals[i] = imag(ratio_wavefunc)/dtheta
		end
		#=
		if isnan(calced_vals[end]) || isinf(calced_vals[end])
			println("Bad Berry: Stopping")
			println(new_wavefunc,", ",wavefunc_data[i],", ",i,", ",det_part)
			throw(DomainError(calced_vals[end],"Irrational"))
		end
		=#
	#end
	
	end
	num_avgd = number_slices - length(findall(calced_vals.==0.0))
	println("Averaged over $num_avgd Slices")
	#
	if num_avgd == number_slices
		final_val = mean(calced_vals)
		val_error = std(calced_vals)
	else
		final_val = Inf
		val_error = 0.0
	end
	#
	return final_val,val_error,calced_vals#,count_two,count_origin
end



#=
particles = 20
filling = 2/5
rm1 = sqrt(2*particles/filling)
#x_rads = [0.01*rm1 + j*(1.29*rm1)/10 for j in 0:9]
#acc_mat = get_full_acc_matrix(particles)
#full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
#full_deriv_ords = get_deriv_orders_matrix(particles)
#println("Got Presets")
flux_type = "CF"
low_fluxtype = "cf"
#mc_steps = 1000
which_np_getberry = 3
top = 10
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
n_berry,p_berry = np_vals[which_np_getberry]     
time_length = 1
for i in 1:top
		starting_val = rand([1:(40100-time_length);])
		
		#
		radii_data_1q = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n_berry,p_berry,30+i,1,true)
		radii_data_1q[3] = radii_data_1q[3][starting_val:starting_val+time_length-1]
		radii_data_1q[2] = radii_data_1q[2][:,starting_val:starting_val+time_length-1]
		corr_length_1q = 1
		rads_1q[i,which_np_getberry] = real(radii_data_1q[1][2][1])
		berry_calc_1q = get_expval(flux_type,n_berry,p_berry,radii_data_1q,1,-0.001,corr_length_1q,log_form)#,acc_mat,full_pasc_tri,full_deriv_ords)
		berries_1q[i,which_np_getberry] = berry_calc_1q[1]
		errors_1q[i,which_np_getberry] = berry_calc_1q[2]
		#=
		
		radii_data_2q = read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n_berry,p_berry,i+2,2,true)
		radii_data_2q[3] = radii_data_2q[3][starting_val:starting_val+time_length-1]
		radii_data_2q[2] = radii_data_2q[2][:,starting_val:starting_val+time_length-1]
		corr_length_2q = 1
		rads_2q[i,which_np_getberry] = real(radii_data_2q[1][2][1])
		berry_calc_2q = get_expval(flux_type,n_berry,p_berry,radii_data_2q,1,-0.001,corr_length_2q,log_form)#,acc_mat,full_pasc_tri,full_deriv_ords)
		berries_2q[i,which_np_getberry] = berry_calc_2q[1]
		errors_2q[i,which_np_getberry] = berry_calc_2q[2]
		=#
		
end
#
theory_1q = (rads_1q[:,which_np_getberry].^2).*(filling/4)
#interaction_part_2q = berries_2q[:,which_np_getberry] .+ theory_1q
#errorbar(rads_2q[:,which_np_getberry],berries_2q[:,which_np_getberry],yerr=[errors_2q[:,which_np_getberry],errors_2q[:,which_np_getberry]],fmt="-o",label="2Q")
errorbar(rads_1q[:,which_np_getberry],-berries_1q[:,which_np_getberry],yerr=[errors_1q[:,which_np_getberry],errors_1q[:,which_np_getberry]],fmt="-o",label="1Q")
plot(rads_1q[:,which_np_getberry],theory_1q,label="TH")
#errorbar(rads_2q[:,which_np_getberry],interaction_part_2q,yerr=[errors_2q[:,which_np_getberry],errors_2q[:,which_np_getberry]],fmt="-o",label="2Q")
#plot(rads_2q[:,which_np_getberry],[filling for i in 1:top],label="TH")
legend()
#println(berries_1q[4,2]/rads_1q[4,2]^2,", ",errors_1q[4,2]/rads_1q[4,2]^2)

=#

#
# 16 parts 2 qparts 1000 mcs 1/3 rad 7 ver 1, has particle stuck near origin qpart causing bad berry calc
# rad 2 ver 2 has first particle stuck too, same rad 5 and rad 9
#=
which_nps_plotting = 2
n_plot,p_plot = np_vals[which_nps_plotting]
th_coeff_vals = [1,1,1]#[2,1,1]
th_coeff = th_coeff_vals[which_nps_plotting]
fill_denom = 2*p_plot*n_plot - 1
rm = sqrt(2*particles*fill_denom/n_plot)

#plot(rads_2q[:,which_nps_plotting],origin_counts[:,which_nps_plotting],"-p",label="OG")
#plot(rads_2q[:,which_nps_plotting],two_counts[:,which_nps_plotting],"-p",label="Two")
#legend()

#berries_1q_n10_all = [0.0 0.0 0.0; 0.0 -1.0257318045919919 0.0; 0.0 -1.7646445517636962 0.0; 0.0 -3.677582158381534 0.0; 0.0 -8.925330663309417 0.0]

errorbar(rads_1q[2:5,which_nps_plotting],-berries_1q[2:5,which_nps_plotting],yerr=[errors_1q[2:5,which_nps_plotting],errors_1q[2:5,which_nps_plotting]],fmt="-o")
#plot(rads_1q[2:5,which_nps_plotting],-berries_1q_n10_all[2:5,which_nps_plotting],"-p")

theory_here_1q = rads_1q[2:5,which_nps_plotting].*rads_1q[2:5,which_nps_plotting]./2#(th_coeff*(2*p_plot*n_plot-1))

plot(rads_1q[2:5,which_nps_plotting],theory_here_1q .+ 1/2,label="1/2")
plot(rads_1q[2:5,which_nps_plotting],theory_here_1q.*(2/3) .+ 1/3,label="1/3")
legend()
xlabel("QP Radius")
title(latexstring("Berry Phase = \$ \\vert \\eta \\vert^2 \\nu \$"))
=#
#theory_here_2q_m = rads_2q[:,which_nps_plotting].*rads_2q[:,which_nps_plotting]./(th_coeff*(2*p*n+1)) - [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,which_nps_plotting])]
#theory_here_2q_p = rads_2q[:,which_nps_plotting].*rads_2q[:,which_nps_plotting]./(th_coeff*(2*p*n+1)) + [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,which_nps_plotting])]
#figure()
#plot(rads_1q[:,which_nps_plotting]./rm,theory_here_1q,"-k",label="TH")# $n / $fill_denom")
#scatter(rads_2q[:,which_nps_plotting],theory_here_2q_m,c="r",label="~TH 2Q M")# $n / $fill_denom")
#scatter(rads_2q[:,which_nps_plotting],theory_here_2q_p,c="g",label="~TH 2Q P")# $n / $fill_denom")
#errorbar(rads_1q[:,which_nps_plotting]./rm,-berries_1q[:,which_nps_plotting],yerr=[errors_1q[:,which_nps_plotting],errors_1q[:,which_nps_plotting]],fmt="-o",c="r",label="1Q")# $n / $fill_denom")
#errorbar(rads_2q[:,which_nps_plotting]./rm,berries_2q[:,which_nps_plotting],yerr=[errors_2q[:,which_nps_plotting],errors_2q[:,which_nps_plotting]],fmt="-o",label="2Q")# $n / $fill_denom")
#xlabel("Quasiparticle Radius")
#ylabel("Berry Phase")
#title(latexstring("Berry Phase for \$ \\nu = $n_plot/$fill_denom \$, N=$particles"))
#legend()
#=
# compare shift in berry phase with 2QP
#diff = -berries_1q[:,j]+berries_2q[:,j]
diff_vals = theory_here_1q .+ berries_2q[:,which_nps_plotting]#berries_1q[:,which_nps_plotting]-berries_2q[:,which_nps_plotting]
#diff_errors = sqrt.(errors_1q[:,which_nps_plotting].^2 + errors_2q[:,which_nps_plotting].^2)
diff_errors = sqrt.(errors_2q[:,which_nps_plotting].^2)
theory_diff = [1/3,1/3]#[2*p_plot/(th_coeff*(2*p_plot*n_plot + 1)),-2*p_plot/(th_coeff*(2*p_plot*n_plot + 1))]
figure()
errorbar(rads_2q[:,which_nps_plotting],diff_vals,yerr=[diff_errors,diff_errors],fmt="-o")
plot(rads_2q[:,which_nps_plotting],[theory_diff[1] for i in 1:length(rads_2q[:,which_nps_plotting])],label="TH P")
plot(rads_2q[:,which_nps_plotting],[theory_diff[2] for i in 1:length(rads_2q[:,which_nps_plotting])],label="TH M")
legend()
xlabel("Quasiparticle Radius")
title(latexstring("Shift from Enclosed CFQP at \$ \\nu=$n_plot/$fill_denom \$"))
=#







"fin"
