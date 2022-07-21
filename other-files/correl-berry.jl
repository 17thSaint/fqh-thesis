include("berry-cf.jl")

function rewrite_data_corr_length(corr_len,wavefunc_data,pos_data)
	new_length_exact = Int(floor((length(wavefunc_data)-1)/corr_len)+1)
	new_length = Int(round(new_length_exact/10)*10)
	new_pos_data = fill(0.0+im*0.0,(size(pos_data)[1],new_length))
	new_wavefunc_data = [0.0+im*0.0 for i in 1:new_length]
	for i in 0:new_length-1
		new_wavefunc_data[i+1] = wavefunc_data[i*corr_len+1]
		new_pos_data[:,i+1] = pos_data[:,i*corr_len+1]
	end
	return new_pos_data,new_wavefunc_data
end

particles = 10
filling = 1/3
rm1 = sqrt(2*particles/filling)
x_rads = [0.01*rm1 + j*(1.29*rm1)/10 for j in 0:9]
acc_mat = get_full_acc_matrix(particles)
full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
full_deriv_ords = get_deriv_orders_matrix(particles)
#println("Got Presets")
flux_type = "RFA"
low_fluxtype = "rfa"
#mc_steps = 1000
which_np_getberry = 2
np_vals = [[1,1],[1,2],[2,1]]
log_form = true
n_berry,p_berry = np_vals[which_np_getberry]
i = 4
mc_steps = 100000

first_1q_data = read_CF_hdf5("$low_fluxtype-data/indiv-comb-$low_fluxtype-data",flux_type,mc_steps,particles,n_berry,p_berry,1,i,1,log_form)
qpart_rad = real(first_1q_data[1][2][1]) + 0.1 - 0.1
full_pos_data_1q = first_1q_data[2] .+ 1.0 .- 1.0
full_wavefunc_data_1q = first_1q_data[3] .+ 1.0 .- 1.0
for j in 2:100
which_file = j

local_1q_data = read_CF_hdf5("$low_fluxtype-data/indiv-comb-$low_fluxtype-data",flux_type,mc_steps,particles,n_berry,p_berry,which_file,i,1,log_form)
#read_comb_CF_hdf5("$low_fluxtype-data",flux_type,particles,n_berry,p_berry,i,1,true)
global full_pos_data_1q = cat(full_pos_data_1q .+ 1.0 .- 1.0,local_1q_data[2],dims=2)
append!(full_wavefunc_data_1q,local_1q_data[3])

end
#radii_data_2q = read_CF_hdf5("$low_fluxtype-data/indiv-comb-$low_fluxtype-data",flux_type,mc_steps,particles,n_berry,p_berry,which,i,2,log_form)
#compl = collect(Iterators.flatten([radii_data_1q[2][i,:] for i in 1:particles]))
#figure()
#hist2D(real.(compl),-imag.(compl),bins=100)
#
all_1q_data = [[1,[qpart_rad+0.0*im]],full_pos_data_1q,full_wavefunc_data_1q]
corr_data = get_autocorr_length(full_wavefunc_data_1q,1)
corr_length = corr_data[1]
indep_1q_data = rewrite_data_corr_length(corr_length,all_1q_data[3],all_1q_data[2])
all_errs_consec = [0.0 for j in 1:6]
all_errs_indep = [0.0 for j in 1:6]
times = [10,50,100,250,500,750]     
for t in 1:length(times)
	top = times[t]
	println(times[t])
	
	for k in 1:10
	global all_1q_data = [[1,[qpart_rad+0.0*im]],full_pos_data_1q,full_wavefunc_data_1q]
	all_indep_data_1q = [[1,[qpart_rad+0.0*im]],indep_1q_data[1],indep_1q_data[2]]
	max = length(all_indep_data_1q[3])-top
	starting_val = rand([1:max;])
	
	all_indep_data_1q[3] = all_indep_data_1q[3][starting_val:starting_val+top-1]
	all_indep_data_1q[2] = all_indep_data_1q[2][:,starting_val:starting_val+top-1]
	
	all_1q_data[3] = all_1q_data[3][starting_val:starting_val+top-1]
	all_1q_data[2] = all_1q_data[2][:,starting_val:starting_val+top-1]
	#println("Corr Length: $corr_length_1q")
	#rads_1q = real(radii_data_1q[1][2][1])
	results_indep = get_expval(flux_type,n_berry,p_berry,all_indep_data_1q,1,-0.001,1,log_form,acc_mat,full_pasc_tri,full_deriv_ords)
	results_consec = get_expval(flux_type,n_berry,p_berry,all_1q_data,1,-0.001,1,log_form,acc_mat,full_pasc_tri,full_deriv_ords)
	println(-100*results_indep[2]/results_indep[1])
	#println("Rez: ",results[1]/qpart_rad^2,", ",results[2]/qpart_rad^2)
	all_errs_indep[t] += -100*results_indep[2]/results_indep[1]/10
	all_errs_consec[t] += -100*results_consec[2]/results_consec[1]/10
	end
	println(all_errs_consec[t])
end
#
plot(times,all_errs_consec,label="Consec")
plot(times,all_errs_indep,label="Indep")
legend()
	
#	berry_calc_1q_consec = get_expval(flux_type,n_berry,p_berry,consec_data,1,-0.001,1,log_form,acc_mat,full_pasc_tri,full_deriv_ords)
	
	#println("Consec: ",berry_calc_1q_consec[1]/qpart_rad^2,", ",berry_calc_1q_consec[2]/qpart_rad^2)
	
	#=
	radii_data_2q[3] = radii_data_2q[3][starting_val:starting_val+top-1]
	radii_data_2q[2] = radii_data_2q[2][:,starting_val:starting_val+top-1]
	corr_length_2q = 1#get_autocorr_length(radii_data_1q[3],1)[1]
	#println("Corr Length: $corr_length_1q")
	rads_2q[i,which_np_getberry] = real(radii_data_2q[1][2][1])
	berry_calc_2q = get_expval(flux_type,n_berry,p_berry,radii_data_2q,1,-0.001,corr_length_1q,log_form,acc_mat,full_pasc_tri,full_deriv_ords)
	=#	

#


"fin"
