include("cf-wavefunc.jl")
using Statistics,Dates,PyPlot,Test

#=
function get_times(particles)
	max_rej = Int(floor(0.2*particles))
	all_parts = [i for i in 2:particles]
	this_config = start_rand_config(particles,1,1)
	rej_parts = [rand(all_parts) for i in 1:max_rej]
	while length(unique(rej_parts)) < max_rej
		append!(rej_parts,[rand(all_parts)])
	end
	unique!(rej_parts)
	allowed_parts = deleteat!(all_parts.+0,sort(rej_parts).-1)

	time_rej = @time get_Jis(this_config,1,rej_parts)
	time_acc = @time get_Jis_acc(this_config,1,allowed_parts)
	return time_rej,time_acc
end
=#
#=
function get_wavefunc_times(particles)
	this_config = start_rand_config(particles,1,1)
	acc_mat = get_allowed_sets_matrix(particles)
	println("Made Acc Mat")
	full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
	full_derivs = get_deriv_orders_matrix(particles)
	all_times = [0.0 for i in 1:10]
	compile_first = get_rf_wavefunc(this_config,acc_mat,full_pasc_tri,full_derivs,[0,[0]],true)
	for i in 1:10
		time_start = now()
		get_rf_wavefunc(this_config,acc_mat,full_pasc_tri,full_derivs,[0,[0]],true)
		time_end = now()
		all_times[i] = (time_end-time_start).value
	end
	avg_time = mean(all_times)
	errors = std(all_times)
	return avg_time,errors
end

parts_vals = [i for i in 2:10]
time_vals = [0.0 for i in 1:length(parts_vals)]
time_errors = [0.0 for i in 1:length(parts_vals)]
for i in 1:length(parts_vals)
	parts = parts_vals[i]
	println(parts)
	rezz = get_wavefunc_times(parts)
	time_vals[i] = rezz[1]
	time_errors[i] = rezz[2]
end
errorbar(parts_vals,time_vals,yerr=[time_errors,time_errors],fmt="-o")

particles = 9
this_config = start_rand_config(particles,1,1)
full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
full_derivs = get_deriv_orders_matrix(particles)
for i in 1:10
println("Inside")
@time get_rf_wavefunc_internal(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
println("External")
@time begin
acc_mat = get_allowed_sets_matrix(particles)
get_rf_wavefunc(this_config,acc_mat,full_pasc_tri,full_derivs,[0,[0]],true)
end
end
=#
#=
part_choices = [2,4,6,8]
time_vals_ext = [0.0 for i in 1:length(part_choices)]
time_vals_int = [0.0 for i in 1:length(part_choices)]
time_vals_sing = [0.0 for i in 1:length(part_choices)]
time_vals_column = [0.0 for i in 1:length(part_choices)]
error_vals_ext = [0.0 for i in 1:length(part_choices)]
error_vals_int = [0.0 for i in 1:length(part_choices)]
error_vals_sing = [0.0 for i in 1:length(part_choices)]
error_vals_column = [0.0 for i in 1:length(part_choices)]
for i in 1:length(part_choices)
particles = part_choices[i]
my_orders = [p for p in 1:particles-1]
this_config = start_rand_config(particles,1,1)
full_pasc_tri = [get_pascals_triangle(k)[2] for k in 1:particles]
full_derivs = get_deriv_orders_matrix(particles)
acc_mat = Matrix{Vector{Vector{Int}}}(undef,particles,particles-1)
#
for which_part in 1:particles
	acc_mat[which_part,:] = read_acc_matrix_data("acc-matrix-data",particles,which_part,my_orders)
end
#
println(particles)
time_vals_sing_local = [0.0 for j in 1:5]
time_vals_ext_local = [0.0 for j in 1:5]
time_vals_column_local = [0.0 for j in 1:5]
time_vals_int_local = [0.0 for j in 1:5]
#read_acc_matrix_data("acc-matrix-data",particles,1,1)
#
get_rf_wavefunc_indivelemcall(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
get_rf_wavefunc_columncall(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
get_rf_wavefunc(this_config,acc_mat,full_pasc_tri,full_derivs,[0,[0]],true)
get_rf_wavefunc_internal(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
#
for j in 1:5

time_start_sing = now()
ind = get_rf_wavefunc_indivelemcall(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
#for ps in 1:particles
#read_acc_matrix_data("acc-matrix-data",particles,ps,1)
#end
time_end_sing = now()
time_vals_sing_local[j] = (time_end_sing - time_start_sing).value
#
time_start_ext = now()
ext = get_rf_wavefunc(this_config,acc_mat,full_pasc_tri,full_derivs,[0,[0]],true)
time_end_ext = now()
time_vals_ext_local[j] = (time_end_ext - time_start_ext).value

time_start_int = now()
int = get_rf_wavefunc_internal(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
time_end_int = now()
time_vals_int_local[j] = (time_end_int - time_start_int).value
#
time_start_column = now()
col = get_rf_wavefunc_columncall(this_config,full_pasc_tri,full_derivs,[0,[0]],true)
#read_acc_matrix_data("acc-matrix-data",particles,my_parts,1)
time_end_column = now()
time_vals_column_local[j] = (time_end_column - time_start_column).value

@test isapprox(ind,ext,atol=sqrt(eps()))
@test isapprox(ind,col,atol=sqrt(eps()))
@test isapprox(int,col,atol=sqrt(eps()))
end
time_vals_sing[i] = mean(time_vals_sing_local)
error_vals_sing[i] = std(time_vals_sing_local)
time_vals_ext[i] = mean(time_vals_ext_local)
error_vals_ext[i] = std(time_vals_ext_local)
time_vals_column[i] = mean(time_vals_column_local)
error_vals_column[i] = std(time_vals_column_local)
time_vals_int[i] = mean(time_vals_int_local)
error_vals_int[i] = std(time_vals_int_local)

end

errorbar(part_choices,time_vals_sing,yerr=[error_vals_sing,error_vals_sing],fmt="-o",label="Single")
errorbar(part_choices,time_vals_ext,yerr=[error_vals_ext,error_vals_ext],fmt="-o",label="External")
errorbar(part_choices,time_vals_column,yerr=[error_vals_column,error_vals_column],fmt="-o",label="Column")
errorbar(part_choices,time_vals_int,yerr=[error_vals_int,error_vals_int],fmt="-o",label="Internal")
legend()
xlabel("Particles")
title("Time for Single Wavefunction Calculation")
=#
#=
include("mc-cf.jl")


particles = 6
mc_steps = 100
log_form = true
np_vals = [[1,1],[1,2],[2,1]]
which_np = 1#parse(Int64,ARGS[1])
n,p = np_vals[which_np]
fill_denom = 2*n*p + 1
filling = n/(2*p*n+1)
rm = sqrt(2*particles/filling)
step_size = 0.4 + 0.175*rm
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
k = 5#parse(Int64,ARGS[1])
rad_choice = k
x_rad = x_rads[rad_choice]
qpart_choices = [[0,[0.0]],[1,[x_rad+im*0.0]],[2,[x_rad+im*0.0,0.0+im*0.0]]]
qpart_selected = 2#parse(Int64,ARGS[1])
qpart = qpart_choices[qpart_selected]
main("RFA",n,p,5,particles,step_size,k,qpart,log_form)
main("RFA",n,p,5,particles,step_size,k,qpart,log_form,true)

howmany = 10
time_vals_one_local = [0.0 for j in 1:howmany]
time_vals_both_local = [0.0 for j in 1:howmany]
for i in 1:howmany
	println(i/howmany)
	time_start_one = now()
	main("RFA",n,p,mc_steps,particles,step_size,k,qpart,log_form)
	time_end_one = now()
	time_vals_one_local[i] = (time_end_one - time_start_one).value
	time_start_both = now()
	main("RFA",n,p,mc_steps,particles,step_size,k,qpart,log_form,true)
	time_end_both = now()
	time_vals_both_local[i] = (time_end_both - time_start_both).value
end
one_val = mean(time_vals_one_local)
one_error = std(time_vals_one_local)
both_val = mean(time_vals_both_local)
both_error = std(time_vals_both_local)
println("One: $one_val +/- $one_error, Both: $both_val +/- $both_error")
=#
#=
parts = [4,5,6,7,8]
times_parts = [0.0 for i in 1:length(parts)]
errors_parts = [0.0 for i in 1:length(parts)]
for j in 1:length(parts)
particles = parts[j]
acc_mat = get_full_acc_matrix(particles)
all_pascs = [get_pascals_triangle(i)[2] for i in 1:particles]
all_derivs = get_deriv_orders_matrix(particles)
this_config = start_rand_config(particles,1,1)
howmany = 10
time_vals_one_local = [0.0 for j in 1:howmany]
get_rf_wavefunc(this_config,acc_mat,all_pascs,all_derivs,[0,[0.0]],true)
for i in 1:howmany
	time_start_one = now()
	get_rf_wavefunc(this_config,acc_mat,all_pascs,all_derivs,[0,[0.0]],true)
	time_end_one = now()
	time_vals_one_local[i] = (time_end_one - time_start_one).value
end
one_val = mean(time_vals_one_local)
one_error = std(time_vals_one_local)

time_full = particles*one_val
error_full = particles*one_error

times_parts[j] = time_full
errors_parts[j] = error_full
end

#errorbar(parts,times_parts,yerr=[errors_parts,errors_parts],fmt="-o")
power_start = 7.6
power_end = 7.8
powers = [power_start + i*(power_end-power_start)/10 for i in 0:10]
howclose = [0.0 for i in 0:10]
for i in 1:length(powers)
	pow = powers[i]
	rez = times_parts./(parts.^pow)
	perc_error = 100*std(rez)/mean(rez)
	howclose[i] = perc_error
end
plot(powers,howclose)
=#
function estim_time(num_parts,mc_steps)
	n4 = 0.4*mc_steps
	estim_time = ((num_parts/4)^7.7)*n4/1000
	return estim_time
end
#=
howmany = 10
parts = [4,5,6,7,8,9]
times_new = [0.0 for i in 1:length(parts)]
errs_new = [0.0 for i in 1:length(parts)]
times_prev = [0.0 for i in 1:length(parts)]
errs_prev = [0.0 for i in 1:length(parts)]
for particles in parts
println(particles)
input = [true,particles-1,get_all_acc_sets(particles-2,1,particles)]
times_prev_local = [0.0 for i in 1:howmany]
times_new_local = [0.0 for i in 1:howmany]
errs_prev_local = [0.0 for i in 1:howmany]
errs_new_local = [0.0 for i in 1:howmany]
get_allowed_sets_matrix(particles,true)
get_allowed_sets_matrix(particles)
#get_all_acc_sets(particles-1,1,particles,input)
#get_all_acc_sets(particles-1,1,particles)
for i in 1:howmany
	time_start_prev = now()
	get_allowed_sets_matrix(particles,true)
	#get_all_acc_sets(particles-1,1,particles,input)
	time_end_prev = now()
	times_prev_local[i] = (time_end_prev - time_start_prev).value
	#
	time_start_new = now()
	get_allowed_sets_matrix(particles)
	#get_all_acc_sets(particles-1,1,particles)
	time_end_new = now()
	times_new_local[i] = (time_end_new - time_start_new).value
	#
end
times_prev[particles-parts[1]+1] = mean(times_prev_local)
errs_prev[particles-parts[1]+1] = std(times_prev_local)

times_new[particles-parts[1]+1] = mean(times_new_local)
errs_new[particles-parts[1]+1] = std(times_new_local)

end
errorbar(parts,times_new,yerr=[errs_new,errs_new],fmt="-o",label="New")
errorbar(parts,times_prev,yerr=[errs_prev,errs_prev],fmt="-o",label="Prev")
legend()
=#
#figure()
#plot(parts,times_new./times_prev)
#title("Speed-up Factor: Acc Matrix w/ Prev Compressed Data")
power_start = 4.0
power_end = 6.0
powers = [power_start + i*(power_end-power_start)/10 for i in 0:10]
howclose = [0.0 for i in 0:10]
for i in 1:length(powers)
	pow = powers[i]
	rez = (times_new./times_prev)./(parts.^pow)
	perc_error = 100*std(rez)/mean(rez)
	howclose[i] = perc_error
end
plot(powers,howclose)


"fin"
