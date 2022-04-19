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
#
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


"fin"
