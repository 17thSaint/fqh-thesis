using LinearAlgebra

include("cf-wavefunc.jl")

particles = 4
log_form = true
cf_n, cf_p = 1,1
acc_mat = get_full_acc_matrix(particles)
all_pascs = [get_pascals_triangle(i)[2] for i in 1:particles]
all_derivs = get_deriv_orders_matrix(particles)

for i in 1:10
this_config = start_rand_config(particles,cf_n,cf_p)
cf_wavefunc = get_wavefunc_fromlog(this_config,cf_n,cf_p)
rf_wavefunc = get_rf_wavefunc(this_config,acc_mat,all_pascs,all_derivs,[0,[0]],true)
#=
cf_th = fill(0.0+im*0.0,(particles,particles))
cf_th[1,:] = [1,1]
cf_th[2,:] = [this_config[1],this_config[2]]
cf_th[:,1] .*= (this_config[1]-this_config[2])^2
=#
println(real(exp(rf_wavefunc))/real(exp(cf_wavefunc)))
end
