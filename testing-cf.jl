using Test

include("cf-wavefunc.jl")

function get_setup(x)
	return [1,x,2,3]
end
#=
@testset "no overlap" begin
	@test get_wavefunc(get_setup(1),4) ≈ 0.0 atol=0.01
	@test get_wavefunc(get_setup(2),4) ≈ 5.0 atol=0.01
	@test get_wavefunc(get_setup(3),4) ≈ 0.0 atol=0.01
end
=#



coords1 = rand(4)
coords2 = rand(4)


coords1[1] = coords1[3]
one_Ji_2 = get_Jis(coords1,2)
@test !isapprox(one_Ji_2,0.0,atol=sqrt(eps()))
one_Ji_3 = get_Jis(coords1,3)
@test isapprox(one_Ji_3,0.0,atol=sqrt(eps()))

coords2[2] = coords2[1]
two_Ji_2 = get_Jis(coords2,2)
@test isapprox(two_Ji_2,0.0,atol=sqrt(eps()))
two_Ji_3 = get_Jis(coords2,3)
@test !isapprox(two_Ji_3,0.0,atol=sqrt(eps()))


comp_coords1 = rand(4) + im*rand(4)
ep = 10^(-5)
comp_coords2 = comp_coords1 .+ 0
comp_coords2[1] += ep

one_ji = get_Jis(comp_coords1,1)
one_jiprime = get_Jiprime(comp_coords1,1)

two_ji = get_Jis(comp_coords2,1)
two_jiprime = get_Jiprime(comp_coords2,1)

num_jiprime = 2*(two_ji - one_ji)/ep
@test isapprox(num_jiprime,one_jiprime,atol=10^(-3))


coords3 = rand(2) + im*rand(2)
three_det_reg = get_wavefunc(coords3)
three_det_log_approx = exp(get_log_det(coords3))
@test isapprox(three_det_reg,three_det_log_approx,atol=10^(-3))





