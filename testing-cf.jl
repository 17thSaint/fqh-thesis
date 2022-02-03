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

num_jiprime = (two_ji - one_ji)/ep
@test isapprox(num_jiprime,one_jiprime,atol=sqrt(eps()))

@testset "logJi" begin
coords3 = rand(4) + im*rand(4)
three_ji = get_Jis(coords3,1)
log_three_ji = get_logJi(coords3,1)

@test isapprox(three_ji,exp(log_three_ji),atol=sqrt(eps()))
end

