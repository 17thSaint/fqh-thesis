using Test

include("cf-wavefunc.jl")

function get_setup(x)
	return [1,x,2,3]
end
@testset "no overlap" begin
	@test get_wavefunc(get_setup(1),4) ≈ 0.0 atol=0.01
	@test get_wavefunc(get_setup(2),4) ≈ 0.0 atol=0.01
	@test get_wavefunc(get_setup(3),4) ≈ 0.0 atol=0.01
end
