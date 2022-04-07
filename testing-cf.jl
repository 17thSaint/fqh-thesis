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


n = parse(Int64,ARGS[1])
p = parse(Int64,ARGS[2])
particles = parse(Int64,ARGS[3])


#=
@testset "compare no p" begin
	for i in 1:4
		Ji_nop = get_Jis_nop(coords1,i)
		Ji_p = get_Jis(coords1,i,1)
		@test isapprox(Ji_nop,Ji_p,atol=sqrt(eps()))
		
		Jiprime_nop = get_Jiprime_nop(coords1,i)
		Jiprime_p = get_Jiprime(coords1,i,1)
		@test isapprox(Jiprime_nop,Jiprime_p,atol=sqrt(eps()))
		
		Ji2prime_nop = get_Ji2prime_nop(coords1,i)
		Ji2prime_p = get_Ji2prime(coords1,i,1)
		@test isapprox(Ji2prime_nop,Ji2prime_p,atol=sqrt(eps()))
	end
	for i in 1:4
	for j in 1:4
	nop_elem = get_elem_projection_nop(coords1,i,j,n)
	p_elem = get_elem_projection(coords1,i,j,n,1)
	@test isapprox(nop_elem,p_elem,atol=sqrt(eps()))
	end
	end
end;
=#
if false
@testset "classics" begin

	coords3 = 10*(rand(particles) + im*rand(particles))
	coords4 = 10*(rand(particles) + im*rand(particles))
	coords1 = rand(particles)
	coords2 = rand(particles)
	comp_coords1 = rand(particles) + im*rand(particles)
	ep = 10^(-5)
	comp_coords2 = comp_coords1 .+ 0.0
	comp_coords2[1] += ep
	#
	# testing that when particles overlap Ji goes to zero
	coords1[1] = coords1[3]
	one_Ji_2 = get_Jis(coords1,2,p)
    	#one_Ji_2_nop = get_Jis_nop(coords1,1)
	one_wavefunc = abs2(get_wavefunc(coords1,n,p))
	#one_log_wavefunc = abs2(get_wavefunc_fromlog(coords1,n,p))
	#@test !isapprox(one_wavefunc,0.0,atol=sqrt(eps()))
	#@test !isapprox(one_log_wavefunc,0.0,atol=sqrt(eps()))
	#@test !isapprox(one_Ji_2,0.0,atol=sqrt(eps()))
	one_Ji_3 = get_Jis(coords1,3,p)
	@test isapprox(one_Ji_3,0.0,atol=sqrt(eps()))

	coords2[2] = coords2[1]
	two_Ji_2 = get_Jis(coords2,2,p)
	@test isapprox(two_Ji_2,0.0,atol=sqrt(eps()))
	two_Ji_3 = get_Jis(coords2,3,p)
	@test !isapprox(two_Ji_3,0.0,atol=sqrt(eps()))
	
	# testing that ji prime is the numerical derivative of ji
	one_ji = get_Jis(comp_coords1,1,p)
	one_jiprime = get_Jiprime(comp_coords1,1,p)
	two_ji = get_Jis(comp_coords2,1,p)
	two_jiprime = get_Jiprime(comp_coords2,1,p)
	num_jiprime = (two_ji - one_ji)/ep
	@test isapprox(num_jiprime,one_jiprime,atol=10^(-3))
	
	# testing that ji2prime is numerical derivative of jiprime
	one_jiprime = get_Jiprime(comp_coords1,1,p)
	one_ji2prime = get_Ji2prime(comp_coords1,1,p)
	two_jiprime = get_Jiprime(comp_coords2,1,p)
	two_ji2prime = get_Ji2prime(comp_coords2,1,p)
	num_ji2prime = (two_jiprime - one_jiprime)/ep
	@test isapprox(num_ji2prime,one_ji2prime,atol=10^(-3))
	

	# testing that jiprime is same for reg and log versions
	reg_jiprime = get_Jiprime(coords3,1,p)
	log_jiprime = get_logJiprime(coords3,1,p)
	@test isapprox(reg_jiprime,exp(log_jiprime),atol=10^(-3))
	#
	
	# testing log addition
	rand1 = 10*(rand(Float64)+im*rand(Float64))
	rand2 = 10*(rand(Float64)+im*rand(Float64) )
	logadd = get_log_add(rand1,rand2)
	regadd = log(exp(rand1) + exp(rand2))
	@test isapprox(real(logadd),real(regadd),atol=10^(-2))
	if imag(logadd) > 0 && imag(regadd) > 0
		modpi = (abs(imag(logadd))-abs(imag(regadd)))/pi
	elseif imag(logadd) < 0 && imag(regadd) < 0
		modpi = (abs(imag(logadd))-abs(imag(regadd)))/pi
	else
		modpi = (abs(imag(logadd))+abs(imag(regadd)))/pi
	end
	rounded = round(modpi,digits=0)
	#println(logadd,", ",regadd)
	@test isapprox(modpi,rounded,atol=10^(-1))
	
	# testing that ji2prime is the same reg vs log
	reg_ji2prime = log(get_Ji2prime(coords3,1,p))
	log_ji2prime = get_logJi2prime(coords3,1,p)
	if imag(log_ji2prime) > 0 && imag(reg_ji2prime) > 0
		modpi_ji2prime = (imag(log_ji2prime) - imag(reg_ji2prime))/pi
	elseif imag(log_ji2prime) < 0 && imag(reg_ji2prime) < 0
		modpi_ji2prime = (imag(log_ji2prime) - imag(reg_ji2prime))/pi
	else
		modpi_ji2prime = (abs(imag(log_ji2prime)) + abs(imag(reg_ji2prime)))/pi
	end
	rounded_ji2prime = round(modpi_ji2prime,digits=0)
	@test isapprox(real(reg_ji2prime),real(log_ji2prime),atol=10^(-2))
	@test isapprox(modpi_ji2prime,rounded_ji2prime,atol=10^(-1))
	
	# testing element of matrix same for reg and log
	for i in 1:particles
	for j in 1:particles
	log_element = get_log_elem_proj(coords3,i,j,n,p)
	reg_element = get_elem_projection(coords3,i,j,n,p)
	@test isapprox(reg_element/exp(log_element),1.0,atol=10^(-2))
	end
	end
	
	# testing wavefunction is same reg and log
	reg_wavefunc = get_wavefunc(coords3,n,p)
	log_wavefunc = get_wavefunc_fromlog(coords3,n,p)
	@test isapprox(reg_wavefunc/exp(log_wavefunc),1.0,atol=10^(-3))
	#
	coords_qpart_test = 100 .*start_rand_config(particles,n,p)
	qpart_test = [2,100 .*start_rand_config(2,n,p)]
	# testing wavefunction is same for reg and log with quasiparticles
	#qpart_test = [2,[rand(Float64)+im*rand(Float64),rand(Float64)+im*rand(Float64)]]
	reg_wavefunc = get_wavefunc(coords_qpart_test,n,p,qpart_test)
	log_wavefunc = get_wavefunc_fromlog(coords_qpart_test,n,p,qpart_test)
	@test isapprox(reg_wavefunc,exp(log_wavefunc),atol=10^(-3))
	
	# testing particle overlapping quasiparticle is zero reg form
	coords_qpart_test[1] == qpart_test[2][1] + eps()
	qpart_overlap_wavefunc = get_wavefunc(coords_qpart_test,n,p,qpart_test)
	@test isapprox(abs2(qpart_overlap_wavefunc),0.0,atol=10^(-4))
	
	# testing particle overlapping quasiparticle is zero log form
	qpart_overlap_wavefunc_log = get_wavefunc_fromlog(coords_qpart_test,n,p,qpart_test)
	@test isapprox(exp(2*real(qpart_overlap_wavefunc_log)),0.0,atol=10^(-5))
	
end;
end

if true
@testset "ji-derivs" begin
	local_config = 10 .*(rand(particles) + im.*rand(particles))
	part = 1
	
	nest_deriv_1 = get_nth_deriv_Ji(local_config,part,1)
	og_deriv_1 = get_Jiprime(local_config,part,1)
	@test isapprox(og_deriv_1,nest_deriv_1,atol=sqrt(eps()))
	
	nest_deriv_2 = get_nth_deriv_Ji(local_config,part,2)
	og_deriv_2 = get_Ji2prime(local_config,part,1)
	@test isapprox(og_deriv_2,nest_deriv_2,atol=sqrt(eps()))
	
	nest_log_deriv_1 = get_nth_deriv_Ji(local_config,part,1,true)
	og_log_deriv_1 = get_logJiprime(local_config,part,1)
	@test isapprox(og_deriv_1,nest_deriv_1,atol=sqrt(eps()))
	
	nest_log_deriv_2 = get_nth_deriv_Ji(local_config,part,2,true)
	og_log_deriv_2 = get_logJi2prime(local_config,part,1)
	@test isapprox(og_log_deriv_2,nest_log_deriv_2,atol=sqrt(eps()))
end;
end

if false && n < 2
@testset "origin-div" begin
	rm = sqrt(2*12*(2*p*n+1)/n)
	data_count = 50
	xs = [-0.5*rm + i*(2*0.5*rm)/data_count for i in 1:data_count]
	coords_config = 10*(rand(12) + im*rand(12))
	laugh_wavefunc = fill(0.0,(data_count,data_count))
	cf_wavefunc = fill(0.0,(data_count,data_count))
	for i in 1:length(xs)
		local_x = xs[i]
		for j in 1:length(xs)
			local_y = xs[j]
			#append!(xs_plot,[local_x])
			#append!(ys_plot,[local_y])
			coords_config[1] = local_x - im*local_y
			
			laugh_wavefunc[i,j] = -prob_wavefunc_laughlin(coords_config,2*p*n+1)
			cf_wavefunc[i,j] = 2*real(get_wavefunc_fromlog(coords_config,n,p))
		end
	end
	normed_laugh_wavefunc = laugh_wavefunc./maximum(laugh_wavefunc)
	normed_cf_wavefunc = cf_wavefunc./maximum(cf_wavefunc)
	diff = abs.(normed_laugh_wavefunc-normed_cf_wavefunc)
	for i in 1:length(xs)
		for j in 1:length(xs)
			@test isapprox(diff[i,j],0.0,atol=10^(-5))
		end
	end
end;
end
#






