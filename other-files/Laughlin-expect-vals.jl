function prob_wavefunc(config, m)
	full = 0
	for j = 1:size(config)[1]
		for i = 1:j-1
			full += -2*m*abs(log( dist_btw(config[i],config[j]) ) )
		end
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)
	end
	return full
end

function start_rand_config(num_parts,m)
	rm = sqrt(2*num_parts*m)
	[[rand(Float64)*rand(-1:2:1)*rm,rand(Float64)*rand(-1:2:1)*rm] for i = 1:num_parts]
end

using Cubature
function integrand_norm(x)
	overlap_cut = 1
	config = [[1.0,1.0]]
	polar_config = [ [sqrt(config[i][1]^2+config[i][2]^2),atan(config[i][2]/config[i][1])] for i in 1:size(config)[1]]
	for i in 1:size(config)[1]
		#if x[1] < sqrt(config[i][1]^2+config[i][2]^2) + 0.001 && x[1] > sqrt(config[i][1]^2+config[i][2]^2) - 0.001 && x[2] < atan(config[i][2]/config[i][1]) + 0.001 && x[2] > atan(config[i][2]/config[i][1]) - 0.001
		if x[1] == polar_config[i][1] && x[2] == polar_config[i][2]
			return 0
		end
	end
	return (x[1]^1)*exp(-prob_wavefunc(push!(polar_config,[x[1],x[2]]),3))
end

function integrand_norm_each(x)
	config = [[1.0,0.75]]
	return integrand_norm(x,config)
end


steps = 100
start_rad = 100
end_rad = 1077
normalization = [0.0 for i in 1:steps]
for i in 0:steps-1
	normalization[i+1] = hcubature(integrand_norm,[0,0],[start_rad + i*(end_rad-start_rad)/steps,2*pi])[1]
	#println(normalization[i+1])
end

#println((maximum(normalization)+minimum(normalization))/2)
using PyPlot
plot(normalization)
