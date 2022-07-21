using PyPlot,Statistics,CurveFit

function start_rand_config(num_parts,n,p)
	filling = n/(2*p*n+1)
	rm = sqrt(2*num_parts/filling)
	config = [rand(Float64)*rand(-1:2:1)*rm - im*rand(Float64)*rand(-1:2:1)*rm for i in 1:num_parts]
	return config
end

#
np_vals = [[1,1],[1,2],[1,3],[2,1],[2,2],[2,3]]
added_vals = [0.0 for i in 1:6]
coeff_vals = [0.0 for i in 1:6]
for f in 1:6
	n,p = np_vals[f]
	filling = n/(2*p*n+1)
	bot = 2*p*n+1
	parts = [100,200,500,1000,1500]
	avg_dists = [0.0 for i in 1:length(parts)]
	errors =[0.0 for i in 1:length(parts)]
	rms = [0.0 for i in 1:length(parts)]
	for k in 1:length(parts)
	println(k)
	particles = parts[k]
	rm = sqrt(2*particles/filling)
	rand_config = start_rand_config(particles,n,p)
	all_dists = []
	for i in 1:particles
		local_seps = []
		for j in 1:particles
			if i == j
				continue
			end
			append!(local_seps,[abs(rand_config[i] - rand_config[j])])
		end
		append!(all_dists,[mean(sort(local_seps)[1:Int(0.05*particles)])])
	end
	avg_dists[k] = mean(all_dists)
	errors[k] = std(all_dists)
	rms[k] = rm
	end
	#
	vals = linear_fit(rms,avg_dists)
	exp_avg = [rms[i]*vals[2]+vals[1] for i in 1:length(rms)]
	println("$n/$bot:", vals)
	added_vals[f] = vals[1]
	coeff_vals[f] = vals[2]
	plot(rms,exp_avg)
	errorbar(rms,avg_dists,yerr=[errors,errors],fmt="-o",label="$n/$bot")
end
#
added_avg = mean(added_vals)
added_errors = std(added_vals)
coeff_avg = mean(coeff_vals)
coeff_errors = std(coeff_vals)
println("AVG Sep = (",added_avg," +- ",added_errors,") + (",coeff_avg,"+-",coeff_errors,")Rm")
#

"fin"
