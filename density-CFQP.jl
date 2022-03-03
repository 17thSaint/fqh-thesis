using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")

function get_density(particles,n,p,hdf5_data,data_count)
	rm = sqrt(2*particles*(2*p*n+1)/n)
	xs = [i*1.3*rm/data_count + 0.5*1.3*rm/data_count for i in 0:data_count-1]
	qpart_data,pos_data,wavefunc_data = hdf5_data
	dens_val = [0.0 for i in 1:data_count]
	errors = [0.0 for i in 1:data_count]
	for k in 1:data_count
	    edge = xs[k]
	    local_denss = [0.0 for w in 1:length(wavefunc_data)]
	    for i in 1:length(wavefunc_data)
		radii = abs.(pos_data[:,i])
		local_count = 0
		for j in 1:length(radii)
		    if radii[j] <= edge
		        local_count += 1
		    end
		end
		local_denss[i] = local_count/(pi*edge^2)
	    end
	    dens_val[k] = mean(local_denss)
	    errors[k] = std(local_denss)
	end
	return xs,dens_val,errors
end


particles = 8
mc_steps = 100000
qpart_count = 1
np_vals = [[1,1],[1,2],[2,1]]
which = 4
j = 1
n,p = np_vals[j]
fill_denom = 2*n*p + 1
hdf5_data = read_CF_hdf5("cf-data",mc_steps,particles,n,p,which,qpart_count)
dens_data = get_density(particles,n,p,hdf5_data,50)
qpart_location = hdf5_data[1][2][1]

plot(dens_data[1],dens_data[2],"-p",label="$n/$fill_denom")
scatter([real(qpart_location)],[dens_data[2][1]],c="r")
legend()
#println(qpart_location)

"fin"
