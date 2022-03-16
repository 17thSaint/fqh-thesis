using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")

function get_density(particles,n,p,hdf5_data,data_count)
	rm = sqrt(2*particles*(2*p*n+1)/n)
	xs = [i*1.3*rm/data_count + 0.5*1.3*rm/data_count for i in 0:data_count-1]
	qpart_data,pos_data,wavefunc_data = hdf5_data
	change_counts = [0.0 for i in 1:data_count]
	total_possible_counts = particles*length(wavefunc_data)
	#errors = [0.0 for i in 1:data_count]
	previous_count = 0
	for k in 1:data_count
	    edge = xs[k]
	    #
	    local_count = 0
	    for i in 1:length(wavefunc_data)
	    	xpos = real.(pos_data[1:end,i])
		ypos = -imag.(pos_data[1:end,i])
		for j in 1:length(xpos)
		    if xpos[j] <= edge && xpos[j] >= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
		        local_count += 1
		    end
		end
	    end
	    
	    change_counts[k] = (local_count - previous_count)/total_possible_counts
	    previous_count = local_count
	    #=
	    local_denss = [0.0 for q in 1:length(wavefunc_data)]
	    for i in 1:length(wavefunc_data)
		xpos = real.(pos_data[:,i])
		ypos = -imag.(pos_data[:,i])
		local_count = 0
		for j in 1:length(xpos)
		    if xpos[j] <= edge && xpos[j] >= 0.0 && ypos[j] <= 0.1*rm && ypos[j] > -0.1*rm
		        local_count += 1
		    end
		end
		local_denss[i] = local_count#/(pi*edge^2)
	    end
	    dens_val[k] = mean(local_denss)
	    errors[k] = std(local_denss)
	    =#
	end
	return xs,change_counts
end

#=
particles = 16
mc_steps = 100000
qpart_count = 2
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = np_vals[which]
#
fill_denom = 2*n*p + 1
rm = sqrt(2*particles*fill_denom/n)
#for i in 5:5
rad_choice = 5
#hdf5_data = read_CF_hdf5("cf-data",mc_steps,particles,n,p,1,rad_choice,qpart_count,false)
hdf5_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,qpart_count,true)
only_no_overlaps_data = [hdf5_data[1],hdf5_data[2][:,1001:end],hdf5_data[3][1001:end]]
#array_data = Iterators.flatten([hdf5_data[2][i,:] for i in 2:particles])#Iterators.flatten([only_no_overlaps_data[2][i,:] for i in 2:particles])
#hist2D(real.(array_data),-imag.(array_data),bins=100)
#
dens_data = get_density(particles,n,p,only_no_overlaps_data,100)
qpart_location = hdf5_data[1][2][1]
#
plot(dens_data[1]./rm,dens_data[2],label="$rad_choice")#"$n/$fill_denom")
scatter([real(qpart_location)]./rm,[dens_data[2][1]],label="$rad_choice")
legend()
#end
=#

"fin"
