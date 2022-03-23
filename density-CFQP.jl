using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")
include("mc-cf.jl")

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
	end
	return xs,change_counts
end

#
particles = 12
mc_steps = 10000
np_vals = [[1,1],[1,2],[2,1]]
which = 1
n,p = 1,0#np_vals[which]
fill_denom = 2*n*p + 1
rm = sqrt(2*particles*fill_denom/n)
#=
x_rads = [0.01*rm + j*(1.29*rm)/10 for j in 0:9]
rad_choice = 5
x_rad = x_rads[rad_choice]
=#
qpart = [0,[1.0+im*0.0]]
#qp2_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,2,true)
#qp1_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_choice,1,true)
#
xs_plot = []
ys_plot = []
data_count = 25
starting_config = start_rand_config(particles,n,p)#1 .*qp1_data[2][:,3800]
parts_xs = real.(starting_config[2:end])
parts_ys = -imag.(starting_config[2:end])
#scatter3D(parts_xs,parts_ys,[1.0 for i in 1:particles-1],c="r")
xs = [-0.001*rm + i*(2*0.001*rm)/data_count for i in 1:data_count]
exp_prob = []
reg_prob = []
logjis = []
part1 = []
givematrix = [fill(0.0+im*0.0,(particles,particles))]
for i in 1:length(xs)
	local_x = xs[i]
	for j in 1:length(xs)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		radius = abs(starting_config[1])
		#if radius < 0.0001
		#println(radius)
		#end
		#append!(part1,[log(starting_config[1])])
		#append!(logjis,[get_logJi(starting_config,1,p)])
		#
		vals = get_wavefunc_fromlog(starting_config,n,p,qpart)
		exp_prob_local = 2*real(vals)
		#givematrix[1] = vals[2]
		append!(exp_prob,[exp_prob_local])
		#
		#=reg_prob_local = log(abs2(get_wavefunc(starting_config,n,p,qpart)))
		if isnan(reg_prob_local) | isinf(reg_prob_local)
			println("Divergent: $i $j")
		end
		append!(reg_prob,[reg_prob_local])
		=#
	end
end
#

#scatter3D(xs_plot,ys_plot,reg_prob./maximum(reg_prob),label="Reg")
scatter3D(xs_plot,ys_plot,exp_prob./maximum(exp_prob),label="Exp")
#sleep(2.0)
#scatter3D(xs_plot,ys_plot,part1+logjis,label="One")

legend()
xlabel("X")
ylabel("Y")
#






#array_data = Iterators.flatten([qp1_data[2][i,:] for i in 1:1])#Iterators.flatten([only_no_overlaps_data[2][i,:] for i in 2:particles])
#hist2D(real.(array_data),-imag.(array_data),bins=100)
#
#dens_data_1q = get_density(particles,n,p,qp1_data,200)
#dens_data_2q = get_density(particles,n,p,qp2_data,200)
#qpart_location = qp1_data[1][2][1]
#
#plot(dens_data_1q[1]./rm,dens_data_1q[2],label="1-QP")#"$n/$fill_denom")
#plot(dens_data_2q[1]./rm,dens_data_2q[2],label="2-QP")#"$n/$fill_denom")
#scatter([real(qpart_location)]./rm,[dens_data_1q[2][1]],label="QP")
#legend()
#end
#

"fin"
