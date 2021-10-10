l0 = 1


function start_config(num_part_edge)
	[convert(Vector{Float64},[-0.5*num_part_edge+i,-0.5*num_part_edge+j]) for i = 0:num_part_edge-1 for j = 0:num_part_edge-1]
end

function dist_btw(part_1,part_2)
	return sqrt((part_1[1] - part_2[1])^2 + (part_1[2] - part_2[2])^2)
end


function prob_wavefunc(config, m)
	full = 0
	for j = 1:size(config)[1]
		for i = 1:j-1
			full += -2*m*abs(log( dist_btw(config[i],config[j]) ) )
		end
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)/(l0^2)
	end
	return full
end

using PyPlot
x = [i for i = -5:5 for j = 0:10]
y = [i for j = 0:10 for i = -5:5]
scatter3D(x,y,[exp(-prob_wavefunc([[0.0,0.0],[i,j]],3)) for i=-5:5 for j=-5:5])



#=
function calc_elec_dens(x,y,config)
	elec_dens = 0
	for i = 1:size(config)[1]
		elec_dens += ==(x,config[i][1])/(2*size(config)[1]) + ==(y,config[i][2])/(2*size(config)[1])
	end
	return elec_dens
end

function rho(x,y,config,m,num_part_edge)
	rm = sqrt(2*m)*l0*num_part_edge
	if sqrt(x^2+y^2) <= rm
		return calc_elec_dens(x,y,config) - 1/(2*pi*m*(l0^2))
	else
		return calc_elec_dens(x,y,config)
	end
end



function calc_pair_corr(config,m,num_part_edge)
	rm = sqrt(2*m)*l0*num_part_edge
	inside_count = 0
	for i = 1:num_part_edge^2
		if sqrt(config[i][1]^2+config[i][2]^2) <= rm
			inside_count += 1
		end
	end
	return (exp(-class_hamilt(config,m))) * ((num_part_edge^2)/area - inside_count/(2*pi*m*(l0^2)))
end
=#


function move_particle(config)
	chosen = rand(1:size(config)[1])
	config[chosen][1] += rand(-1:2:1)*rand(Float64)*100
	config[chosen][2] += rand(-1:2:1)*rand(Float64)*100
	return config
end

function acc_rej_move(config,m)     
	start_nrg = prob_wavefunc(config,m)
	new_config = move_particle(config)
	new_nrg = prob_wavefunc(new_config,m)
	if exp(-(new_nrg - start_nrg)) > rand(Float64)
		return new_config, 1
	else
		return config, 0
	end
	
	return "Acceptance Calculation Error"
end

function main(steps,num_part_edge,m)
	start = start_config(num_part_edge)
	running_config = start
	energies = []
	stepped_config_x = fill(0.0,(num_part_edge^2,steps))
	stepped_config_y = fill(0.0,(num_part_edge^2,steps))
	rej_count = 0
	for i = 1:steps
		if i%(steps*0.01) == 0
			println(100*i/steps)
		end
		push!(energies,prob_wavefunc(running_config,m))
		for j = 1:num_part_edge^2
			stuff = acc_rej_move(running_config,m)
			running_config = stuff[1]
			rej_count += stuff[2]/(steps*(num_part_edge^2))
		end
		stepped_config_x[:,i] = [running_config[x][1] for x in 1:(num_part_edge^2)]
		stepped_config_y[:,i] = [running_config[x][2] for x in 1:(num_part_edge^2)]
	end
	println(rej_count)
	return energies,stepped_config_x,stepped_config_y
end

#=
function plot_particles(config)
	plot([config[i][1] for i = 1:size(config)[1]],[config[j][2] for j = 1:size(config)[1]], seriestype = :scatter)
end


#function plot_energy()
m = 3
steps = 500000
data_3, x_pos, y_pos = main(steps,3,m) 
#plot(transpose(x_pos),transpose(y_pos))
plot(1:steps,data_3)
#plot_particles(data_3[2])
	#plot([1:size(data_1[1])[1],1:size(data_2[1])[1],1:size(data_3[1])[1]], [data_1[1],data_2[1],data_3[1]] )
#end

=#



"fin"
