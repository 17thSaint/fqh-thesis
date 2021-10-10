function start_rand_config(num_parts,m)
	rm = sqrt(2*num_parts*m)
	[[rand(Float64)*rand(-1:2:1)*rm,rand(Float64)*rand(-1:2:1)*rm] for i = 1:num_parts]
end

function dist_btw(part_1,part_2)
	return sqrt((part_1[1] - part_2[1])^2 + (part_1[2] - part_2[2])^2)
end

function prob_wavefunc(config, m)
	full = 0
	for j = 1:size(config)[1]
		for i = 1:j-1
			dist = dist_btw(config[i],config[j]) 
			full += 2*m*abs(log( dist ) )
		end
		full += 0.5 * (config[j][1]^2 + config[j][2]^2)
	end
	return full
end

function move_particle(num_parts,chosen,step_size)
	shift_matrix = [[0.0,0.0] for i = 1:num_parts]
	shift_matrix[chosen][1] += rand(-1:2:1)*rand(Float64)*step_size
	shift_matrix[chosen][2] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config,chosen,num_parts,m,step_size)
	start_ham = prob_wavefunc(config,m)
	shift_matrix = move_particle(num_parts,chosen,step_size)
	new_ham = prob_wavefunc(config+shift_matrix,m)
	rand_num = rand(Float64)
	if exp(-(new_ham - start_ham)) > rand_num
		return config+shift_matrix, 1, start_ham, new_ham, rand_num, config
	else
		return config, 0, start_ham, new_ham, rand_num
	end
	
	return "Acceptance Calculation Error"
end

#include("ProgressTimer.jl")

function main(steps,num_parts,m,step_size)
	running_config = start_rand_config(num_parts,m)
	samp_freq = 4
	therm_time = 0.1*steps
	time_config_x = fill(0.0,(num_parts,steps))
	time_config_y = fill(0.0,(num_parts,steps))
	#rej_rate = fill(0.0,(num_parts,Int(count_end*0.001)))
	#rej_count = [0.0 for w in 1:num_parts]
	#energies = fill(0.0,steps*num_parts)
	index = 1
	for i in 1:steps*samp_freq
		if i%(steps*0.05) == 0
			println(100*i/steps)
		end
		for j in 1:num_parts
			movement = acc_rej_move(running_config,j,num_parts,m,step_size)
			running_config = movement[1]
			if dist_btw(running_config[1],running_config[2]) < 0.01
				println(movement)
				error("Ouch, too close!")
			end
			#rej_count[j] += movement[2]/1000
			#energies[(i-1)*num_parts+j] = movement[4]
		end
		
		if i > therm_time && i%samp_freq == 0
			time_config_x[:,index] = [running_config[x][1] for x in 1:num_parts]
			time_config_y[:,index] = [running_config[x][2] for x in 1:num_parts]
			index += 1
		end
		
		#if i%1000 == 0
		#	rej_rate[:,Int(i*0.001)] = rej_count
		#end
		#if i%(steps*0.001) == 0
		#	rej_count = [0.0 for k in 1:num_parts]
		#end
		#energies[i] = prob_wavefunc(running_config,m)
	end
	
	#println(rej_count)
	return time_config_x,time_config_y#,energies
end

using PyPlot
function plotting_trajs(data)
	plot(transpose(data[1]),transpose(data[2]))	
end


using DelimitedFiles
function save_data(data,num_parts,m)
	cd("MC-Data/MC-Data/")
	writedlm("xpos-P$num_parts-M$m.csv",data[1],',')
	writedlm("ypos-P$num_parts-M$m.csv",data[2],',')
	cd("../../")
end

mc_steps = 5000
step_size = 0.1
m = 3
#rm = sqrt(2*m*particles)
particles = 4

main_data = main(mc_steps,particles,m,step_size)



#=
for j in 2:8
	if j == 4
		continue
	end
	particles = j
	for i in 1:3
		m = i
		main_data = main(mc_steps,particles,m,step_size)
		save_data(main_data,particles,m)
		println("Saved M=",m,", ","P=",particles)
	end
end
=#


#=
for i in 1:particles
	println("Particle=",i)
	plot_ref = hist2D(main_data[1][i,:],main_data[2][i,:],100)
	savefig("part-$i-M$m-position")
end
=#

#hist(mag_dist[2],100)


#k = 100
#println(main_data[1],main_data[2])
#println([[main_data[1][(k-1)*particles+l],main_data[2][(k-1)*particles+l]] for l in 1:particles])

#rm = sqrt(2*m*particles)
#plot(transpose(main_data[1]),transpose(main_data[2]),-rm:rm,[sqrt(rm^2-x^2) for x=-rm:rm],-rm:rm,[-sqrt(rm^2-x^2) for x=-rm:rm])


#=
using CurveFit, LambertW, Statistics
wait_times = [[0.0 for i in 1:20] for j in 1:8]
for j in 2:9
	println("Particles=",j)
	for i in 1:20
		main_data = main(mc_steps,j,m,step_size)
		fit = log_fit([i*1.0 for i in 1:mc_steps], main_data[3])
		wait_times[j-1][i] = 20/lambertw(20*exp(fit[1]/fit[2]))
	end
	println("Mean=",mean(wait_times[j-1]),", ","STD=",std(wait_times[j-1]))
end
y0b = [fit[1]+fit[2]*log(i) for i in 1:mc_steps]
plot([i*1.0 for i in 1:mc_steps],main_data[3],[i*1.0 for i in 1:mc_steps],y0b)

=#




"fin"
