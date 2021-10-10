using DelimitedFiles
function read_data(num_parts,m)
	cd("mc-data")
	xpos = readdlm("xpos-P$num_parts-M$m.csv",',',Float64)
	ypos = readdlm("ypos-P$num_parts-M$m.csv",',',Float64)
	cd("..")
	#println("Done Reading M=",m,", ","P=",num_parts)
	return [xpos,ypos]
end

function save_hist_data(data,num_parts,m)
	cd("mc-data")
	writedlm("histo-P$num_parts-M$m.csv",data,',')
	cd("..")
end

function get_complex_pos(particles,mc_steps,m,data)
	z = fill(0.0+im*0.0,(particles,mc_steps*particles))
	for i in 1:particles
		z[i,:] = data[1][i,:]+im*data[2][i,:]		
	end
	return z
end

function get_mag_dist_btw(particles,mc_steps,m,data)
	choices = factorial(particles)/(2*factorial(particles-2))
	mag_dist_btw = fill(0.0,(Int(choices),particles*mc_steps))
	z = get_complex_pos(particles,mc_steps,m,data)
	why = 0
	for i in 1:particles
		for j in 1:i-1
			why += 1
			mag_dist_btw[why,:] = sqrt.(abs.(z[i,:]-z[j,:]))
		end
	end
	avg_mag_dist = fill(0.0,particles*mc_steps)
	for i in 1:Int(choices)
		avg_mag_dist += mag_dist_btw[i,:]./choices
	end
	#println("Done Mag Dist M=",m,", ","P=",particles)
	return mag_dist_btw, avg_mag_dist, choices
end

using PyPlot
function hist_sep(data,choices,bins,m,particles)  # makes histogram of average particle separation
	rm = sqrt(2*m*particles)		    # for second half of simulation with max of rm
	bin_count = fill(0,2*bins)
	x_bins = fill(0.0,2*bins-2)
	for j in 1:Int(2*bins-2)
		if j%2 == 0
			x_bins[j] = 0.5*j*rm/bins
		else
			x_bins[j] = 0.5*(j+1)*rm/bins
		end
	end
	x_bins = append!([0.0],x_bins)
	append!(x_bins,[rm])
	for i in 1:choices
		for j in 1:bins
			for k in Int(mc_steps/2):mc_steps
				if data[i,k] < j*rm/bins && data[i,k] >= (j-1)*rm/bins
					bin_count[2*j-1] += 1
					bin_count[2*j] += 1
				end
			end
		end
	end
	#println("Done Histo M=",m,", ","P=",particles)
	return x_bins,bin_count
end


bins = 100
mc_steps = 2000000

#main_data = [read_data(j,i) for i in 1:3 for j in 2:8]
#mag_dist = [get_mag_dist_btw(j,mc_steps,i,main_data[i]) for i in 1:3 for j in 2:8]
#hist_data = [hist_sep(mag_dist[i][1],Int(mag_dist[i][3]),bins,i,j) for i in 1:3 for j in 2:8]

function make_all(particles,m)
	main_data = read_data(particles,m)
	println("Done Reading M=",m,", ","P=",particles)
	mag_dist = get_mag_dist_btw(particles,2000000,m,main_data)
	println("Done Mag Dist M=",m,", ","P=",particles)
	hist_data = hist_sep(mag_dist[1],Int(mag_dist[3]),100,m,particles)
	println("Done Histo M=",m,", ","P=",particles)
	save_hist_data(hist_data,particles,m)
	println("Done Saving M=",m,", ","P=",particles)
end

for i in 1:3
	for j in 2:8
		make_all(j,i)
	end
end

#=
for i in 1:3
	plot(hist_data[i][1],hist_data[i][2],label="$i")
	legend()
end
=#

# LLL m=1 means all electrons are in same energy level
# B/C Pauli exclusion they need to be more spacially separated
# which is why histogram shows more distribution of particle separation
# makes sense that peak is at magnetic length l0=1




"fin"
