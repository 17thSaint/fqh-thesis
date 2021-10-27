using HDF5,PyPlot,CurveFit,Statistics
function read_hdf5_data(num_parts,m,data_type,folder,version)
	cd("..")
	cd("$folder")
	#file = h5open("acc-rate-data.hdf5","r")
	if folder == "mc-data"
		file = h5open("$data_type-mc2000000-notherm.hdf5","r")
		data = read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
	elseif folder == "qhole-data"
		file = h5open("$data_type-pos-mc-2000000-p-$num_parts-m-$m-qhole-$version.hdf5","r")
		qhole_data = read(file["metadata"],"qhole_position")
		data = [read(file["all-data"],"deets"),qhole_data]
	else
		println("Problem with Folder Name")
	end
	cd("..")
	cd("Codes")
	return data
end

function dist_btw(part_1,part_2)
	return sqrt((part_1[1] - part_2[1])^2 + (part_1[2] - part_2[2])^2)
end

function get_expval(particles,full_data_x,full_data_y,dtheta,rad,step)
	qhole_location = full_data_x[2][1] + im*full_data_x[2][2]
	chord = 2*rad*sin(dtheta/2)*exp(im*0.5*(pi+dtheta))
	qhole_shifted = qhole_location + chord
	sliced_data = [full_data_x[1][:,Int(i*step)] + im.*full_data_y[1][:,Int(i*step)] for i in 1:Int(1800000/step)]
	exp_val = 0
	distances = []
	for j in 1:length(sliced_data)
		exp_val_step = 1
		for i in 1:particles
			#dist = abs(sliced_data[j][i] - qhole_location)
			#append!(distances,[dist])
			#if dist < 0.002
			#	println("Time=",j,", Particle=",i,", Dist=",dist)	
			#end
			
			part = (sliced_data[j][i] - qhole_shifted)/(sliced_data[j][i] - qhole_location)
			exp_val_step *= part
		end
		exp_val += exp_val_step/length(sliced_data)
	end
	return imag(exp_val),real(exp_val)
end

function get_linear_fit(xdata,ydata)
	curve = linear_fit(xdata,ydata)
	return curve[1],curve[2]
end

function get_calc_berry(rad,xdata,ydata,steps,particles,m_final)
	calc_vals = [ [[0.0 for k in 1:10] for j in 1:length(xdata[1])] for i in 1:m_final ]
	thetas = [0.0 for i in 1:10]
	for k in 1:m_final
		for i in 1:10
			thet = -0.11 + (i-1)*0.2/10
			thetas[i] = thet
			for j in 1:length(xdats[1])
				rezz = get_expval(particles,xdata[k][j],ydata[k][j],thet,rad,steps)
				calc_vals[k][j][i] = rezz[1]
			end
		end
	end
	return calc_vals,thetas
end

function get_phase_from_ideal_dtheta(xdata,ydata,berry_vals,thetas,m_final)
	choices = Int(0.5*factorial(length(xdata[1]))/factorial(length(xdata[1])-2))
	ideal_thets = []
	params_a = [[0.0 for i in 1:length(xdata[1])] for j in 1:3]
	params_b = [[0.0 for i in 1:length(xdata[1])] for j in 1:3]
	for l in 1:m_final
		ideal = 0
		for i in 1:length(xdata)
			param = get_linear_fit(thetas,berry_vals[l][i])
			params_a[l][i] = param[1]
			params_b[l][i] = param[2]
			for j in 1:i-1
				ideal += (params_a[l][i]-params_a[l][j])/(params_b[l][j]-params_b[l][i])/choices
			end
		end
		append!(ideal_thets,[ideal])
	end
	
	phase = params_a + params_b.*ideal_thets
	return phase,ideal_thets,params_a,params_b
end

particles = 20
steps = 1000
rad = 0.001
xdats = [ [read_hdf5_data(particles,i,"x","qhole-data","mk2"),read_hdf5_data(particles,i,"x","qhole-data","rnd"), read_hdf5_data(particles,i,"x","qhole-data","og")] for i in 1:3]
ydats = [ [read_hdf5_data(particles,i,"y","qhole-data","mk2"),read_hdf5_data(particles,i,"y","qhole-data","rnd"), read_hdf5_data(particles,i,"y","qhole-data","og")] for i in 1:3]

th_vals = [-(rad^2)/(4*i) for i in 1:3]
berries = get_calc_berry(rad,xdats,ydats,steps,particles,3)
phase_data = get_phase_from_ideal_dtheta(xdats,ydats,berries[1],berries[2],3)
plot(1:3,phase_data[1])
plot(1:3,th_vals,label="Theory")
legend()






# no distances less than 0.1 with 1000 steps for all data sets and all m's
# take time autocorrelation function from ModSim to see if 1000 is good, might be the issue






"fin"
