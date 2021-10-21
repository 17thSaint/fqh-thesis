using HDF5
function read_hdf5_data(num_parts,m,data_type,folder)
	cd("..")
	cd("$folder")
	#file = h5open("acc-rate-data.hdf5","r")
	if folder == "mc-data"
		file = h5open("$data_type-mc2000000-notherm.hdf5","r")
		data = read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
	elseif folder == "qhole-data"
		file = h5open("$data_type-pos-mc-2000000-p-$num_parts-m-$m-qhole-mk2.hdf5","r")
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

function get_expval(particles,full_data_x,full_data_y,dtheta,rad)
	qhole_location = full_data_x[2][1] + im*full_data_x[2][2]
	chord = 2*rad*sin(dtheta/2)*exp(im*0.5*(pi+dtheta))
	qhole_shifted = qhole_location + chord
	final_config = [ last(full_data_x[1][i,:]) + im*last(full_data_y[1][i,:]) for i in 1:particles]
	#final_config = [ full_data_x[1][i,length(full_data_x[1][i,:])-10] + im*full_data_y[1][i,length(full_data_y[1][i,:])-10] for i in 1:particles]
	exp_val = 0
	for i in 1:particles
		part = (final_config[i] - qhole_shifted)/(final_config[i] - qhole_location)
		exp_val += part/particles
	end
	a_theta = (exp_val-1)/dtheta
	return a_theta
end

using PyPlot
particles = 20
start_dtheta = 0
end_dtheta = 2*pi
radius = 0.01
count = 20
vals_dtheta = [start_dtheta+i*(end_dtheta-start_dtheta)/count for i in 0:count-1]
atheta_dtheta = [[0.0+im*0.0,0.0+im*0.0,0.0+im*0.0] for i in 1:count]
for i in 1:3
	m = i
	data_x = read_hdf5_data(particles,m,"x","qhole-data")
	data_y = read_hdf5_data(particles,m,"y","qhole-data")
	for j in 1:count
		println("M=",m,", ","Count=",j)
		thet = start_dtheta+(j-1)*(end_dtheta-start_dtheta)/count
		atheta_dtheta[j][i] = get_expval(particles,data_x,data_y,thet,radius)
	end
	plot(vals_dtheta,[imag(atheta_dtheta[j][i]) for j in 1:count],label="Im M=$i")
	plot(vals_dtheta,[real(atheta_dtheta[j][i]) for j in 1:count],label="Re M=$i")
	plot(vals_dtheta,[-(radius^2)/(4*m) for i in 1:count],label="TH M=$m")
end
legend()




"fin"