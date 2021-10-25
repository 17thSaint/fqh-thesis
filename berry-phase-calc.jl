using HDF5,PyPlot,CurveFit
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
	#distances = []
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
	return imag(exp_val)
end

particles = 20
steps = 1000
xdats = [ [read_hdf5_data(particles,i,"x","qhole-data","mk2"),read_hdf5_data(particles,i,"x","qhole-data","rnd"), read_hdf5_data(particles,i,"x","qhole-data","og")] for i in 1:3]
ydats = [ [read_hdf5_data(particles,i,"y","qhole-data","mk2"),read_hdf5_data(particles,i,"y","qhole-data","rnd"), read_hdf5_data(particles,i,"y","qhole-data","og")] for i in 1:3]

calc_vals = [[0.0+im*0.0,0.0+im*0.0,0.0+im*0.0] for i in 1:10]
rad = 0.001
thetas = [0.0 for i in 1:10]
th_vals = [-(rad^2)/(4*i) for i in 1:3]
for i in 1:10
	thet = 0.001 + (i-1)*(0.1/20)
	thetas[i] = thet
	for j in 1:3
		calc_vals[i][j] = get_expval(particles,xdats[3][j],ydats[3][j],thet,rad,steps)
	end
end
plot(thetas,transpose(calc_vals))

# take time autocorrelation function from ModSim to see if 1000 is good, might be the issue






"fin"
