using HDF5, PyPlot
function read_hdf5_data(num_parts,m,data_type)
	cd("..")
	cd("qhole-data")
	file = h5open("$data_type-pos-mc-2000000-p-$num_parts-m-$m.hdf5","r")
	#file = h5open("acc-rate-data.hdf5","r")
	#data = read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
	data = read(file["all-data"],"deets")
	cd("..")
	cd("Codes")
	return data
end

particles = 5
m = 3
rm = sqrt(2*m*particles)
here_data_x = read_hdf5_data(particles,m,"x")
here_data_y = read_hdf5_data(particles,m,"y")

for i in 1:particles
	hist2D(here_data_x[i,:], here_data_y[i,:],bins = 100,range = [[-rm,rm],[-rm,rm]])
	sleep(2)
end








"fin"
