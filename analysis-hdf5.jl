using HDF5, PyPlot
function read_hdf5_data(num_parts,m,data_type)
	cd("mc-data")
	#file = h5open("$data_type-mc2000000-notherm.hdf5","r")
	file = h5open("acc-rate-data.hdf5","r")
	data = read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
	cd("..")
	return data
end


plot([hist_m3_p8[1,:] hist_m3_p7[1,:] hist_m3_p6[1,:] hist_m3_p5[1,:] hist_m3_p4[1,:] hist_m3_p3[1,:] hist_m3_p2[1,:]],[hist_m3_p8[2,:] hist_m3_p7[2,:] hist_m3_p6[2,:] hist_m3_p5[2,:] hist_m3_p4[2,:] hist_m3_p3[2,:] hist_m3_p2[2,:]],label=("8","7","6","5","4","3","2"))
legend()
xlabel("Avg Separation")
title("Latter-Half Time Histo for M=3 Range Particle #")




"fin"
