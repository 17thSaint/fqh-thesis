using HDF5,PyPlot,Statistics
#=
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
=#
include("cf-wavefunc.jl")
include("read-CF-data.jl")
function prob_wavefunc(config, m, qhole, num_parts)
	full = 0
	for j = 1:num_parts
		for i = 1:j-1
			dist = sqrt((config[1][i]-config[1][j])^2 + (config[2][i]-config[2][j])^2)
			full += -2*m*log( dist )
		end
		if qhole[1] > 0
			full += -2*log( sqrt((config[1][j]-qhole[1])^2 + (config[2][j]-qhole[2])^2))
		end
		full += 0.5 * (config[1][j]^2 + config[2][j]^2)
	end
	return full
end

function get_energies_alltimes(xs,ys,m,qhole,num_parts,dt)
	time_counts = Int(size(xs)[2]/dt)
	energies = [0.0 for i in 1:time_counts]
	for i in 1:time_counts
		#if i % (0.05*time_counts) == 0
		#	println("Getting Energies: ",(i-1)*100/time_counts,"%")
		#end
		time = Int(1+(i-1)*dt)
		input_config = 
		local_energy = prob_wavefunc(input_config,m,qhole,num_parts)
		energies[i] = local_energy
	end
	return energies
end

function auto_correlation(energies, delta_t)
    average_energy = mean(energies)
    
    points = Int(floor(length(energies)-delta_t))
    
    energy_fluctuations = energies.-average_energy
    
    autocorrelation_top = [0.0 for i in 1:points]
    autocorrelation_bottom = mean(energy_fluctuations.^2)
    
    for i in 1:points
        autocorrelation_top[i] = energy_fluctuations[i]*energy_fluctuations[i+delta_t]
    end
    
    return (mean(autocorrelation_top)/autocorrelation_bottom)
end

function get_autocorr_length(wavefunc_data,samp_freq)
	energy = 2 .*real.(wavefunc_data)
	full_length = length(wavefunc_data)
	#println("Full Length = $full_length")
	len = 1000
	dts = [1+(i-1)*1 for i in 1:Int(0.1*len)-1]
	autocorr = [0.0 for i in 1:Int(0.1*len)-1]
	for i in 1:Int(0.1*len)-1
		autocorr[i] = auto_correlation(energy,dts[i])
	end
	check_below_tol = autocorr .< [0.01 for i in 1:length(autocorr)]
	#println(autocorr)
	corr_length = samp_freq*dts[findall(check_below_tol)[1]]
	return corr_length,dts,autocorr
end

#=
particles = 16
n,p = 1,1
rad_count = 5
qpart_count = 2
log_form = true
hdf5_data = read_comb_CF_hdf5("cf-data",particles,n,p,rad_count,qpart_count,log_form)
rezz = get_autocorr_length(hdf5_data[3],1)
println(rezz[1])
=#



"fin"
