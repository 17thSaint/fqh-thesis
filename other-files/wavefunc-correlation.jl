using HDF5,Statistics,PyPlot

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

include("read-CF-data.jl")
particles = 8
mc_steps = 100000
step_size = 0.5
n = 2
p = 1
qpart = [1,[0.0+im*0.0]]
full_dats = read_CF_hdf5("Codes",mc_steps,particles,n,p,1,qpart)
len = length(full_dats[end])
upper = 100
dts = [i for i in 1:upper]
autocorr = [0.0 for i in 1:upper]
for i in 1:upper
	#if i % 0.05*(Int(0.1*len)-1) == 0
	#	println(100*i/(Int(0.1*len)-1))
	#end
	wavefunc_data = abs2.(full_dats[end])
	autocorr[i] = auto_correlation(wavefunc_data,dts[i])
end
plot(dts,autocorr)








"fin"
