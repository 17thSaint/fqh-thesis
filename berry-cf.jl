using HDF5,Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")


function get_expval(hdf5_data,chosen_qpart,dtheta)
	qpart_data,pos_data,wavefunc_data = hdf5_data
	len = length(wavefunc_data)
	single_qpart_location = qpart_data[2][chosen_qpart]
	qpart_shifted = single_qpart_location*exp(im*dtheta)
	qpart_data[2][chosen_qpart] = qpart_shifted
	calced_vals = [0.0 for i in 1:len]
	for i in 1:len # averaged over a bunch of time slices separated by sampling frequency
			# sampling frequency found from time autocorrelation: ~100
		local_config = pos_data[:,i]
		new_wavefunc = get_wavefunc(local_config,n,p,qpart_data)
		ratio_wavefunc = new_wavefunc / wavefunc_data[i]
		calced_vals[i] = imag(ratio_wavefunc)/dtheta
	end
	final_val = mean(calced_vals)
	val_error = std(calced_vals)
	return final_val,val_error
end

particles = 8
mc_steps = 100000
n = 1
p = 1
fill_denom = 2*n*p + 1
filling = round(n/(2*1*n+1),digits=3)
rm = sqrt(2*particles/filling)
top = 10
rads = [0.0 for i in 1:top]
berries = [0.0 for i in 1:top]
errors = [0.0 for i in 1:top]
for i in 1:top
	radii_data = read_CF_hdf5("Codes",mc_steps,particles,n,p,i,1)
	rads[i] = real(radii_data[1][2][1])
	berry_calc = get_expval(radii_data,1,-0.001)
	berries[i] = berry_calc[1]
	errors[i] = berry_calc[2]
end
theory_here = rads.*rads./(2 .*(p*2+1))
plot(rads./rm,theory_here,"-r",label="TH")
errorbar(rads./rm,-berries,yerr=[errors,errors],fmt="-ko",label="EXP")
legend()
xlabel("Quasiparticle Radius r/rm")
title(latexstring("Calculated Berry Phase for \$ \\nu = $n / $fill_denom \$"))








"fin"
