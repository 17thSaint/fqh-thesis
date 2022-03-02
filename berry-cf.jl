using HDF5,Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")


function get_expval(n,p,hdf5_data,chosen_qpart,dtheta)
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
top = 10
#=
rads = fill(0.0,(top,3))
berries = fill(0.0,(top,3))
errors = fill(0.0,(top,3))
np_vals = [[1,1],[1,2],[2,1]]
for j in 1:3
        n,p = np_vals[j]
        for i in 1:top 
            radii_data = read_CF_hdf5("cf-data",mc_steps,particles,n,p,i,1)
            rads[i,j] = real(radii_data[1][2][1])
            berry_calc = get_expval(n,p,radii_data,1,-0.001)
            berries[i,j] = berry_calc[1]
            errors[i,j] = berry_calc[2]
        end
end
=#
j = 3
n,p = np_vals[j]
fill_denom = 2*p*n + 1
rm = sqrt(2*particles*fill_denom/n)
theory_here = rads[:,j].*rads[:,j]./(1*(2*p*n+1))
plot(rads[:,j]./rm,theory_here,label="~TH")# $n / $fill_denom")
errorbar(rads[:,j]./rm,-berries[:,j],yerr=[errors[:,j],errors[:,j]],fmt="-o",label="EXP")# $n / $fill_denom")
legend()
xlabel("Quasiparticle Radius r/rm")
title(latexstring("Calculated Berry Phase for \$ \\nu = $n / $fill_denom \$"))







"fin"
