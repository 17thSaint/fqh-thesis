using Statistics,PyPlot

include("read-CF-data.jl")
include("cf-wavefunc.jl")

function get_expval(n,p,hdf5_data,chosen_qpart,dtheta,log_form=false)
	qpart_data,pos_data,wavefunc_data = hdf5_data
	len = length(wavefunc_data)
	single_qpart_location = qpart_data[2][chosen_qpart]
	qpart_shifted = single_qpart_location*exp(im*dtheta)
	qpart_data[2][chosen_qpart] = qpart_shifted
	calced_vals = [0.0 for i in 1:len]
	for i in 1:len # averaged over a bunch of time slices separated by sampling frequency
			# sampling frequency found from time autocorrelation: ~100
		local_config = pos_data[:,i]
		
		if log_form
			new_wavefunc = get_wavefunc_fromlog(local_config,n,p,qpart_data)
			ratio_wavefunc = new_wavefunc - wavefunc_data[i]
			calced_vals[i] = exp(real(ratio_wavefunc))*sin(imag(ratio_wavefunc))/dtheta
			#println(calced_vals[i])
		else
			new_wavefunc = get_wavefunc(local_config,n,p,qpart_data)
			ratio_wavefunc = new_wavefunc / wavefunc_data[i]
			calced_vals[i] = imag(ratio_wavefunc)/dtheta
		end
	end
	final_val = mean(calced_vals)
	val_error = std(calced_vals)
	return final_val,val_error
end

particles = 16
mc_steps = 1000
top = 10
#
rads_1q = fill(0.0,(top,3))
berries_1q = fill(0.0,(top,3))
errors_1q = fill(0.0,(top,3))
rads_2q = fill(0.0,(top,3))
berries_2q = fill(0.0,(top,3))
errors_2q = fill(0.0,(top,3))
np_vals = [[1,1],[1,2],[2,1]]
for j in 1:1
        n,p = np_vals[j]
        for i in 1:top 
            #=
            radii_data_1q = read_CF_hdf5("cf-data",mc_steps,particles,n,p,i,1)
            rads_1q[i,j] = real(radii_data_1q[1][2][1])
            berry_calc_1q = get_expval(n,p,radii_data_1q,1,-0.001)
            berries_1q[i,j] = berry_calc_1q[1]
            errors_1q[i,j] = berry_calc_1q[2]
            =#
            radii_data_2q = read_CF_hdf5("cf-data",mc_steps,particles,n,p,i,2,true)
            rads_2q[i,j] = real(radii_data_2q[1][2][1])
            berry_calc_2q = get_expval(n,p,radii_data_2q,1,-0.001,true)
            berries_2q[i,j] = berry_calc_2q[1]
            errors_2q[i,j] = berry_calc_2q[2]
        end
end
#


j = 1
n,p = np_vals[j]
th_coeff_vals = [2,1,1]
th_coeff = th_coeff_vals[j]
fill_denom = 2*p*n + 1
#rm = sqrt(2*particles*fill_denom/n)

# direct compare berry phase
#theory_here_1q = rads_1q[:,j].*rads_1q[:,j]./(th_coeff*(2*p*n+1))
#
theory_here_2q_m = rads_2q[:,j].*rads_2q[:,j]./(th_coeff*(2*p*n+1)) - [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,j])]
theory_here_2q_p = rads_2q[:,j].*rads_2q[:,j]./(th_coeff*(2*p*n+1)) + [2*p/(2*p*n + 1) for i in 1:length(rads_2q[:,j])]
#plot(rads_1q[:,j],theory_here_1q,label="~TH 1Q")# $n / $fill_denom")
plot(rads_2q[:,j],theory_here_2q_m,label="~TH 2Q M")# $n / $fill_denom")
plot(rads_2q[:,j],theory_here_2q_p,label="~TH 2Q P")# $n / $fill_denom")
#errorbar(rads_1q[:,j],-berries_1q[:,j],yerr=[errors_1q[:,j],errors_1q[:,j]],fmt="-o",label="1Q")# $n / $fill_denom")
errorbar(rads_2q[:,j],-berries_2q[:,j],yerr=[errors_2q[:,j],errors_2q[:,j]],fmt="-o",label="2Q")# $n / $fill_denom")
xlabel("Quasiparticle Radius")
ylabel("Berry Phase")
title(latexstring("Berry Phase for \$ \\nu = $n/$fill_denom \$"))
#

#= compare shift in berry phase with 2QP
#diff = -berries_1q[:,j]+berries_2q[:,j]
diff = -theory_here_1q[:,j]+berries_2q[:,j]
#diff_errors = sqrt.(errors_1q[:,j].^2 + errors_2q[:,j].^2)
diff_errors = sqrt.(errors_2q[:,j].^2)
theory_diff = [2*p/(2*p*n + 1),-2*p/(2*p*n + 1)]
errorbar(rads_2q[:,j],diff,yerr=[diff_errors,diff_errors],fmt="-o")
plot(rads_2q[:,j],[theory_diff[1] for i in 1:length(rads_2q[:,j])],label="TH P")
plot(rads_2q[:,j],[theory_diff[2] for i in 1:length(rads_2q[:,j])],label="TH M")
legend()
xlabel("Quasiparticle Radius")
title(latexstring("Shift from Enclosed CFQP at \$ \\nu=$n/$fill_denom \$"))
=#







"fin"
