include("cf-wavefunc.jl")
using PyPlot

function get_logdet_diff(log_matrix,exp_shift)
	parts = length(exp_shift)
	multip_log_matrix = fill(0.0+im*0.0,(parts,parts))
	#exp_shift = [rand(Float64)+im*rand(Float64) for j in 1:parts]
	for j in 1:parts
	multip_log_matrix[:,j] = log_matrix[:,j] .+ exp_shift[j]
	end
	logdet_multip = get_log_det(multip_log_matrix)
	#println(multip_log_matrix,", ",log_matrix)
	logdet_reg = get_log_det(log_matrix)
	expected_shift = sum(exp_shift)
	diff = logdet_multip - (logdet_reg + expected_shift)
	#println(round(imag(diff)/pi,digits=0),", ",imag(diff)/pi)
	#if !isapprox(round(imag(diff)/pi,digits=0),imag(diff)/pi,atol=10^(-1))
	#	println(imag(diff)/pi)
	#end
	return diff
end

#part_vals = [i for i in 2:50]
diff_vals = [0.0+0.0*im for i in 1:100]#length(part_vals)]
for i in 1:100#length(part_vals)
	parts = 12#part_vals[i]
	log_matrix = 10 .*(rand(parts,parts) + im*rand(parts,parts))
	log_matrix[1,:] = [0.0 for j in 1:parts]
	exp_shift = [rand(Float64)*rand((-1,1))+im*rand(Float64)*rand((-1,1)) for j in 1:parts]
	diff_local = get_logdet_diff(log_matrix,exp_shift)
	diff_vals[i] = diff_local
end

plot(real.(diff_vals),label="Re")
#plot(part_vals,imag.(diff_vals)./pi,"-p",label="Im")
#legend()
