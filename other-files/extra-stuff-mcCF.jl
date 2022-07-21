#write_pos_data_hdf5(mc_steps,particles,n,step_size,qpart,rezz,1)
#=
plot(real(transpose(rezz)),imag(transpose(rezz)))
for i in 1:length(qpart[2])
	scatter([real(qpart[2][i])],[imag(qpart[2][i])])
end
xs = collect(Iterators.flatten([real(rezz[i,:]) for i in 1:particles]))
ys = collect(Iterators.flatten([-imag(rezz[i,:]) for i in 1:particles]))
hist2D(xs,ys,bins=100)
title(latexstring("Histogram with Quasihole and \$ \\nu = $filling \$"))

radii = collect(Iterators.flatten([abs2.(rezz[i,:]) for i in 1:particles]))
hist(radii,bins=100)
title_string = ["without Quasiparticle","with Quasihole at Origin"][qpart[1]+1]
title(latexstring("Radial Distribution $title_string, \$ \\nu = $n / $fill_denom \$"))
=#
#=
starting_config = start_rand_config(particles,n)
data_count = 20
xs = [-1.2*sqrt(2*particles*5/2) + i*(2*1.2*sqrt(2*particles*5/2))/data_count for i in 0:data_count]
xs_plot = []
ys_plot = []
probs = []
for i in 1:length(xs)
	local_x = xs[i]
	for j in 1:length(xs)
		#println(i,", ",j)
		local_y = xs[j]
		append!(xs_plot,[local_x])
		append!(ys_plot,[local_y])
		starting_config[1] = local_x - im*local_y
		prob_local = abs2(get_wavefunc(starting_config,n,qpart)[1])
		append!(probs,[prob_local])
	end
end

scatter3D(xs_plot,ys_plot,probs)
=#




# mostly RFA stuff

#=
corr_length_1q = 1#get_autocorr_length(radii_data_1q[3],1)[1]
#rads_1q[k] = x_rad
berry_calc_1q = get_expval(flux_type,n,p,[qpart,rezz[2],rezz[3]],1,-0.001,corr_length_1q,log_form,allowed_sets_matrix,full_pasc_tri,full_derivs)
berries_1q[i][k] = -berry_calc_1q[1]
errors_1q[i][k] = berry_calc_1q[2]
end
#println(berries_1q/rads_1q^2,", ",errors_1q/rads_1q^2)
#end
#

#berries_1q_6000 = [1.0459178640024025, 1.1423713786935825, 1.2809715371589483, 1.5945364110032794, 2.0085205380906896, 2.8916456297541804, 3.8306735264144267, 4.841223214894849, 5.6464704765704, 6.2129571005747675]
#errors_1q_6000 = [0.28946957606106744, 0.39294989129750857, 0.5859301762328589, 0.8841494089487836, 1.3809291752042423, 1.6818686320896157, 1.8717397644123184, 1.7927787436346159, 1.4859387283793097, 1.1089982366003532]
#berries_1q_500 = [1.0244685971603706, 1.1615802648283324, 1.3680703283025386, 1.5271108092684424, 1.8789722879208561, 3.0835231336995275, 3.678369921518863, 4.6646025897091805, 5.731212552482283, 6.017636449898667]
#errors_1q_500 = [0.28585543410645375, 0.40606235727980555, 0.5719843400401192, 0.7821428039698142, 1.393477330611598, 1.6302203676280527, 2.243530006359698, 1.9169303167864815, 1.3202948431510777, 1.404690487409908]
#errorbar(x_rads./rm,berries_1q_500./(x_rads.^2),yerr=[errors_1q_500./(x_rads.^2),errors_1q_500./(x_rads.^2)],fmt="-o",label="500")
#errorbar(x_rads./rm .+ 0.01,berries_1q_6000./(x_rads.^2),yerr=[errors_1q_6000./(x_rads.^2),errors_1q_6000./(x_rads.^2)],fmt="-o",label="6000")
#berries_1q_10000cf_n6 = [0.4430958911515426, 0.5637476101008279, 0.7010102434088423, 0.9087386374433143, 1.1540642910795673, 1.5417194017512101, 2.043677254994068, 2.7214967129768834, 3.590417674697592, 4.412568389462677]
#errors_1q_10000cf_n6 = [0.43275606894659213, 0.5203035750462965, 0.5766207984522211, 0.7008642034130591, 0.8292270869505539, 1.0629264766566975, 1.3058581693682856, 1.527302028559309, 1.5604472695881084, 1.5027006110006698]
#berries_1q_1000cf_n10 = [0.9471940168573781, 1.1448401610106105, 1.4554199856122059, 1.80600807791446, 2.230924874734249, 2.6271163621854763, 3.133084318538919, 4.101642949420335, 5.156529259928382, 6.7291089889530165]
#errors_1q_1000cf_n10 = [0.6787765860464676, 0.7503441702282811, 0.7984004526163007, 0.90672826605486, 0.895155714486338, 0.9111859145773838, 1.4334602648046986, 1.6314911785554747, 2.132204560833004, 2.1241187406734867]
#berries_1q_1000cf_n18 = [1.6159801233191649, 2.271501402769542, 2.8892093741168514, 3.613579921976991, 4.625803078730288, 5.446384288827975, 6.390984752478777, 7.183930105809739, 8.371508650032807, 10.851608477237045]
#errors_1q_1000cf_n18 = [0.786953658219406, 1.032310047365125, 1.19529589871609, 1.273033587521588, 1.3600759024408857, 1.4932618053231148, 1.708227704707404, 1.8540649803637126, 2.1220912087909314, 3.1320408270013615]
#errorbar(x_rads6./rm6,berries_1q_10000cf_n6./(x_rads6.^2),yerr=[errors_1q_10000cf_n6./(x_rads6.^2),errors_1q_10000cf_n6./(x_rads6.^2)],fmt="-o",label="N6")
#errorbar(x_rads10./rm10 .+ 0.001,berries_1q_1000cf_n10./(x_rads10.^2),yerr=[errors_1q_1000cf_n10./(x_rads10.^2),errors_1q_1000cf_n10./(x_rads10.^2)],fmt="-o",label="N10")
#errorbar(x_rads./rm .- 0.001,berries_1q_1000cf_n18./(x_rads.^2),yerr=[errors_1q_1000cf_n18./(x_rads.^2),errors_1q_1000cf_n18./(x_rads.^2)],fmt="-o",label="N18")

errorbar(x_rads,berries_1q[i]./(x_rads.^2),yerr=[errors_1q[i]./(x_rads.^2),errors_1q[i]./(x_rads.^2)],fmt="-o",label="N$particles")
end
plot([0.3,0.75],[1/4,1/4],label="1/4")
plot([0.3,0.75],[1/3,1/3],label="1/3")
legend()
xlabel("QP Radius/rm (1/3)")
title("Fractional Charge vs QP Radius (MC=500+100Therm)")
#
function get_xrad(particles,filling,which)
	rm = sqrt(2*particles/filling)
	x_rads = [0.3*rm + j*0.5*rm/10 for j in 0:9]
	return x_rads[which]
end

if false
for i in 5:length(parts)
particles = parts[i]
x_rads = [get_xrad(particles,1/3,j) for j in 1:10]
rm = sqrt(2*particles*3)
errorbar(x_rads,berries_1q[i]./(x_rads.^2),yerr=[errors_1q[i]./(x_rads.^2),errors_1q[i]./(x_rads.^2)],fmt="-o",label="N$particles")
end
#plot([0.3,0.75],[1/4,1/4],label="1/4")
#plot([0.3,0.75],[1/3,1/3],label="1/3")
legend()
xlabel("QP Radius/rm (1/3)")
title("Fractional Charge vs QP Radius (MC=500+100Therm)")
end

if false
figure()
for j in 1:10
	bers = [berries_1q[i][j]/get_xrad(parts[i],1/3,j)^2 for i in 1:length(parts)]
	ers = [errors_1q[i][j]/get_xrad(parts[i],1/3,j)^2 for i in 1:length(parts)]
	errorbar(parts,bers,yerr=[ers,ers],fmt="-ko")
end
plot([parts[1],parts[end]],[1/4,1/4],label="1/4")
plot([parts[1],parts[end]],[1/3,1/3],label="1/3")
legend()
xlabel("Number of Particles")
title("Fractional Charge vs Ne (range Radii)")
end

=#

#println(berries_1q,", ",errors_1q)

#=comb_dats = read_comb_CF_hdf5("rfa-data","RFA",particles,n,p,1,qpart[1],true)
compl = collect(Iterators.flatten([rezz[2][i,:] for i in 1:particles]))
figure()
hist2D(real.(compl),-imag.(compl),bins=100)

=#
#
#=
top = 0.5*rm
selected_dats = []
for i in 1:length(rezz[3])
	for j in 1:particles
		rad = abs(rezz[2][j,i])
		if rad <= top
			append!(selected_dats,[rezz[2][j,i]])
		end
	end
end
=#


#
#dub_qpart = [1,[qp2]]#qpart_choices[3]
#no_qpart = qpart_choices[1]
#sing_qpart = qpart_choices[2]
#qpl = qpart[2][1]


#=
flat_data = Iterators.flatten([rezz[2][i,:] for i in 1:particles])
figure()
hist2D(real.(flat_data),-imag.(flat_data),bins=100,range=[[0.5*rm,0.75*rm],[-0.25*rm,0.25*rm]])
title("$qpl")
=#

#=
tots = 10
top = 10.5
rm_1 = sqrt(2*particles*3)
#
#given_locs = [0.0+im*0.0 for i in 0:tots]
#sim_locs = [0.0+im*0.0 for i in 0:tots]
#for k in 1:tots+1
#println(k/(tots+1))
#this_loc = (top - (k-1)*top*2/tots) + 0.0*im
#given_locs[k] = this_loc
#sing_qpart = [1,[(0.5+im*0.0)*rm_1]]
nqp_qpart = [0,[0.0]]
#dub_qpart = [2,[(0.5+im*0.0)*rm_1,0.0+im*0.0]]
n,p = 1,1
#
data_count = 100
allowed_sets_matrix = get_full_acc_matrix(particles)
full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
full_derivs = get_deriv_orders_matrix(particles)
#xs = [real(qpl)-0.2*rm + i*(2*0.2*rm)/data_count for i in 1:data_count]
coords_config = [0.0+0.01*im]#start_rand_config(particles,n,p)#10*(rand(particles) + im*rand(particles))
for i in 0:particles - 2
	append!(coords_config,[1.0*rm_1*exp(im*i*1*pi/(particles-1))])
end
if true
xs = [-(top+0.05)*rm_1 + 2*(top+0.05)*rm_1*i/data_count + 0.001*rm_1 for i in 1:data_count]
#ys = [-0.25*rm + 2*0.25*rm*i/data_count + 0.001*rm for i in 1:data_count]
rf_wavefunc_nqp = fill(0.0+im*0.0,(data_count,data_count))
#rf_wavefunc_sqp = fill(0.0+im*0.0,(data_count,data_count))
#rf_wavefunc_dqp = fill(0.0+im*0.0,(data_count,data_count))
cf_wavefunc_nqp = fill(0.0+im*0.0,(data_count,data_count))
#cf_wavefunc_sqp = fill(0.0+im*0.0,(data_count,data_count))
#cf_wavefunc_dqp = fill(0.0+im*0.0,(data_count,data_count))
for i in 1:length(xs)
	#println(i/length(xs))
	local_x = xs[i]
	for j in 1:length(xs)
		local_y = xs[j]
		coords_config[1] = local_x - im*local_y
			
		#cf_wavefunc_sqp[j,i] = get_wavefunc_fromlog(coords_config,n,p,sing_qpart)
		#cf_wavefunc_dqp[j,i] = get_wavefunc_fromlog(coords_config,n,p,dub_qpart)
		#cf_wavefunc_nqp[j,i] = get_wavefunc_fromlog(coords_config,n,p,nqp_qpart)
		#rf_wavefunc_sqp[j,i] = get_rf_wavefunc(coords_config,allowed_sets_matrix,full_pasc_tri,full_derivs,sing_qpart,true)
		#rf_wavefunc_dqp[i,j]  = get_rf_wavefunc(coords_config,allowed_sets_matrix,full_pasc_tri,full_derivs,dub_qpart,true)
		rf_wavefunc_nqp[j,i]  = get_rf_wavefunc(coords_config,allowed_sets_matrix,full_pasc_tri,full_derivs,nqp_qpart,true)
	end
end
end
function get_flux_plot(wfn,coords_config,title_string)
	figure()
	#subplot(1,2,1)
	#imshow(2*real.(wfn)./maximum(2*real.(wfn)))
	#title(title_string)
	#colorbar()
	#subplot(1,2,2)
	imshow(mod.(imag.(wfn),2*pi),cmap="bwr",extent=[-top,top,-top,top].*rm_1)
	title(title_string)
	scatter(real.(coords_config[2:end]),imag.(coords_config[2:end]))
end

#diff = 2 .*real.(cf_wavefunc_nqp) - 2 .*real.(rf_wavefunc_nqp)


if true
#get_flux_plot(rf_wavefunc_sqp,coords_config,"RF QP")
get_flux_plot(rf_wavefunc_nqp,coords_config,"RF No QP")
#get_flux_plot(rf_wavefunc_dqp,coords_config,"RF 2 QP")
end
if false
#get_flux_plot(cf_wavefunc_sqp,coords_config,"CF QP")
get_flux_plot(cf_wavefunc_nqp,coords_config,"CF No QP")
#get_flux_plot(cf_wavefunc_dqp,coords_config,"CF 2 QP")
end


#=
figure()
imshow(rf_wavefunc_nqp./maximum(rf_wavefunc_nqp))
title("RF QP No")
colorbar()
figure()
imshow(cf_wavefunc_sqp./maximum(cf_wavefunc_sqp))
title("CF QP")
colorbar()
figure()
imshow(cf_wavefunc_nqp./maximum(cf_wavefunc_nqp))
title("CF QP No")
colorbar()
#figure()
#imshow(rf_wavefunc_dqp./maximum(rf_wavefunc_dqp))
#title("RF QP Conj")
#colorbar()
=#

#=
matrix_elements = findfirst(rf_wavefunc_sqp.==minimum(rf_wavefunc_sqp))
sim_loc = (xs[matrix_elements[2]] - im*xs[matrix_elements[1]])/rm
sim_locs[k] = sim_loc

end
#
figure()
plot(real.(given_locs),imag(given_locs),label="Given")
plot(real.(sim_locs),imag(sim_locs),"-p",label="Sim")
legend()

figure()
plot(abs.(sim_locs - given_locs))
title("Radius Change")
=#

#=
if false
figure()
imshow(rf_wavefunc_nqp./maximum(rf_wavefunc_nqp))
title("RF No QP")
colorbar()
end
if false
figure()
imshow(rf_wavefunc_sqp./maximum(rf_wavefunc_sqp))
title("RF QP")
colorbar()
end
if false
figure()
imshow(rf_wavefunc_dqp./maximum(rf_wavefunc_dqp))
title("RF 2QP")
colorbar()
end
if false
figure()
imshow(rf_wavefunc_sqp./maximum(rf_wavefunc_sqp)-rf_wavefunc_nqp./maximum(rf_wavefunc_nqp))
title("RF Diff 0-1")
colorbar()
end
if false
figure()
imshow(rf_wavefunc_dqp./maximum(rf_wavefunc_dqp)-rf_wavefunc_sqp./maximum(rf_wavefunc_sqp))
title("RF Diff 1-2")
colorbar()
end
if true
figure()
imshow(cf_wavefunc_nqp./maximum(cf_wavefunc_nqp))
title("CF No QP")
colorbar()
end
if true
figure()
imshow(cf_wavefunc_sqp./maximum(cf_wavefunc_sqp))
title("CF QP")
colorbar()
end
if true
figure()
imshow(cf_wavefunc_sqp./maximum(cf_wavefunc_sqp)-cf_wavefunc_nqp./maximum(cf_wavefunc_nqp))
title("CF Diff 0-1")
colorbar()
end
if false
figure()
imshow(cf_wavefunc_dqp./maximum(cf_wavefunc_dqp))
title("CF 2 QP")
colorbar()
end
if false
figure()
imshow(cf_wavefunc_sqp./maximum(cf_wavefunc_sqp)-cf_wavefunc_dqp./maximum(cf_wavefunc_dqp))
title("CF Diff 1-2")
colorbar()
end
=#

=#

#berry_n10_1q = [1.9708481669564184, 2.8380103721445114,3.7483697121119963, 5.881828355572036, 8.0668164423497, 9.113358439422495, 10.172428196169864, 10.728654041915975, 10.860118261537325,11.137260821449043]
#errors_n10_1q = [0.973671327030824, 1.244224454689576, 1.5169039542783764, 2.494056560690706, 2.650817158558882, 2.2579962916267493, 1.4529633155250201, 1.1727768098615257, 1.017772891241818, 0.7979520454548175]
