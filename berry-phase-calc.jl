using HDF5,PyPlot,Statistics
#using CurveFit
#using Profile

function read_hdf5_data(num_parts,m,data_type,folder,version,qcount=1,long=false,q2loc=1)
	cd("..")
	cd("$folder")
	#file = h5open("acc-rate-data.hdf5","r")
	if folder == "mc-data"
		file = h5open("$data_type-mc2000000-notherm.hdf5","r")
		data = read(file["all-data"]["m-$m"]["parts-$num_parts"],"data")
	elseif folder == "qhole-data"
		if qcount <= 1
			file = h5open("$data_type-pos-mc-2000000-p-$num_parts-m-$m-qhole-$version.hdf5","r")
			qhole_data = [read(file["metadata"],"qhole_position")]
			full_poses = read(file["all-data"],"deets")
			if m == 3
				data = [full_poses[:,1:Int(size(full_poses)[2]/100)],qhole_data]
			else
				data = [full_poses,qhole_data]
			end
		elseif !long
			file = h5open("$data_type-pos-mc-2000000-p-$num_parts-m-$m-qhole-$qcount-q2loc-$q2loc-$version.hdf5","r")
			qhole_data = [read(file["metadata"],"qhole_position_$i") for i in 1:qcount]
			full_poses = read(file["all-data"],"deets")
			data = [full_poses,qhole_data]
		else
			file = h5open("$data_type-pos-p-$num_parts-m-$m-qhole-long-$qcount-q2loc-$q2loc-$version.hdf5","r")
			qhole_data = [read(file["metadata"],"qhole_position_$i") for i in 1:qcount]
			full_poses = read(file["all-data"],"deets")
			data = [full_poses,qhole_data]
		end
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

function get_expval(particles,full_data_x,full_data_y,dtheta,step) #for OG include rad
	qhole_location = full_data_x[2][1][1] + im*full_data_x[2][1][2]
	#chord = 2*rad*sin(dtheta/2)*exp(im*0.5*(pi+dtheta))
	#qhole_shifted = qhole_location + chord
	qhole_shifted = qhole_location*exp(im*dtheta)
	sliced_data = [full_data_x[1][:,Int(i*step)] + im.*full_data_y[1][:,Int(i*step)] for i in 1:Int(size(full_data_x[1])[2]/step)]
	exp_val = 0
	distances = []
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
	change_qhole = abs(qhole_location - qhole_shifted)
	return imag(exp_val)/dtheta
end

function get_expval_stats(particles,full_data_x,full_data_y,dtheta,step,length_step)
	qhole_location = full_data_x[2][1][1] + im*full_data_x[2][1][2]
	qhole_shifted = qhole_location*exp(im*dtheta)
	start = Int((step-1)*length_step + 1)
	ending = Int(step*length_step)
    	#println(start,", ",ending)
	sliced_data = [full_data_x[1][:,i] + im.*full_data_y[1][:,i] for i in start:ending]
	exp_val = 0
	distances = []
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
	change_qhole = abs(qhole_location - qhole_shifted)
	return imag(exp_val)/dtheta
end

#=
function get_linear_fit(xdata,ydata)
	curve = linear_fit(xdata,ydata)
	return curve[1],curve[2]
end
=#

function get_calc_almost_berry(xdata,ydata,steps,particles,m_final,stats=false,num_steps=100,theta_start=-0.01,theta_count=1)
	calc_vals = [ [[0.0 for k in 1:theta_count] for j in 1:length(xdata[1])] for i in 1:m_final ]
	std_vals = [ [[0.0 for k in 1:theta_count] for j in 1:length(xdata[1])] for i in 1:m_final ]
	theta_change=abs(2*theta_start)
	thetas = [0.0 for i in 1:theta_count]
	for k in 1:m_final
		for i in 1:theta_count
			thet = theta_start + (i-1)*theta_change/theta_count
			thetas[i] = thet
			
			for j in 1:length(xdata[1])
				println("Calc Berry: K=$k, J=$j")
				if !stats
					rezz = get_expval(particles,xdata[k][j],ydata[k][j],thet,steps)[1]
				else
					length_steps = Int(size(xdata[1][1][1])[2]/num_steps)
					rezz_s = [0.0 for i in 1:num_steps]
					for l in 1:num_steps
						rezz_s[l] = get_expval_stats(particles,xdata[k][j],ydata[k][j],thet,l,length_steps)[1]
					end
					rezz = mean(rezz_s)
					std_vals[k][j][i] = std(rezz_s)
				end
				calc_vals[k][j][i] = rezz
			end
		end
	end
	println("Finished Berry Calc")
	return calc_vals,std_vals
end

function get_phase_from_ideal_dtheta(xdata,ydata,berry_vals,thetas,m_final)
	choices = Int(0.5*factorial(length(xdata[1]))/factorial(length(xdata[1])-2))
	ideal_thets = []
	params_a = [[0.0 for i in 1:length(xdata[1])] for j in 1:3]
	params_b = [[0.0 for i in 1:length(xdata[1])] for j in 1:3]
	for l in 1:m_final
		ideal = 0
		for i in 1:length(xdata[1])
			println("Calc Phase, Set $i")
			param = get_linear_fit(thetas,berry_vals[l][i])
			params_a[l][i] = param[1]
			params_b[l][i] = param[2]
			for j in 1:i-1
				local_ideal = (params_a[l][i]-params_a[l][j])/(params_b[l][j]-params_b[l][i])
				println("I=$i, J=$j, Ideal=$local_ideal")
				ideal += local_ideal/choices
			end
		end
		append!(ideal_thets,[ideal])
	end
	
	phase = params_a + params_b.*ideal_thets
	return phase,ideal_thets,params_a,params_b
end


particles = 20
steps = 1
q_rad_count = 9


# single quasihole
#xdats_1q = [[read_hdf5_data(particles,j,"x","qhole-data",i) for i in 2:q_rad_count+1] for j in 3:2:7]
#ydats_1q = [[read_hdf5_data(particles,j,"y","qhole-data",i) for i in 2:q_rad_count+1] for j in 3:2:7]

#=
j = 2
comp_x_1q = collect(Iterators.flatten([xdats_1q[1][j][1][i,:] for i in 1:particles]))
comp_y_1q = collect(Iterators.flatten([ydats_1q[1][j][1][i,:] for i in 1:particles]))
figure()
hist2D(comp_x_1q,-comp_y_1q,bins=40)
colorbar()
=#
#=
berries_1q = get_calc_almost_berry(xdats_1q,ydats_1q,steps,particles,3,true)
berry_phase_1q = [ [berries_1q[1][i][j][1] for j in 1:q_rad_count] for i in 1:3]
#berry_phase_1q_errors = [ [berries_1q[2][i][j][1] for j in 1:q_rad_count] for i in 1:3]
radii_1q = [ [xdats_1q[i][j][2][1][1] for j in 1:q_rad_count] for i in 1:3]
theory_1q = [radii_1q[i].*radii_1q[i]./(2 .*(i*2+1)) for i in 1:3]


for i in 1:3
	denom_filling = 2*i+1
	plot(radii_1q[i],berry_phase_1q[i],"-p",label="1/$denom_filling")
	if i == 3
		plot(radii_1q[i],theory_1q[i],"-k",label="TH")
	else
		plot(radii_1q[i],theory_1q[i],"-k")
	end
end
=#
#plot(radii_1q[2],berry_phase_1q[2],"-p")
#plot(radii_1q[2],theory_1q[2])

#
#= two quasiholes, radius to 1.5*rm
xdats_2q = [[read_hdf5_data(particles,j,"x","qhole-data",i,2,false,1) for i in 1:q_rad_count] for j in 3:2:7]
ydats_2q = [[read_hdf5_data(particles,j,"y","qhole-data",i,2,false,1) for i in 1:q_rad_count] for j in 3:2:7]
berries_2q = get_calc_almost_berry(xdats_2q,ydats_2q,steps,particles,3,true)
berry_phase_origin = [ [berries_2q[1][i][j][1] for j in 1:q_rad_count] for i in 1:3]
berry_phase_origin_errors = [ [berries_2q[2][i][j][1] for j in 1:q_rad_count] for i in 1:3]
radii_origin = [ [xdats_2q[i][j][2][1][1] for j in 1:q_rad_count] for i in 1:3]
theory_origin = [radii_origin[i].*radii_origin[i]./(2 .*(i*2+1)) for i in 1:3]
=#

#= two quasiholes, more data, radius out to 1*rm
xdats_2q_long = [[read_hdf5_data(particles,j,"x","qhole-data",i,2,true,1) for i in 1:q_rad_count] for j in 3:2:7]
ydats_2q_long = [[read_hdf5_data(particles,j,"y","qhole-data",i,2,true,1) for i in 1:q_rad_count] for j in 3:2:7]
berries_2q_long = get_calc_almost_berry(xdats_2q_long,ydats_2q_long,steps,particles,3,true)
berry_phase_origin_long = [ [berries_2q_long[1][i][j][1] for j in 1:q_rad_count] for i in 1:3]
berry_phase_origin_long_errors = [ [berries_2q_long[2][i][j][1] for j in 1:q_rad_count] for i in 1:3]
radii_2q_long = [ [xdats_2q_long[i][j][2][1][1] for j in 1:q_rad_count] for i in 1:3]
theory_2q_long = [radii_2q_long[i].*radii_2q_long[i]./(2 .*(i*2+1)) for i in 1:3]
=#



#
for i in 1:1
	ms = [3,5,7]
	m_here = ms[i]
	filling = 1/m_here
	rm1 = sqrt(2*particles/filling)
	#plot(radii_origin[i],berry_phase_origin[i],"-p",label="M=$m_here")
	#plot(radii_1q[i],berry_phase_1q[i],"-p",label="M=$m_here")
	errors = berry_phase_origin_long_errors[i]
	calced_m = (theory_2q_long[i]-berry_phase_origin_long[i])#.^(-1)
	#errors_m = ((theory_2q_long[i]-berry_phase_origin_long[i])+errors).^(-1) - calced_m 
	#plot(radii_2q_long[i],berry_phase_origin_long[i],"-p",label="M=$m_here")
	figure()
	plot(radii_2q_long[i]./rm1,[1/m_here for i in 1:9],"-k",label="TH")
	errorbar(radii_2q_long[i]./rm1,calced_m,yerr=errors,c="red",label="1/$m_here")
	ylim((0.0,1.5))
	legend(loc="upper left")
	ylabel("Calculated Filling Factor")
	xlabel("Radius of Quasihole / Rm")
	title("Two QHs Berry Phase Interaction Term, N=20")
end
#
#



#phase_dats = get_phase_from_ideal_dtheta(xdats,ydats,berries[1],berries[2],1)
#plot(phase_dats[1][1])
#phase_data = get_phase_from_ideal_dtheta(xdats,ydats,berries[1],berries[2],3)
#plot(1:3,phase_data[1])
#plot(1:3,th_vals,label="Theory")
#legend()





# no distances less than 0.1 with 1000 steps for all data sets and all m's
# take time autocorrelation function from ModSim to see if 1000 is good, might be the issue


#=  OG running for local rotation of qhole
particles = 20
steps = 1000
rad = 0.001
xdats = [ [read_hdf5_data(particles,i,"x","qhole-data","mk2"),read_hdf5_data(particles,i,"x","qhole-data","rnd"), read_hdf5_data(particles,i,"x","qhole-data","og")] for i in 1:3]
ydats = [ [read_hdf5_data(particles,i,"y","qhole-data","mk2"),read_hdf5_data(particles,i,"y","qhole-data","rnd"), read_hdf5_data(particles,i,"y","qhole-data","og")] for i in 1:3]

th_vals = [-(rad^2)/(4*i) for i in 1:3]
berries = get_calc_berry(rad,xdats,ydats,steps,particles,3)
phase_data = get_phase_from_ideal_dtheta(xdats,ydats,berries[1],berries[2],3)
plot(1:3,phase_data[1])
plot(1:3,th_vals,label="Theory")
legend()
=#

#


"fin"
