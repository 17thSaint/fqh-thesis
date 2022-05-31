#import Pkg; Pkg.add("HDF5")
using HDF5

function make_vecovecs_matrix(given_data::Vector{Vector{Any}})
	matrix_version::Matrix{Int64} = fill(0,(length(given_data),length(given_data[1][2])+1))
	for i in 1:length(given_data)
		flat_given_data = collect(Iterators.flatten(given_data[i]))
		for j in 1:length(flat_given_data)
			matrix_version[i,j] = flat_given_data[j]
		end
	end
	return matrix_version
end

function make_matrix_vecovecs(given_data::Matrix{Int64})
	vec_version::Vector{Vector{Int64}} = [[0 for i in 1:size(given_data)[2]] for j in 1:size(given_data)[1]]
	for i in 1:size(given_data)[2]
		for j in 1:size(given_data)[1]
			vec_version[j][i] = given_data[j,i]
		end
	end
	return vec_version
end

function make_matrix_vecovecs_compressed(given_data::Matrix{Int64})
	vec_version::Vector{Vector{Any}} = [[0,[0 for i in 1:size(given_data)[2]-1]] for j in 1:size(given_data)[1]]
	for i in 1:size(given_data)[1]
		vec_version[i][1] = given_data[i,1]
		vec_version[i][2] = given_data[i,2:end]
	end
	return vec_version
end

function rewrite_acc_data(particles::Int64)
	cd("..")
	cd("acc-matrix-data")
	
	og_file = h5open("AccMat-parts-$particles-old.hdf5","r")
	new_file = h5open("AccMat-parts-$particles.hdf5","w")
	for i in 1:particles-1
		create_group(new_file,"order-$i")
		order_data = new_file["order-$i"]
		for j in 1:particles
			data = read(og_file["part-$j"],"ord-$i")
			order_data["part-$j"] = data
		end
	end
	close(og_file)
	close(new_file)
	cd("..")
	cd("Codes")
end

function write_acc_matrix_data(folder::String,particles::Int64,part::Int64,order::Int64,data::Matrix{Int64})
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if !isfile("AccMat-parts-$particles-compressed.hdf5")
		println("Making New Data File for N=$particles")
		binary_file = h5open("AccMat-parts-$particles-compressed.hdf5","w")
		for i in 1:particles
			create_group(binary_file,"part-$i")
		end
		part_data = binary_file["part-$part"]
		part_data["ord-$order"] = data
	else
		binary_file = h5open("AccMat-parts-$particles-compressed.hdf5","cw")
		part_data = binary_file["part-$part"]
		part_data["ord-$order"] = data
	end
	close(binary_file)
	if folder != "NA"
		cd("..")
		cd("Codes")
	end
end

function read_acc_matrix_data(folder::String,particles::Int64,part,order,which)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	all_data = Vector{Vector{Vector{Any}}}(undef,length(order))
	binary_file = h5open("AccMat-parts-$particles-$which.hdf5","r")
	data = read(binary_file["part-$part"])#read(binary_file["order-$order"])
	index = 0
	if which == "old"
		for i in order
			index += 1
			all_data[index] = make_matrix_vecovecs(data["ord-$i"])
		end
	elseif which == "compressed"
		for i in order
			index += 1
			all_data[index] = make_matrix_vecovecs_compressed(data["ord-$i"])
		end
	end
	close(binary_file)
	if folder != "NA"
		cd("..")
		cd("Codes")
	end
	if typeof(order) == Int64
		return all_data[1]
	end
	return all_data
end

function get_full_acc_matrix(num_parts::Int64,which="compressed")
	acc_mat = Matrix{Any}(undef,num_parts,num_parts-1)
	my_orders = [i for i in 1:num_parts-1]
	for which_part in 1:num_parts
		acc_mat[which_part,:] = read_acc_matrix_data("acc-matrix-data",num_parts,which_part,my_orders,which)
	end
	return acc_mat
end

function compress_acc_set(parts_count::Int64,acc_set::Vector{Vector{Int64}})
	new_acc_set = [ [1,[1 for j in 1:length(acc_set[1])]] for i in 1:binomial(parts_count-1,length(acc_set[1]))]
	index = 1
	for i in unique(acc_set)
		found = findall(x->x==i,acc_set)
		new_acc_set[index] = [length(found),i]
		index += 1
	end
	return new_acc_set
end



#=
chosen_parts = 9#parse(Int64,ARGS[1])
for part in 1:chosen_parts
	for order in 1:chosen_parts-1
		#calced_list = get_all_acc_sets(order,part,chosen_parts)
		#matrix_list = make_vecovecs_matrix(calced_list)
		#old_data = read_acc_matrix_data("acc-matrix-data",chosen_parts,part,order,"old")
		compressed_version = compress_acc_set(chosen_parts,old_data)
		matrix_comp_ver = make_vecovecs_matrix(compressed_version)
		write_acc_matrix_data("acc-matrix-data",chosen_parts,part,order,matrix_comp_ver)
	end
end
=#

"fin"
