#import Pkg; Pkg.add("HDF5")
using HDF5

include("cf-wavefunc.jl")

function make_vecovecs_matrix(given_data::Vector{Vector{Int64}})
	matrix_version::Matrix{Int64} = fill(0,(length(given_data),length(given_data[1])))
	for i in 1:length(given_data)
		for j in 1:length(given_data[1])
			matrix_version[i,j] = given_data[i][j]
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

function write_acc_matrix_data(folder::String,particles::Int64,part::Int64,order::Int64,data::Matrix{Int64})
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	if !isfile("AccMat-parts-$particles.hdf5")
		println("Making New Data File for N=$particles")
		binary_file = h5open("AccMat-parts-$particles.hdf5","w")
		for i in 1:particles
			create_group(binary_file,"part-$i")
		end
		part_data = binary_file["part-$part"]
		part_data["ord-$order"] = data
	else
		binary_file = h5open("AccMat-parts-$particles.hdf5","cw")
		part_data = binary_file["part-$part"]
		part_data["ord-$order"] = data
	end
	close(binary_file)
end

function read_acc_matrix_data(folder::String,particles::Int64,part::Int64,order::Int64)
	if folder != "NA"
		cd("..")
		cd("$folder")
	end
	binary_file = h5open("AccMat-parts-$particles.hdf5","r")
	data = read(binary_file["part-$part"],"ord-$order")
	vecovecs_data = make_matrix_vecovecs(data)
	close(binary_file)
	return vecovecs_data
end




"fin"
