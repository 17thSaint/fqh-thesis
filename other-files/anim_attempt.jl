function sin1(k)
	#return [[main_data[1][(k-1)*particles+l],main_data[2][(k-1)*particles+l]] for l in 1:particles]
	return [[1,2],[3,4]]
end

using Plots
anim = @animate for i in 1:10
	plot(sin1, 0, i)
end

gif(anim, fps=15)
