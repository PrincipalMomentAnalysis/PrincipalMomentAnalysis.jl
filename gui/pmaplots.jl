function plotsimplices(V, sa, G, colorBy, colorDict;
	                   drawPoints=true, drawLines=true, drawTriangles=true,
	                   title="",
	                   opacity=0.3, markerSize=5, lineWidth=2,
	                   shapeBy=nothing, shapeDict=nothing,
	                   width=1536, height=768)
	traces = GenericTrace[]

	# plot each group in different colors
	if drawPoints


		# testing plotting with colorscale instead of dict when colorDict=nothing.
		# TODO: merge the two cases below and implement for lines/triangles too

		if colorDict != nothing
			for cb in unique(sa[!,colorBy])
				ind = findall( sa[!,colorBy].==cb )
				col = colorDict[cb]

				extras = []
				shapeBy!=nothing && shapeDict!=nothing && push!(extras, (marker_symbol=[shapeDict[k] for k in sa[ind,shapeBy]],))
				isempty(extras) || (extras = pairs(extras...))

				points = scatter3d(;x=V[ind,1],y=V[ind,2],z=V[ind,3], mode="markers", marker_color=col, marker_size=markerSize, marker_line_width=0, name=string(cb), extras...)
				push!(traces, points)
			end
		else
			extras = []
			shapeBy!=nothing && shapeDict!=nothing && push!(extras, (marker_symbol=[shapeDict[k] for k in sa[!,shapeBy]],))
			isempty(extras) || (extras = pairs(extras...))

			# points = scatter3d(;x=V[:,1],y=V[:,2],z=V[:,3], mode="markers", marker_color=sa[!,colorBy], marker_size=markerSize, marker_line_width=0, name=string(colorBy), extras...)
			points = scatter3d(;x=V[:,1],y=V[:,2],z=V[:,3], mode="markers", marker=attr(color=sa[!,colorBy], colorscale="Viridis", showscale=true, size=markerSize, line_width=0), name=string(colorBy), extras...)
			push!(traces, points)
		end


	end

	if drawLines
		# plot lines based on graph
		for ci in findall(G)
			r,c = Tuple(ci)
			r>c || continue # just use lower triangular part
			cb = sa[r,colorBy] # use color from the first one, should be same
			col = colorDict[cb]
			ind = [r,c]
			line = scatter3d(;x=V[ind,1],y=V[ind,2],z=V[ind,3], mode="lines", line=attr(color=col, width=lineWidth), showlegend=false)
			push!(traces, line)
		end
	end

	if drawTriangles
		# plot triangles base on graph - TODO: improve code!
		triangleInds = Int[]
		for c=1:size(G,1)
			ind = findall(G[:,c])
			isempty(ind) && continue

			length(ind)<3 && continue # no triangles
			# length(ind)>3 && error("Tetrahedrons not yet supported") # tetrahedrons are above, TODO implement

			# slow and ugly solution
			for tri in subsets(ind,3)
				append!(triangleInds, sort(tri))
			end
		end
		triangleInds = reshape(triangleInds,3,:)
		triangleInds = unique(triangleInds,dims=2) # remove duplicates

		for i=1:size(triangleInds,2)
			ind = triangleInds[:,i]
			cb = sa[ind[1],colorBy] # use color from the first one, should be same
			col = colorDict[cb]

			triangle = mesh3d(;x=V[ind,1],y=V[ind,2],z=V[ind,3], color=col, opacity=opacity, showlegend=false)
			push!(traces, triangle)
		end
	end


	layout = Layout(autosize=false, width=width, height=height, margin=attr(l=0, r=0, b=0, t=65), title=title)
	# layout = Layout(margin=attr(l=0, r=0, b=0, t=65), title=title)
	#plot(traces, layout)
	traces, layout # return plot args rather than plot because of threading issues.
end

_distinguishable_colors(n) = distinguishable_colors(n+1,colorant"white")[2:end]

function colordict(x::AbstractArray)
	k = unique(x)
	Dict(s=>c for (s,c) in zip(k,_distinguishable_colors(length(k))))
end
