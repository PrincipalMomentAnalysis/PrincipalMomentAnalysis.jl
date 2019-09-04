using Colors
import Gadfly; const gf = Gadfly
using IterTools

function plotsimplices_gadfly(V, sa, G, colorBy, colorDict; drawPoints=true, drawLines=true, drawTriangles=true, title="", opacity=0.3, markerSize=5, lineWidth=2, shapeBy=nothing, colorTitle=nothing, shapeTitle=nothing)
	layers = []

	# plot each group in different colors
	if drawPoints
		for cb in unique(sa[colorBy])
			ind = findall( sa[colorBy].==cb )
			col = colorDict[cb]

			extras = []
			shapeBy==nothing || push!(extras, (shape=sa[ind,shapeBy],))
			isempty(extras) || (extras = pairs(extras...))

			points = gf.layer(x=V[ind,1],y=V[ind,2], gf.Geom.point, gf.style(default_color=col, point_size=markerSize); extras... )
			push!(layers, points)
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
			
			line = gf.layer(x=V[ind,1],y=V[ind,2],gf.Geom.path, gf.style(default_color=col, line_width=lineWidth))

			push!(layers, line)
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

			triangle = gf.layer(x=V[ind,1],y=V[ind,2], gf.Geom.polygon(preserve_order=true, fill=true), gf.style(default_color=RGBA{Float32}(col.r,col.g,col.b,opacity)))
			push!(layers, triangle)
		end
	end

	extras = []
	colorTitle==nothing || push!(extras, gf.Guide.manual_color_key(colorTitle, collect(keys(colorDict)), collect(values(colorDict))))
	shapeTitle==nothing || push!(extras, gf.Guide.shapekey(shapeTitle,[""],[]))

	# gf.plot(layers..., gf.Theme(background_color=colorant"white", key_position=:none))
	gf.plot(layers..., extras..., gf.Theme(background_color=colorant"white"))
end



nothing