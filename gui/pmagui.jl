module PMAGUI

using PrincipalMomentAnalysis
using LinearAlgebra
using Statistics
using DataFrames

using Blink
using JSExpr
using DataStructures

using Colors
using PlotlyJS
using Measures

using DelimitedFiles

include("qlucore.jl")
include("pmaplots.jl")
include("pmaplots_gadfly.jl")
include("pca.jl")

struct Dataset
	data::Matrix
	sa::DataFrame
	va::DataFrame
	filepath::String
	errorMsg::String
end
Dataset() = Dataset(zeros(0,0),DataFrame(),DataFrame(),"","")
Dataset(data::Matrix, sa::DataFrame, va::DataFrame,filepath::String) = Dataset(data,sa,va,filepath,"")
Dataset(errorMsg::String) = Dataset(zeros(0,0),DataFrame(),DataFrame(),"",errorMsg)

function opendataset(filepath::String)::Dataset
	isempty(filepath) && return Dataset("Please select file.")
	isfile(filepath) || return Dataset("File not found: \"$filepath\".")
	try
		data,sa,va = Qlucore.read(filepath)
		Dataset(data,sa,va,filepath)
	catch e
		Dataset("Error loading file: \"$filepath\" ($e)")
	end
end


struct Result
	U::Matrix{Float64}
	S::Vector{Float64}
	V::Matrix{Float64}
	G::AbstractMatrix # the graph (needed for plotting)
end
Result() = Result(zeros(0,0),zeros(0),zeros(0,0),zeros(0,0))


function populategui!(w, ds::Dataset)
	infoStr = ""
	sa = String[]
	va = String[]

	if isempty(ds.errorMsg)
		infoStr = """File: "$(ds.filepath)" with $(size(ds.data,2)) samples, $(size(ds.data,1)) variables, $(size(ds.sa,2)) sample annotations and $(size(ds.va,2)) variable annotations."""
		sa = string.(names(ds.sa))
		va = string.(names(ds.va))
	else
		infoStr = ds.errorMsg
	end

	js(w, js"""document.getElementById("info").innerHTML = '$infoStr';""")

	# TODO: simple implementation, optimize
	js(w, js"""document.getElementById("sampleAnnot").options.length = 0;""")
	for a in sa
		js(w, js"""var opt = document.createElement("option");
			       opt.value = $a;
			       opt.text = $a;
			       document.getElementById("sampleAnnot").options.add(opt);""")
	end

	# TODO: simple implementation, optimize
	js(w, js"""document.getElementById("timeAnnot").options.length = 0;""")
	for a in sa
		js(w, js"""var opt = document.createElement("option");
			       opt.value = $a;
			       opt.text = $a;
			       document.getElementById("timeAnnot").options.add(opt);""")
	end

	# TODO: simple implementation, optimize
	js(w, js"""document.getElementById("exportVarAnnot").options.length = 0;""")
	for a in va
		js(w, js"""var opt = document.createElement("option");
			       opt.value = $a;
			       opt.text = $a;
			       document.getElementById("exportVarAnnot").options.add(opt);""")
	end

end


function changesamplemethod!(w, sampleMethod::Symbol)
	enabled = sampleMethod in (:SA,:Time,:NNSA)
	js(w, js"""document.getElementById("sampleAnnot").disabled = $(!enabled);""")

	enabled = sampleMethod == :Time
	js(w, js"""document.getElementById("timeAnnot").disabled = $(!enabled);""")

	enabled = sampleMethod in (:NN,:NNSA)
	js(w, js"""document.getElementById("kNearestNeighbors").disabled = $(!enabled);""")
	js(w, js"""document.getElementById("distNearestNeighbors").disabled = $(!enabled);""")
end

function rundimreduction(ds::Dataset;
                         dimReductionMethod::Symbol, sampleMethod::Symbol, sampleAnnotation::Symbol,
                         timeAnnotation::Symbol, kNearestNeighbors::Int, distNearestNeighbors::Float64)
	isempty(ds.errorMsg) || return
	@assert dimReductionMethod in (:PMA,:PCA)
	@assert sampleMethod in (:SA,:Time,:NN,:NNSA)
	# @assert plotDims in 2:3

	X = zeros(size(ds.data))
	if any(ismissing,ds.data)
		# Replace missing values with mean over samples with nonmissing data
		println("Reconstructing missing values (taking the mean over all nonmissing samples)")
		for i=1:size(X,1)
			m = ismissing.(ds.data[i,:])
			X[i,.!m] .= ds.data[i,.!m]
			X[i,m] .= mean(ds.data[i,.!m])
		end
	else
		X .= ds.data # just copy
	end

	normalizemeanstd!(X)

	G = nothing
	if sampleMethod == :SA
		G = buildgraph(ds.sa[!,sampleAnnotation])
	elseif sampleMethod == :Time
		eltype(ds.sa[!,timeAnnotation])<:Number || @warn "Expected time annotation to contain numbers, got $(eltype(ds.sa[!,timeAnnotation])). Fallback to default sorting."
		G = buildgraph(ds.sa[!,sampleAnnotation], ds.sa[!,timeAnnotation])
	elseif sampleMethod == :NN
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50);
	elseif sampleMethod == :NNSA
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50; groupBy=ds.sa[!,sampleAnnotation]);
	end

	# dim = 3
	dim = min(10, size(X)...)

	if dimReductionMethod==:PMA
		U,S,V = pma(X,G,dim=dim)
	elseif dimReductionMethod==:PCA
		U,S,V = pca(X,dim=dim)
	end
	Result(U,S,V,G)
end

function showplot(ds::Dataset, result::Result, plotDims::Int;
                  dimReductionMethod::Symbol, sampleAnnotation::Symbol,
                  kwargs...)
	colorBy = sampleAnnotation
	colorDict = colordict(ds.sa[!,colorBy])

	println(collect(keys(colorDict)))
	println(collect(values(colorDict)))

	opacity = 0.05
	title = string(splitdir(ds.filepath)[2], " ", string(dimReductionMethod))
	drawTriangles = false#true
	drawLines = true

	if plotDims==3
		markerSize = 5
		lineWidth = 1

		pl = plotsimplices(result.V,ds.sa,result.G,colorBy,colorDict, title=title,
		                   drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                   opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
		                   width=1024, height=768)
		display(pl)
	elseif plotDims==2
		markerSize = 0.9mm
		lineWidth = 0.3mm

		pl = plotsimplices_gadfly(result.V,ds.sa,result.G,colorBy,colorDict, title=title,
		                          drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                          opacity=opacity, markerSize=markerSize, lineWidth=lineWidth)
		display(pl)
	end
end


function exportloadings(ds::Dataset, result::Result, filepath::String, varAnnot::Symbol, columns::Symbol, dims::Int, columnSort::Symbol;
                        dimReductionMethod::Symbol, sampleAnnotation::Symbol,
                        kwargs...)
	@assert columns in (:All, :Selected)
	@assert columnSort in (:Abs, :Descending)

	df = DataFrame()
	df[!,varAnnot] = ds.va[!,varAnnot]

	if columns==:All
		for i=1:min(dims,length(result.S))
			colName = Symbol("$dimReductionMethod$i")
			df[!,colName] = result.U[:,i]
		end
	elseif columns==:Selected
		colName = Symbol("$dimReductionMethod$dims")
		df[!,colName] = result.U[:,dims]

		if columnSort==:Abs
			sort!(df, colName, by=abs, rev=true)
		elseif columnSort==:Descending
			sort!(df, colName, rev=true)
		end
	end

	# TODO: add CSV dependency and save using CSV.write() instead?
	mat = Matrix{String}(undef, size(df,1)+1, size(df,2))
	mat[1,:] .= string.(names(df))
	mat[2:end,:] .= string.(df)
	writedlm(filepath, mat, '\t')
end


function main()
	# init
	dataset = opendataset("")
	dimReductionParams = nothing
	result = nothing
	messageQueue = Queue{Pair{String,Any}}()


	# setup gui
	w = Window(Dict(:width=>512,:height=>768))

	# event listeners
	handle(w, "gedataopen") do args
		println("gedataopen: ", args)
		fn = isempty(args) ? "" : args[1]
		enqueue!(messageQueue, "gedataopen"=>fn)
	end
	handle(w, "samplemethod") do args
		println("samplemethod: ", args)
		enqueue!(messageQueue, "samplemethod"=>args)
	end
	handle(w, "showplot") do args
		println("showplot: ", args)
		enqueue!(messageQueue, "showplot"=>args)
	end
	handle(w, "exportloadings") do args
		println("exportloadings: ", args)
		enqueue!(messageQueue, "exportloadings"=>args)
	end

	doc = read(joinpath(@__DIR__,"content.html"),String)
	body!(w,doc,async=false)
	changesamplemethod!(w,:SA)


	# Message handling loop
	while isopen(w.content.sock) # is there a better way to check if the window is still open?
		if !isempty(messageQueue)
			msg = dequeue!(messageQueue)
			if msg.first == "gedataopen"
				filepath = msg.second
				println("Loading ", filepath)
				dataset = opendataset(filepath)
				dimReductionParams = nothing
				result = nothing
				if isempty(dataset.errorMsg)
					println(filepath, " loaded successfully")
				else
					println("Loading failed with error, ", dataset.errorMsg)
				end
				populategui!(w,dataset)
			elseif msg.first == "samplemethod"
				sampleMethod = Symbol(msg.second)
				println("Changing method to ", sampleMethod)
				changesamplemethod!(w,sampleMethod)
			elseif msg.first == "showplot"
				args = msg.second
				dimReductionMethod = Symbol(args[1])

				params = Dict()
				try
					params = Dict(:dimReductionMethod=>dimReductionMethod,
					              :sampleMethod=>Symbol(args[2]),
					              :sampleAnnotation=>Symbol(args[3]),
					              :timeAnnotation=>Symbol(args[4]),
					              :kNearestNeighbors=>parse(Int,args[5]),
					              :distNearestNeighbors=>parse(Float64,args[6]))
				catch e
					println("Failed to parse parameters ", args, ": ", sprint(showerror, e))
				end

				if params != dimReductionParams
					try
						println("Running ", dimReductionMethod)
						result = rundimreduction(dataset; params...)
						println("Done")
						dimReductionParams = params
					catch e
						println("Failed to run ", dimReductionMethod, ": ", sprint(showerror, e))
						result = nothing
						dimReductionParams = nothing
					end
				end

				if result!=nothing
					try
						println("Showing plot")
						plotDims = parse(Int,args[7])
						showplot(dataset, result, plotDims; dimReductionParams...)
					catch e
						println("Error showing plot: ", sprint(showerror,e))
					end
				end
			elseif msg.first == "exportloadings"
				# TODO: avoid code duplication
				args = msg.second
				dimReductionMethod = Symbol(args[1])

				params = Dict()
				try
					params = Dict(:dimReductionMethod=>dimReductionMethod,
					              :sampleMethod=>Symbol(args[2]),
					              :sampleAnnotation=>Symbol(args[3]),
					              :timeAnnotation=>Symbol(args[4]),
					              :kNearestNeighbors=>parse(Int,args[5]),
					              :distNearestNeighbors=>parse(Float64,args[6]))
				catch e
					println("Failed to parse parameters ", args, ": ", sprint(showerror, e))
				end

				if params != dimReductionParams
					try
						println("Running ", dimReductionMethod)
						result = rundimreduction(dataset; params...)
						println("Done")
						dimReductionParams = params
					catch e
						println("Failed to run ", dimReductionMethod, ": ", sprint(showerror, e))
						result = nothing
						dimReductionParams = nothing
					end
				end

				if result!=nothing
					try
						println("Exporting loadings")
						println(args)
						filepath = args[7]
						varAnnot = Symbol(args[8])
						columns = Symbol(args[9])
						dims = parse(Int,args[10])
						columnSort = Symbol(args[11])
						exportloadings(dataset, result, filepath, varAnnot, columns, dims, columnSort; dimReductionParams...)
						println("Done")
					catch e
						println("Error exporting loadings: ", sprint(showerror,e))
					end
				end
			else
				@warn "Unknown message type: $(msg.first)"
			end
		end

		# yield() # Allow GUI to run
		sleep(0.05) # Allow GUI to run
	end


end


end

PMAGUI.main()