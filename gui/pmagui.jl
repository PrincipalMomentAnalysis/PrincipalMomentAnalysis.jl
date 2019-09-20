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

function populategui!(w, ds::Dataset)
	infoStr = ""
	sa = String[]

	if isempty(ds.errorMsg)
		infoStr = """File: "$(ds.filepath)" with $(size(ds.data,2)) samples, $(size(ds.data,1)) variables, $(size(ds.sa,2)) sample annotations and $(size(ds.va,2)) variable annotations."""
		sa = string.(names(ds.sa))
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

	# js(w, js"""console.log('$infoStr');""")	
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


function runpma(ds::Dataset, sampleMethod::Symbol, sampleAnnotation::Symbol, timeAnnotation::Symbol, kNearestNeighbors::Int, distNearestNeighbors::Float64, plotDims::Int)
	isempty(ds.errorMsg) || return
	@assert sampleMethod in (:SA,:Time,:NN,:NNSA)
	@assert plotDims in 2:3

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
		G = buildgraph(ds.sa[!,sampleAnnotation], ds.sa[!,timeAnnotation])
	elseif sampleMethod == :NN
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50);
	elseif sampleMethod == :NNSA
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50; groupBy=ds.sa[!,sampleAnnotation]);
	end

	dim = plotDims
	UPMA,SPMA,VPMA = pma(X,G,dim=dim)

	colorBy = sampleAnnotation
	colorDict = colordict(ds.sa[!,colorBy])

	println(collect(keys(colorDict)))
	println(collect(values(colorDict)))

	opacity = 0.05
	title = splitdir(ds.filepath)[2]
	drawTriangles = false#true
	drawLines = true

	if plotDims==3
		markerSize = 5
		lineWidth = 1

		plPMA = plotsimplices(VPMA,ds.sa,G,colorBy,colorDict, title="$title PMA",
		                      drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                      opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
		                      width=1024, height=768)
		display(plPMA)
	elseif plotDims==2
		markerSize = 0.9mm
		lineWidth = 0.3mm

		plPMA = plotsimplices_gadfly(VPMA,ds.sa,G,colorBy,colorDict, title="$title PMA",
		                             drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                             opacity=opacity, markerSize=markerSize, lineWidth=lineWidth)
		display(plPMA)
	end
end

function runpca(ds::Dataset, sampleMethod::Symbol, sampleAnnotation::Symbol, timeAnnotation::Symbol, kNearestNeighbors::Int, distNearestNeighbors::Float64, plotDims::Int)
	isempty(ds.errorMsg) || return
	@assert sampleMethod in (:SA,:Time,:NN,:NNSA)
	@assert plotDims in 2:3

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
		G = buildgraph(ds.sa[!,sampleAnnotation], ds.sa[!,timeAnnotation])
	elseif sampleMethod == :NN
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50);
	elseif sampleMethod == :NNSA
		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50; groupBy=ds.sa[!,sampleAnnotation]);
	end

	dim = plotDims
	UPCA,SPCA,VPCA = pca(X,dim=dim)

	colorBy = sampleAnnotation
	colorDict = colordict(ds.sa[!,colorBy])

	println(collect(keys(colorDict)))
	println(collect(values(colorDict)))

	opacity = 0.05
	title = splitdir(ds.filepath)[2]
	drawTriangles = false#true
	drawLines = true

	if plotDims==3
		markerSize = 5
		lineWidth = 1

		plPCA = plotsimplices(VPCA,ds.sa,G,colorBy,colorDict, title="$title PCA",
		                      drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                      opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
		                      width=1024, height=768)
		display(plPCA)
	elseif plotDims==2
		markerSize = 0.9mm
		lineWidth = 0.3mm

		plPCA = plotsimplices_gadfly(VPCA,ds.sa,G,colorBy,colorDict, title="$title PCA",
		                             drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
		                             opacity=opacity, markerSize=markerSize, lineWidth=lineWidth)
		display(plPCA)
	end
end



function main()
	# init
	dataset = opendataset("")
	messageQueue = Queue{Pair{String,Any}}()


	# setup gui
	w = Window(Dict(:width=>512,:height=>512))

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
	handle(w, "runpma") do args
		println("runpma: ", args)
		enqueue!(messageQueue, "runpma"=>args)
	end
	handle(w, "runpca") do args
		println("runpca: ", args)
		enqueue!(messageQueue, "runpca"=>args)
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
			elseif msg.first == "runpma"
				println("Running PMA")
				try
					isempty(dataset.errorMsg) || error(dataset.errorMsg)
					args = msg.second
					sampleMethod = Symbol(args[1])
					sampleAnnotation = Symbol(args[2])
					timeAnnotation = Symbol(args[3])
					kNearestNeighbors = parse(Int,args[4])
					distNearestNeighbors = parse(Float64,args[5])
					plotDims = parse(Int,args[6])
					runpma(dataset,sampleMethod, sampleAnnotation, timeAnnotation, kNearestNeighbors, distNearestNeighbors, plotDims)
				catch e
					println("Failed to run PMA: ", sprint(showerror, e))
				end
			elseif msg.first == "runpca" # TODO: merge with PMA code above?
				println("Running PCA")
				try
					isempty(dataset.errorMsg) || error(dataset.errorMsg)
					args = msg.second
					sampleMethod = Symbol(args[1])
					sampleAnnotation = Symbol(args[2])
					timeAnnotation = Symbol(args[3])
					kNearestNeighbors = parse(Int,args[4])
					distNearestNeighbors = parse(Float64,args[5])
					plotDims = parse(Int,args[6])
					runpca(dataset,sampleMethod, sampleAnnotation, timeAnnotation, kNearestNeighbors, distNearestNeighbors, plotDims)
				catch e
					println("Failed to run PCA: ", sprint(showerror, e))
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