module PMAGUI

using PrincipalMomentAnalysis
using LinearAlgebra
# using Statistics
using DataFrames
using CSV

using Blink
using JSExpr
# using DataStructures

using Colors
using PlotlyJS
# using Measures
using IterTools

# using DelimitedFiles

include("Schedulers.jl")
using .Schedulers

include("qlucore.jl")
include("pmaplots.jl")
#include("pmaplots_gadfly.jl")
include("pca.jl")


struct JobGraph
	scheduler::Scheduler
	sampleStatus::Ref{String}
	fileIOLock::ReentrantLock
	paramIDs::Dict{String,JobID}
	loadSampleID::JobID
	normalizeID::JobID
	setupGraphID::JobID
	dimreductionID::JobID
	makeplotID::JobID
end

struct SampleData
	data::Matrix
	sa::DataFrame
	va::DataFrame
	filepath::String
end
SampleData() = SampleData(zeros(0,0),DataFrame(),DataFrame(),"")

struct ReducedSampleData
	F::Factorization
	sa::DataFrame
	va::DataFrame
end


function loadcsv(filepath::String; delim, nbrSampleAnnots, transpose::Bool=false)
	df = CSV.read(filepath; delim=delim, transpose=transpose)
	sa = df[:, 1:nbrSampleAnnots]
	va = DataFrame(VariableID=names(df)[nbrSampleAnnots+1:end])
	data = convert(Matrix, df[!,nbrSampleAnnots+1:end])'
	@assert eltype(data) <: Union{Number,Missing}
	data,sa,va
end

function loadsample(st, input::Dict{String,Any})::SampleData
	@assert length(input)==3
	filepath = input["filepath"]::String
	nbrSampleAnnots = parse(Int,input["nbrsampleannots"]::String)
	rowsAsSamples   = parse(Bool,input["rowsassamples"])
	filepath == :__INVALID__ && return Nothing
	filepath::String
	isempty(filepath) && return Nothing
	@assert isfile(filepath) "Sample file not found: \"$filepath\""
	ext = lowercase(splitext(filepath)[2])
	if ext==".gedata"
		originalData,sa,va = Qlucore.read(filepath)
	elseif ext in (".csv",".tsv",".txt")
		originalData,sa,va = loadcsv(filepath; delim=ext==".csv" ? ',' : '\t', nbrSampleAnnots=nbrSampleAnnots, transpose=!rowsAsSamples)
	end

	data = zeros(size(originalData))
	if any(ismissing,originalData)
		# Replace missing values with mean over samples with nonmissing data
		println("Reconstructing missing values (taking the mean over all nonmissing samples)")
		for i=1:size(data,1)
			m = ismissing.(originalData[i,:])
			data[i,.!m] .= originalData[i,.!m]
			data[i,m] .= mean(originalData[i,.!m])
		end
	else
		data .= originalData # just copy
	end

	SampleData(data,sa,va,filepath)
end


# callback function
showsampleannotnames(sampleData, toGUI) = put!(toGUI, :displaysampleannotnames=>names(sampleData.sa))


function normalizesample(st, input::Dict{String,Any})
	@assert length(input)==2
	sampleData = input["sampledata"]
	method = input["method"]
	@assert method in ("None", "Mean=0", "Mean=0,Std=1")
	X = sampleData.data
	if method == "Mean=0,Std=1"
		X = normalizemeanstd(X)
	elseif method == "Mean=0"
		X = normalizemean(X)
	end
	SampleData(X,sampleData.sa,sampleData.va,sampleData.filepath)
end


function setupgraph(st, input::Dict{String,Any})
	@assert length(input)==6
	sampleData  = input["sampledata"]
	method      = Symbol(input["method"])
	sampleAnnot = Symbol(input["sampleannot"])
	timeAnnot   = Symbol(input["timeannot"])
	kNN         = parse(Int,input["knearestneighbors"])
	distNN      = parse(Float64, input["distnearestneighbors"])
	@assert method in (:SA,:Time,:NN,:NNSA)

	G = nothing
	if method == :SA
		G = buildgraph(sampleData.sa[!,sampleAnnot])
	elseif method == :Time
		eltype(sampleData.sa[!,timeAnnot])<:Number || @warn "Expected time annotation to contain numbers, got $(eltype(sampleData.sa[!,timeAnnot])). Fallback to default sorting."
		G = buildgraph(sampleData.sa[!,sampleAnnot], sampleData.sa[!,timeAnnot])
	elseif method == :NN
		G = neighborhoodgraph(sampleData.data, kNN, distNN, 50);
	elseif method == :NNSA
		G = neighborhoodgraph(sampleData.data, kNN, distNN, 50; groupBy=sampleData.sa[!,sampleAnnot]);
	end
	G
end


function dimreduction(st, input::Dict{String,Any})
	@assert length(input)==3
	sampleData  = input["sampledata"]
	sampleGraph = input["samplegraph"]
	method      = Symbol(input["method"])

	X = sampleData.data::Matrix{Float64}

	# dim = 3
	dim = min(10, size(X)...)

	if method==:PMA
		F = pma(X, sampleGraph, nsv=dim)
	elseif method==:PCA
		F = svdbyeigen(X,nsv=dim)
	end
	ReducedSampleData(F, sampleData.sa, sampleData.va)
end

function makeplot(st, input::Dict{String,Any})
	@assert length(input)==11
	reduced            = input["reduced"]::ReducedSampleData
	dimReductionMethod = Symbol(input["dimreductionmethod"])
	sampleAnnot        = Symbol(input["sampleannot"])
	sampleGraph        = input["samplegraph"]
	plotDims           = parse(Int,input["plotdims"])
	plotWidth          = parse(Int,input["plotwidth"])
	plotHeight         = parse(Int,input["plotheight"])
	markerSize         = parse(Float64,input["markersize"])
	showPoints         = parse(Bool,input["showpoints"])
	showLines          = parse(Bool,input["showlines"])
	showTriangles      = parse(Bool,input["showtriangles"])


	@assert plotDims==3 "Only 3 plotting dims supported for now"

	opacity = 0.05
	title = dimReductionMethod

	colorBy = Symbol(sampleAnnot)

	# TODO: handle missing values in sample annotations?

	colorDict = nothing
	if !(eltype(reduced.sa[!,colorBy]) <: Real)
		colorDict = colordict(reduced.sa[!,colorBy])
	end

	plotArgs = nothing
	if plotDims==3
		lineWidth = 1
		plotArgs = plotsimplices(reduced.F.V,reduced.sa,sampleGraph,colorBy,colorDict, title=title,
		                         drawPoints=showPoints, drawLines=showLines, drawTriangles=showTriangles,
		                         opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
		                         width=plotWidth, height=plotHeight)
	end
	plotArgs
end

showplot(plotArgs, toGUI::Channel) = put!(toGUI, :displayplot=>plotArgs)


function samplestatus(jg::JobGraph)
	#jg.scheduler, jg.paramIDs["samplefilepath"] == :__INVALID__ && return "Please load sample."
	status = jobstatus(jg.scheduler, jg.loadSampleID)
	status==:done && return "Sample loaded."
	status in (:waiting,:running) && return "Loading sample."
	"Please load sample."
end

function setsamplestatus(jg::JobGraph, toGUI::Channel)
	sampleStatus = samplestatus(jg)
	sampleStatus!=jg.sampleStatus[] && put!(toGUI, :samplestatus=>sampleStatus)
	jg.sampleStatus[] = sampleStatus
end


function JobGraph()
	scheduler = Scheduler()
	# scheduler = Scheduler(threaded=false) # For DEBUG
	sampleIDs = Dict{String,Tuple{JobID,JobID}}()
	annotIDs  = Dict{String,Tuple{JobID,JobID}}()

	# Data Nodes (i.e. parameters chosen in the GUI)
	paramIDs = Dict{String,JobID}()

	loadSampleID = createjob!(loadsample, scheduler, "loadsample")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"samplefilepath")=>loadSampleID, "filepath")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"loadnbrsampleannots")=>loadSampleID, "nbrsampleannots")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"loadrowsassamples")=>loadSampleID, "rowsassamples")

	normalizeID = createjob!(normalizesample, scheduler, "normalizesample")
	add_dependency!(scheduler, loadSampleID=>normalizeID, "sampledata")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"normalizemethod")=>normalizeID, "method")

	setupGraphID = createjob!(setupgraph, scheduler, "setupgraph")
	add_dependency!(scheduler, normalizeID=>setupGraphID, "sampledata")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"samplegraphmethod")=>setupGraphID, "method")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"sampleannot")=>setupGraphID, "sampleannot")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"timeannot")=>setupGraphID, "timeannot")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"knearestneighbors")=>setupGraphID, "knearestneighbors")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"distnearestneighbors")=>setupGraphID, "distnearestneighbors")

	dimreductionID = createjob!(dimreduction, scheduler, "dimreduction")
	# add_dependency!(scheduler, loadSampleID=>dimreductionID, "sampledata")
	add_dependency!(scheduler, normalizeID=>dimreductionID, "sampledata")
	add_dependency!(scheduler, setupGraphID=>dimreductionID, "samplegraph")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"dimreductionmethod")=>dimreductionID, "method")

	makeplotID = createjob!(makeplot, scheduler, "makeplot")
	add_dependency!(scheduler, dimreductionID=>makeplotID, "reduced")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"dimreductionmethod")=>makeplotID, "dimreductionmethod")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"sampleannot")=>makeplotID, "sampleannot")
	add_dependency!(scheduler, setupGraphID=>makeplotID, "samplegraph")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"plotdims")=>makeplotID, "plotdims")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"plotwidth")=>makeplotID, "plotwidth")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"plotheight")=>makeplotID, "plotheight")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"markersize")=>makeplotID, "markersize")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"showpoints")=>makeplotID, "showpoints")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"showlines")=>makeplotID, "showlines")
	add_dependency!(scheduler, getparamjobid(scheduler,paramIDs,"showtriangles")=>makeplotID, "showtriangles")


	JobGraph(scheduler,
	         Ref(""),
	         ReentrantLock(),
	         paramIDs,
	         loadSampleID,
	         normalizeID,
	         setupGraphID,
	         dimreductionID,
	         makeplotID)
end

getparamjobid(s::Scheduler, paramIDs::Dict{String,JobID}, name::String, create::Bool=true) = create ? get!(paramIDs,name,createjob!(s, :__INVALID__, name)) : paramIDs[name]
getparamjobid(jg::JobGraph, name::String, args...) = getparamjobid(jg.scheduler, jg.paramIDs,name,args...)
setparam(jg::JobGraph, name::String, value) = setresult!(jg.scheduler, jg.paramIDs[name], value)


let jobGraph = JobGraph()
	global getjobgraph() = return jobGraph
end



function process_thread(fromGUI::Channel, toGUI::Channel)
	try
		# setup dependency graph
		jg = getjobgraph()
		scheduler = jg.scheduler
		lastSchedulerTime = UInt64(0)


		while true
			# @info "[Processing] tick"
			timeNow = time_ns()
			if isready(fromGUI)
				try
					msg = take!(fromGUI)
					msgName = msg.first
					msgArgs = msg.second
					@info "[Processing] Got message: $msgName $msgArgs"

					if msgName == :exit
						break
					elseif msgName == :cancel
						@info "[Processing] Cancelling all future events."
						cancelall!(scheduler)
					elseif msgName == :setvalue
						varName = msgArgs[1]
						value = msgArgs[2]
						if haskey(jg.paramIDs, varName)
							setresult!(scheduler, jg.paramIDs[varName], value)
						else
							@warn "Unknown variable name: $varName"
						end
					elseif msgName == :loadsample
						# schedule!(scheduler, jg.loadSampleID)
						schedule!(x->showsampleannotnames(x,toGUI), scheduler, jg.loadSampleID)
					# elseif msgName == :dimreduction
					# 	schedule!(scheduler, jg.dimreductionID)
					elseif msgName == :showplot
						schedule!(x->showplot(x,toGUI), scheduler, jg.makeplotID)
					else
						@warn "Unknown message type: $(msgName)"
					end
				catch e
					@warn "[Processing] Error processing GUI message."
					showerror(stdout, e, catch_backtrace())
				end
				setsamplestatus(jg, toGUI)
			elseif wantstorun(scheduler) || (isactive(scheduler) && (timeNow-lastSchedulerTime)/1e9 > 5.0)
				lastSchedulerTime = timeNow
				try
					step!(scheduler)
					status = statusstring(scheduler)
					@info "Job status: $status"
				catch e
					@warn "[Processing] Error processing event."
					showerror(stdout, e, catch_backtrace())
				end
				setsamplestatus(jg, toGUI)
			else
				sleep(0.05)
			end
		end

	catch e
		@warn "[Processing] Fatal error."
		showerror(stdout, e, catch_backtrace())
	end
	@info "[Processing] Exiting thread."
	put!(toGUI, :exited=>nothing)
end


function main()
	# This is the GUI thread

	@info "[PMAGUI] Using $(Threads.nthreads()) of $(Sys.CPU_THREADS) available threads."
	Threads.nthreads() == 1 && @warn "[PMAGUI] Threading not enabled, please set the environment variable JULIA_NUM_THREADS to the desired number of threads."

	# init
	fromGUI = Channel{Pair{Symbol,Vector}}(Inf)
	toGUI   = Channel{Pair{Symbol,Any}}(Inf)

	# start processing thread
	processingThreadRunning = true
	Threads.@spawn process_thread(fromGUI, toGUI)

	# setup gui
	w = Window(Dict(:width=>512,:height=>768))

	# event listeners
	handle(w, "msg") do args
		msgName = Symbol(args[1])
		msgArgs = args[2:end]
		@info "[GUI] sending message: $msgName $(join(msgArgs,", "))"
		processingThreadRunning && put!(fromGUI, msgName=>msgArgs)
	end

	doc = read(joinpath(@__DIR__,"content.html"),String)
	body!(w,doc,async=false)
	js(w, js"init()")

	while isopen(w.content.sock) # is there a better way to check if the window is still open?
		# @info "[GUI] tick"

		if isready(toGUI)
			msg = take!(toGUI)
			msgName = msg.first
			msgArgs = msg.second
			@info "[GUI] got message: $msgName"

			if msgName == :samplestatus
				js(w, js"""setSampleStatus($msgArgs)""")
			elseif msgName == :displayplot
				display(plot(msgArgs...))
			elseif msgName == :displaysampleannotnames
				js(w, js"""setSampleAnnotNames($msgArgs)""")
			elseif msgName == :exited
				processingThreadRunning = false
			end
		end
		# yield() # Allow GUI to run
		sleep(0.05) # Allow GUI to run
	end

	@info "[GUI] Waiting for scheduler thread to finish."
	processingThreadRunning && put!(fromGUI, :exit=>[])
	# wait until all threads have exited
	while processingThreadRunning
		msg = take!(toGUI)
		msgName = msg.first
		msgArgs = msg.second
		@info "[GUI] got message: $msgName"
		msgName == :exited && (processingThreadRunning = false)
		sleep(0.05)
	end
	@info "[GUI] Scheduler thread finished."
end


# struct SampleData
# 	data::Matrix
# 	sa::DataFrame
# 	va::DataFrame
# 	filepath::String
# 	errorMsg::String
# end
# SampleData() = SampleData(zeros(0,0),DataFrame(),DataFrame(),"","")
# SampleData(data::Matrix, sa::DataFrame, va::DataFrame,filepath::String) = SampleData(data,sa,va,filepath,"")
# SampleData(errorMsg::String) = SampleData(zeros(0,0),DataFrame(),DataFrame(),"",errorMsg)

# function opendataset(filepath::String)::SampleData
# 	isempty(filepath) && return SampleData("Please select file.")
# 	isfile(filepath) || return SampleData("File not found: \"$filepath\".")
# 	try
# 		data,sa,va = Qlucore.read(filepath)
# 		SampleData(data,sa,va,filepath)
# 	catch e
# 		SampleData("Error loading file: \"$filepath\" ($e)")
# 	end
# end


# struct Result
# 	U::Matrix{Float64}
# 	S::Vector{Float64}
# 	V::Matrix{Float64}
# 	G::AbstractMatrix # the graph (needed for plotting)
# end
# Result() = Result(zeros(0,0),zeros(0),zeros(0,0),zeros(0,0))

# function rundimreduction(ds::SampleData;
#                          dimReductionMethod::Symbol, sampleMethod::Symbol, sampleAnnotation::Symbol,
#                          timeAnnotation::Symbol, kNearestNeighbors::Int, distNearestNeighbors::Float64)
# 	isempty(ds.errorMsg) || return
# 	@assert dimReductionMethod in (:PMA,:PCA)
# 	@assert sampleMethod in (:SA,:Time,:NN,:NNSA)
# 	# @assert plotDims in 2:3

# 	X = zeros(size(ds.data))
# 	if any(ismissing,ds.data)
# 		# Replace missing values with mean over samples with nonmissing data
# 		println("Reconstructing missing values (taking the mean over all nonmissing samples)")
# 		for i=1:size(X,1)
# 			m = ismissing.(ds.data[i,:])
# 			X[i,.!m] .= ds.data[i,.!m]
# 			X[i,m] .= mean(ds.data[i,.!m])
# 		end
# 	else
# 		X .= ds.data # just copy
# 	end

# 	normalizemeanstd!(X)

# 	G = nothing
# 	if sampleMethod == :SA
# 		G = buildgraph(ds.sa[!,sampleAnnotation])
# 	elseif sampleMethod == :Time
# 		eltype(ds.sa[!,timeAnnotation])<:Number || @warn "Expected time annotation to contain numbers, got $(eltype(ds.sa[!,timeAnnotation])). Fallback to default sorting."
# 		G = buildgraph(ds.sa[!,sampleAnnotation], ds.sa[!,timeAnnotation])
# 	elseif sampleMethod == :NN
# 		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50);
# 	elseif sampleMethod == :NNSA
# 		G = neighborhoodgraph(X,kNearestNeighbors,distNearestNeighbors,50; groupBy=ds.sa[!,sampleAnnotation]);
# 	end

# 	# dim = 3
# 	dim = min(10, size(X)...)

# 	if dimReductionMethod==:PMA
# 		U,S,V = pma(X,G,dim=dim)
# 	elseif dimReductionMethod==:PCA
# 		U,S,V = pca(X,dim=dim)
# 	end
# 	Result(U,S,V,G)
# end

# function showplot(ds::SampleData, result::Result, plotDims::Int;
#                   dimReductionMethod::Symbol, sampleAnnotation::Symbol,
#                   kwargs...)
# 	colorBy = sampleAnnotation
# 	colorDict = colordict(ds.sa[!,colorBy])

# 	println(collect(keys(colorDict)))
# 	println(collect(values(colorDict)))

# 	opacity = 0.05
# 	title = string(splitdir(ds.filepath)[2], " ", string(dimReductionMethod))
# 	drawTriangles = false#true
# 	drawLines = true

# 	if plotDims==3
# 		markerSize = 5
# 		lineWidth = 1

# 		pl = plotsimplices(result.V,ds.sa,result.G,colorBy,colorDict, title=title,
# 		                   drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
# 		                   opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
# 		                   width=1024, height=768)
# 		display(pl)
# 	elseif plotDims==2
# 		markerSize = 0.9mm
# 		lineWidth = 0.3mm

# 		pl = plotsimplices_gadfly(result.V,ds.sa,result.G,colorBy,colorDict, title=title,
# 		                          drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
# 		                          opacity=opacity, markerSize=markerSize, lineWidth=lineWidth)
# 		display(pl)
# 	end
# end


# function exportloadings(ds::SampleData, result::Result, filepath::String, varAnnot::Symbol, columns::Symbol, dims::Int, columnSort::Symbol;
#                         dimReductionMethod::Symbol, sampleAnnotation::Symbol,
#                         kwargs...)
# 	@assert columns in (:All, :Selected)
# 	@assert columnSort in (:Abs, :Descending)

# 	df = DataFrame()
# 	df[!,varAnnot] = ds.va[!,varAnnot]

# 	if columns==:All
# 		for i=1:min(dims,length(result.S))
# 			colName = Symbol("$dimReductionMethod$i")
# 			df[!,colName] = result.U[:,i]
# 		end
# 	elseif columns==:Selected
# 		colName = Symbol("$dimReductionMethod$dims")
# 		df[!,colName] = result.U[:,dims]

# 		if columnSort==:Abs
# 			sort!(df, colName, by=abs, rev=true)
# 		elseif columnSort==:Descending
# 			sort!(df, colName, rev=true)
# 		end
# 	end

# 	# TODO: add CSV dependency and save using CSV.write() instead?
# 	mat = Matrix{String}(undef, size(df,1)+1, size(df,2))
# 	mat[1,:] .= string.(names(df))
# 	mat[2:end,:] .= string.(df)
# 	writedlm(filepath, mat, '\t')
# end


# function main()
# 	# init
# 	dataset = opendataset("")
# 	dimReductionParams = nothing
# 	result = nothing
# 	messageQueue = Queue{Pair{String,Any}}()


# 	# setup gui
# 	w = Window(Dict(:width=>512,:height=>768))

# 	# event listeners
# 	handle(w, "gedataopen") do args
# 		println("gedataopen: ", args)
# 		fn = isempty(args) ? "" : args[1]
# 		enqueue!(messageQueue, "gedataopen"=>fn)
# 	end
# 	handle(w, "samplemethod") do args
# 		println("samplemethod: ", args)
# 		enqueue!(messageQueue, "samplemethod"=>args)
# 	end
# 	handle(w, "showplot") do args
# 		println("showplot: ", args)
# 		enqueue!(messageQueue, "showplot"=>args)
# 	end
# 	handle(w, "exportloadings") do args
# 		println("exportloadings: ", args)
# 		enqueue!(messageQueue, "exportloadings"=>args)
# 	end

# 	doc = read(joinpath(@__DIR__,"content.html"),String)
# 	body!(w,doc,async=false)
# 	changesamplemethod!(w,:SA)


# 	# Message handling loop
# 	while isopen(w.content.sock) # is there a better way to check if the window is still open?
# 		if !isempty(messageQueue)
# 			msg = dequeue!(messageQueue)
# 			if msg.first == "gedataopen"
# 				filepath = msg.second
# 				println("Loading ", filepath)
# 				dataset = opendataset(filepath)
# 				dimReductionParams = nothing
# 				result = nothing
# 				if isempty(dataset.errorMsg)
# 					println(filepath, " loaded successfully")
# 				else
# 					println("Loading failed with error, ", dataset.errorMsg)
# 				end
# 				populategui!(w,dataset)
# 			elseif msg.first == "samplemethod"
# 				sampleMethod = Symbol(msg.second)
# 				println("Changing method to ", sampleMethod)
# 				changesamplemethod!(w,sampleMethod)
# 			elseif msg.first == "showplot"
# 				args = msg.second
# 				dimReductionMethod = Symbol(args[1])

# 				params = Dict()
# 				try
# 					params = Dict(:dimReductionMethod=>dimReductionMethod,
# 					              :sampleMethod=>Symbol(args[2]),
# 					              :sampleAnnotation=>Symbol(args[3]),
# 					              :timeAnnotation=>Symbol(args[4]),
# 					              :kNearestNeighbors=>parse(Int,args[5]),
# 					              :distNearestNeighbors=>parse(Float64,args[6]))
# 				catch e
# 					println("Failed to parse parameters ", args, ": ", sprint(showerror, e))
# 				end

# 				if params != dimReductionParams
# 					try
# 						println("Running ", dimReductionMethod)
# 						result = rundimreduction(dataset; params...)
# 						println("Done")
# 						dimReductionParams = params
# 					catch e
# 						println("Failed to run ", dimReductionMethod, ": ", sprint(showerror, e))
# 						result = nothing
# 						dimReductionParams = nothing
# 					end
# 				end

# 				if result!=nothing
# 					try
# 						println("Showing plot")
# 						plotDims = parse(Int,args[7])
# 						showplot(dataset, result, plotDims; dimReductionParams...)
# 					catch e
# 						println("Error showing plot: ", sprint(showerror,e))
# 					end
# 				end
# 			elseif msg.first == "exportloadings"
# 				# TODO: avoid code duplication
# 				args = msg.second
# 				dimReductionMethod = Symbol(args[1])

# 				params = Dict()
# 				try
# 					params = Dict(:dimReductionMethod=>dimReductionMethod,
# 					              :sampleMethod=>Symbol(args[2]),
# 					              :sampleAnnotation=>Symbol(args[3]),
# 					              :timeAnnotation=>Symbol(args[4]),
# 					              :kNearestNeighbors=>parse(Int,args[5]),
# 					              :distNearestNeighbors=>parse(Float64,args[6]))
# 				catch e
# 					println("Failed to parse parameters ", args, ": ", sprint(showerror, e))
# 				end

# 				if params != dimReductionParams
# 					try
# 						println("Running ", dimReductionMethod)
# 						result = rundimreduction(dataset; params...)
# 						println("Done")
# 						dimReductionParams = params
# 					catch e
# 						println("Failed to run ", dimReductionMethod, ": ", sprint(showerror, e))
# 						result = nothing
# 						dimReductionParams = nothing
# 					end
# 				end

# 				if result!=nothing
# 					try
# 						println("Exporting loadings")
# 						println(args)
# 						filepath = args[7]
# 						varAnnot = Symbol(args[8])
# 						columns = Symbol(args[9])
# 						dims = parse(Int,args[10])
# 						columnSort = Symbol(args[11])
# 						exportloadings(dataset, result, filepath, varAnnot, columns, dims, columnSort; dimReductionParams...)
# 						println("Done")
# 					catch e
# 						println("Error exporting loadings: ", sprint(showerror,e))
# 					end
# 				end
# 			else
# 				@warn "Unknown message type: $(msg.first)"
# 			end
# 		end

# 		# yield() # Allow GUI to run
# 		sleep(0.05) # Allow GUI to run
# 	end


# end


end

PMAGUI.main()
