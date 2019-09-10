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

include("qlucore.jl")
include("pmaplots.jl")
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


function changesamplemethod!(w, sampleMethod::String)
	enabled = sampleMethod in ("SA","Time","NNSA")
	js(w, js"""document.getElementById("sampleAnnot").disabled = $(!enabled);""")

	enabled = sampleMethod == "Time"
	js(w, js"""document.getElementById("timeAnnot").disabled = $(!enabled);""")

	enabled = sampleMethod in ("NN","NNSA")
	js(w, js"""document.getElementById("kNearestNeighbors").disabled = $(!enabled);""")
	js(w, js"""document.getElementById("distNearestNeighbors").disabled = $(!enabled);""")
end



function runpma(ds::Dataset, sampleMethod::String, sampleAnnotation::String, timeAnnotation::String, kNearestNeighbors::String, distNearestNeighbors::String)
	isempty(ds.errorMsg) || return
	@assert sampleMethod in ("SA","Time","NN","NNSA")

	sampleAnnotation = Symbol(sampleAnnotation)

	# TODO: handle missing values and avoid making too many matrix copies
	X = copy(ds.data)
	X = convert(Matrix{Float64},X)
	normalizemeanstd!(X)

	G = buildgraph(ds.sa[!,sampleAnnotation])
	dim = 3

	UPMA,SPMA,VPMA = pma(X,G,dim=dim)

	colorBy = sampleAnnotation
	colorDict = colordict(ds.sa[!,colorBy])

	println(collect(keys(colorDict)))
	println(collect(values(colorDict)))

	markerSize = 5
	lineWidth = 1
	opacity = 0.05
	title = splitdir(ds.filepath)[2]
	drawTriangles = false#true
	drawLines = true

	plPMA = plotsimplices(VPMA,ds.sa,G,colorBy,colorDict, title="$title PMA",
	                      drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
	                      opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
	                      width=1024, height=768)
	display(plPMA)
end

function runpca(ds::Dataset, sampleMethod::String, sampleAnnotation::String, timeAnnotation::String, kNearestNeighbors::String, distNearestNeighbors::String)
	isempty(ds.errorMsg) || return
	@assert sampleMethod in ("SA","Time","NN","NNSA")

	sampleAnnotation = Symbol(sampleAnnotation)

	# TODO: handle missing values and avoid making too many matrix copies
	X = copy(ds.data)
	X = convert(Matrix{Float64},X)
	normalizemeanstd!(X)

	G = buildgraph(ds.sa[!,sampleAnnotation])
	dim = 3

	UPCA,SPCA,VPCA = pca(X,dim=dim)

	colorBy = sampleAnnotation
	colorDict = colordict(ds.sa[!,colorBy])

	println(collect(keys(colorDict)))
	println(collect(values(colorDict)))

	markerSize = 5
	lineWidth = 1
	opacity = 0.05
	title = splitdir(ds.filepath)[2]
	drawTriangles = false#true
	drawLines = true

	plPCA = plotsimplices(VPCA,ds.sa,G,colorBy,colorDict, title="$title PCA",
	                      drawTriangles=drawTriangles, drawLines=drawLines, drawPoints=true,
	                      opacity=opacity, markerSize=markerSize, lineWidth=lineWidth,
	                      width=1024, height=768)
	display(plPCA)
end


function main()
	# init
	dataset = Dataset()
	messageQueue = Queue{Pair{String,Any}}()


	# setup gui
	w = Window(Dict(:width=>512,:height=>384))

	# doc = """<button onclick='Blink.msg("gedataopen", "hej")'>Load Qlucore .gedata</button>"""
	# doc  = """<button onclick="var electron = require('electron'); var fp = electron.remote.dialog.showOpenDialog({properties: ['openFile']}); if (!(fp===undefined)) { Blink.msg('gedataopen',fp) }">Open Qlucore .gedata file</button>"""
	# doc  = """<div>
	#               <button onclick='var electron = require("electron"); var fp = electron.remote.dialog.showOpenDialog({properties: ["openFile"]}); if (!(fp===undefined)) { Blink.msg("gedataopen",fp) }'>Open Qlucore .gedata file</button>
	#               <p id="info">Please select file.<p/>
	#               <div>
	#                   <p>Choose sample annotation:</p>
	#                   <select id="sampleAnnot"></select>
	#                   <p>Choose time annotation:</p>
	#                   <select id="timeAnnot"></select>
	#                   <p>Number of Nearest Neighbors:</p>
	#                   <input type="number" id="kNearestNeighbors" min="0" max="1000000">
	#                   <p>Nearest Neighbors distance:</p>
	#                   <input type="range" id="distNearestNeighbors" min="0" max="1" step="0.001" disabled>
	#               </div>
	#               <button onclick='Blink.msg("runpma", document.getElementById("sampleAnnot").value)'>Run PMA</button>
	#           </div>
	#        """
	doc  = """<div>
	              <button onclick='var electron = require("electron"); var fp = electron.remote.dialog.showOpenDialog({properties: ["openFile"]}); if (!(fp===undefined)) { Blink.msg("gedataopen",fp) }'>Open Qlucore .gedata file</button>
	              <p id="info">Please select file.<p/>
	              <div>
	                  <p>Choose Method:</p>
	                  <div>
	                      <input type="radio" name="samplemethod" value="methodSA"   onclick='Blink.msg("samplemethod","SA")'   checked>Sample Annotation<br>
	                      <input type="radio" name="samplemethod" value="methodTime" onclick='Blink.msg("samplemethod","Time")'        >Time Series<br>
	                      <input type="radio" name="samplemethod" value="methodNN"   onclick='Blink.msg("samplemethod","NN")'          >Nearest Neighbor<br>
	                      <input type="radio" name="samplemethod" value="methodNNSA" onclick='Blink.msg("samplemethod","NNSA")'        >Nearest Neighbor within groups<br>
	                  </div>
	                  <p>Choose sample annotation:</p>
	                  <select id="sampleAnnot"></select>
	                  <p>Choose time annotation:</p>
	                  <select id="timeAnnot"></select>
	                  <p>Number of Nearest Neighbors:</p>
	                  <input type="number" id="kNearestNeighbors" min="0" max="1000000" value="0">
	                  <p>Nearest Neighbors distance:</p>
	                  <input type="range" id="distNearestNeighbors" min="0" max="1" step="0.001" value="0">
	              </div>
	              <button onclick='Blink.msg("runpma", ["SA", document.getElementById("sampleAnnot").value, document.getElementById("timeAnnot").value, document.getElementById("kNearestNeighbors").value, document.getElementById("distNearestNeighbors").value])'>Run PMA</button>
	              <button onclick='Blink.msg("runpca", ["SA", document.getElementById("sampleAnnot").value, document.getElementById("timeAnnot").value, document.getElementById("kNearestNeighbors").value, document.getElementById("distNearestNeighbors").value])'>Run PCA</button>
	          </div>
	       """

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

	# handle(w, "closed") do args...
	# 	println("Closed!")
	# end

	body!(w,doc,async=false)
	changesamplemethod!(w,"SA")

	
	# Message handling loop
	# while true
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
				sampleMethod = msg.second
				println("Changing method to ", sampleMethod)
				changesamplemethod!(w,sampleMethod)
			elseif msg.first == "runpma"
				println("Running PMA ")
				runpma(dataset,msg.second...)
			elseif msg.first == "runpca"
				println("Running PCA ")
				runpca(dataset,msg.second...)
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