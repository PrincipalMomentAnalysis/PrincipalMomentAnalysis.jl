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
	js(w, js"""document.getElementById("sampleannot").options.length = 0;""")
	for a in sa
		js(w, js"""var opt = document.createElement("option");
			       opt.value = $a;
			       opt.text = $a;
			       document.getElementById("sampleannot").options.add(opt);""")
	end

	# js(w, js"""console.log('$infoStr');""")	
end


function runpma(ds::Dataset, sampleAnnotation::String)
	isempty(ds.errorMsg) || return

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
	                      opacity=opacity, markerSize=markerSize, lineWidth=lineWidth)
	display(plPMA)
end



function main()
	# init
	dataset = Dataset()
	messageQueue = Queue{Pair{String,Any}}()


	# setup gui
	w = Window(Dict(:width=>512,:height=>384))

	# doc = """<button onclick='Blink.msg("gedataopen", "hej")'>Load Qlucore .gedata</button>"""
	# doc  = """<button onclick="var electron = require('electron'); var fp = electron.remote.dialog.showOpenDialog({properties: ['openFile']}); if (!(fp===undefined)) { Blink.msg('gedataopen',fp) }">Open Qlucore .gedata file</button>"""
	doc  = """<div>
	              <button onclick='var electron = require("electron"); var fp = electron.remote.dialog.showOpenDialog({properties: ["openFile"]}); if (!(fp===undefined)) { Blink.msg("gedataopen",fp) }'>Open Qlucore .gedata file</button>
	              <p id="info">Please select file.<p/>
	              <div>
	                  <p>Choose sample annotation:</p>
	                  <select id="sampleannot"></select>
	              </div>
	              <button onclick='Blink.msg("runpma", document.getElementById("sampleannot").value)'>Run PMA</button>
	          </div>
	       """

	# event listeners
	handle(w, "gedataopen") do args
		println("gedataopen: ", args)
		fn = isempty(args) ? "" : args[1]
		enqueue!(messageQueue, "gedataopen"=>fn)
	end

	# event listeners
	handle(w, "runpma") do args
		println("runpma: ", args)
		enqueue!(messageQueue, "runpma"=>args)
	end

	# handle(w, "closed") do args...
	# 	println("Closed!")
	# end

	body!(w,doc,async=false)

	
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
			elseif msg.first == "runpma"
				println("Running PMA ")
				runpma(dataset,msg.second)
			else
				@warn "Unknown message type: $(msg.first)"
			end
		end

		yield() # Allow GUI to run
	end


end


end

PMAGUI.main()