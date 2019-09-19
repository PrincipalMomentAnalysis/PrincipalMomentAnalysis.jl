module Qlucore

using DataFrames
using Missings


function parsefloatmissing(x::AbstractString)::Union{Missing,Float64}
	isempty(x) && return missing
	parse(Float64,x)
end


function read(filename::String)
	lowercase(splitext(filename)[2])==".gedata" || @warn "Expected $filename to be a .gedata file, parsing might fail."
	open(filename) do io
		read(io)
	end
end

function read(io::IO)
	line = readline(io)
	m = match(r"^qlucore\tgedata\tversion (.+)$", line)
	@assert m!=nothing "Invalid header line: $line"
	version = VersionNumber(m.captures[1])
	expectedVersions = (v"1.0",v"1.1")
	version in expectedVersions || @warn "Expected .gedata versions $expectedVersion, got $version, parsing might fail."

	nbrExpectedEmptyLines = version>=v"1.1" ? 2 : 1
	for i=1:nbrExpectedEmptyLines
		line = readline(io)
		@assert isempty(line) "Expected empty line, got $line"
	end

	# sample info
	line = readline(io)
	m = match(r"^(\d+)\tsamples\twith\t(\d+)\tannotations$", line)
	@assert m!=nothing "Invalid sample info line: $line"
	nbrSamples = parse(Int, m.captures[1])
	nbrSampleAnnotations = parse(Int, m.captures[2])

	# variable info
	line = readline(io)
	m = match(r"^(\d+)\tvariables\twith\t(\d+)\tannotations$", line)
	@assert m!=nothing "Invalid variables info line: $line"
	nbrVariables = parse(Int, m.captures[1])
	nbrVariableAnnotations = parse(Int, m.captures[2])

	# one empty lines
	line = readline(io)
	@assert isempty(line) "Expected empty line, got $line"

	# sample annotations
	sa = DataFrame()
	saNamesFound = Set{Symbol}()
	for i=1:nbrSampleAnnotations
		line = readline(io)
		splitLine = split(line,'\t')
		@assert length(splitLine)==nbrVariableAnnotations+1+nbrSamples "Could not parse sample annotations. Wrong number of columns."
		@assert all(i->isempty(splitLine[i]),1:nbrVariableAnnotations) "Could not parse sample annotations. Unexpected data before sample annotation."

		saName = Symbol(splitLine[nbrVariableAnnotations+1])
		@assert !in(saName, saNamesFound) "Could not parse sample annotations. Annotation \"$saName\" found twice."
		push!(saNamesFound,saName)
		sa[!,saName] = string.(splitLine[nbrVariableAnnotations+2:end])
	end

	# variable annotations and data
	data = Matrix{Union{Missing,Float64}}(missing,nbrVariables,nbrSamples) # consider parsing into transposed data matrix and transpose before returning
	va = DataFrame()

	# variable annotations header
	line = readline(io)
	splitLine = split(line,'\t')
	@assert length(splitLine)==nbrVariableAnnotations "Could not parse variable annotations header. Wrong number of columns."
	vaNames = Symbol.(splitLine)
	vaNamesFound = Set{Symbol}()

	# init variable annotations
	for vaName in vaNames
		@assert !in(vaName, vaNamesFound) "Could not parse variable annotations. Annotation \"$vaName\" found twice."
		push!(vaNamesFound,vaName)
		va[!,vaName] = Vector{String}(undef,nbrVariables)
	end

	# for each line with variable annotations and data
	# for i=1:nbrVariables
	for (i,vaRow) in enumerate(eachrow(va))
		line = readline(io)
		splitLine = split(line,'\t')
		@assert length(splitLine)==nbrVariableAnnotations+1+nbrSamples "Could not parse data line. Wrong number of columns."

		vaRow .= splitLine[1:nbrVariableAnnotations]
		@assert isempty(splitLine[nbrVariableAnnotations+1]) "Could not parse data line. Expected empty column between variable annotations and data."

		data[i,:] .= parsefloatmissing.(splitLine[nbrVariableAnnotations+2:end])
	end

	data,sa,va,version
end



end