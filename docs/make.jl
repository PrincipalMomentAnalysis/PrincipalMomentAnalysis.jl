push!(LOAD_PATH,"../src/")
using Documenter
using PrincipalMomentAnalysis

makedocs(
	sitename="PrincipalMomentAnalysis.jl",
	modules = [PrincipalMomentAnalysis],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
    	"Home" => "index.md",
    	"Tutorial" => "tutorial.md",
    	"Reference" => "reference.md",
    ]
)
