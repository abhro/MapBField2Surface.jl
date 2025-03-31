using MapBField2Surface
using Documenter

DocMeta.setdocmeta!(MapBField2Surface, :DocTestSetup, :(using MapBField2Surface); recursive=true)

pages = ["Home" => "index.md"]

makedocs(;
    modules = [MapBField2Surface],
    authors = "Abhro R. and contributors",
    sitename = "MapBField2Surface.jl",
    format = Documenter.HTML(),
    pages,
)

deploydocs(;
    repo = "github.com/abhro/MapBField2Surface.jl",
    devbranch = "main",
)
