using Backboner
using Documenter

DocMeta.setdocmeta!(Backboner, :DocTestSetup, :(using VectorizedKmers); recursive=true)

makedocs(;
    modules = [Backboner],
    sitename = "Backboner.jl",
    doctest = false,
    pages = [
        "Home" => "index.md",
    ],
    authors = "Anton Oresten",
    checkdocs = :all,
)

deploydocs(;
    repo = "github.com/MurrellGroup/Backboner.jl.git",
    push_preview = true,
    deps = nothing,
    make = nothing,
)