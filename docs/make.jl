using AstrodynamicsUtilities
using Documenter

DocMeta.setdocmeta!(AstrodynamicsUtilities, :DocTestSetup, :(using AstrodynamicsUtilities); recursive=true)

makedocs(;
    modules=[AstrodynamicsUtilities],
    authors="Dylan LaSalle",
    sitename="AstrodynamicsUtilities.jl",
    format=Documenter.HTML(;
        canonical="https://dlasalle37.github.io/AstrodynamicsUtilities.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dlasalle37/AstrodynamicsUtilities.jl",
    devbranch="master",
)
