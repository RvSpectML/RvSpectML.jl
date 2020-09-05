using RvSpectML
using Documenter

makedocs(;
    modules=[RvSpectML],
    authors="Eric Ford and contributors",
    repo="https://github.com/eford/RvSpectML.jl/blob/{commit}{path}#L{line}",
    sitename="RvSpectML.jl",
#    latest="main",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://eford.github.io/RvSpectML.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Internals" => "internals.md",
        "Index" => "longlist.md",
    ],
)

deploydocs(;
    repo="github.com/eford/RvSpectML.jl",
    devbranch="main" 
)
