using RvSpectML
using Documenter

makedocs(;
    modules=[RvSpectML],
    authors="Eric Ford and contributors",
    repo="https://github.com/RvSpectML/RvSpectML.jl/blob/{commit}{path}#L{line}",
    sitename="RvSpectML.jl",
#    latest="main",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/RvSpectML.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => [
            "Overview" => "contents.md",
            "Modules" => "modules.md",
            "Functions" => "functions.md",
            "Types" => "types.md"
            ],
        "Internals" => "internals.md",
        "Index" => "longlist.md",
    ],
    checkdocs=:none,
    #checkdocs=:exports,
)

deploydocs(;
    repo="github.com/RvSpectML/RvSpectML.jl",
)
