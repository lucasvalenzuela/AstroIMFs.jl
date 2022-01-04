using AstroIMFs
using Documenter

DocMeta.setdocmeta!(AstroIMFs, :DocTestSetup, :(using AstroIMFs); recursive=true)

makedocs(;
    modules=[AstroIMFs],
    authors="Lucas Valenzuela",
    repo="https://github.com/lucasvalenzuela/AstroIMFs.jl/blob/{commit}{path}#{line}",
    sitename="AstroIMFs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lucasvalenzuela.github.io/AstroIMFs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lucasvalenzuela/AstroIMFs.jl",
    devbranch="main",
)
