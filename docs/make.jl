using LastHomework
using Documenter

DocMeta.setdocmeta!(LastHomework, :DocTestSetup, :(using LastHomework); recursive=true)

makedocs(;
    modules=[LastHomework],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/singularitti/LastHomework.jl/blob/{commit}{path}#{line}",
    sitename="LastHomework.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/LastHomework.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/LastHomework.jl",
    devbranch="main",
)
