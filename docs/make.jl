using Lanczos
using Documenter

DocMeta.setdocmeta!(Lanczos, :DocTestSetup, :(using Lanczos); recursive=true)

makedocs(;
    modules=[Lanczos],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="Lanczos.jl",
    format=Documenter.HTML(;
        canonical="https://singularitti.github.io/Lanczos.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/Lanczos.jl",
    devbranch="main",
)
