using PoissonSchrodingerSimulator
using Documenter

DocMeta.setdocmeta!(PoissonSchrodingerSimulator, :DocTestSetup, :(using PoissonSchrodingerSimulator); recursive=true)

makedocs(;
    modules=[PoissonSchrodingerSimulator],
    authors="singularitti <singularitti@outlook.com> and contributors",
    sitename="PoissonSchrodingerSimulator.jl",
    format=Documenter.HTML(;
        canonical="https://singularitti.github.io/PoissonSchrodingerSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/PoissonSchrodingerSimulator.jl",
    devbranch="main",
)
