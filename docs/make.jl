using Documenter, PoissonRandom

include("pages.jl")

makedocs(
    sitename="PoissonRandom.jl",
    authors="Chris Rackauckas",
    modules=[PoissonRandom],
    clean=true,doctest=false,
    format = Documenter.HTML(analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://poissonrandom.sciml.ai/stable/"),
    pages=pages
)

deploydocs(
   repo = "github.com/SciML/PoissonRandom.jl.git";
   push_preview = true
)
