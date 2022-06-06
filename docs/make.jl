using Documenter, QuasiMonteCarlo

include("pages.jl")

makedocs(
    sitename="QuasiMonteCarlo.jl",
    authors="Chris Rackauckas",
    modules=[PoissonRandom],
    clean=true,doctest=false,
    format = Documenter.HTML(analytics = "UA-90474609-3",
                             assets = ["assets/favicon.ico"],
                             canonical="https://quasimontecarlo.sciml.ai/stable/"),
    pages=pages
)

deploydocs(
   repo = "github.com/SciML/QuasiMonteCarlo.jl.git";
   push_preview = true
)
