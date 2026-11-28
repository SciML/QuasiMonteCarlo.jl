using Documenter, QuasiMonteCarlo

# Doctests of internal functions name the package, so load it before each one.
DocMeta.setdocmeta!(QuasiMonteCarlo, :DocTestSetup, :(using QuasiMonteCarlo); recursive = true)

include("pages.jl")

makedocs(
    sitename = "QuasiMonteCarlo.jl",
    authors = "Chris Rackauckas",
    modules = [QuasiMonteCarlo],
    clean = true, doctest = true, linkcheck = true,
    checkdocs = :exports,
    linkcheck_ignore = [
        "https://www.sciencedirect.com/science/article/pii/S0010465509003087",
        "https://artowen.su.domains/mc/qmcstuff.pdf",
    ],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/QuasiMonteCarlo/stable/"
    ),
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/QuasiMonteCarlo.jl.git";
    push_preview = true
)
