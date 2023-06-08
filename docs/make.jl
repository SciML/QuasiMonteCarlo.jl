using Documenter, QuasiMonteCarlo

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "QuasiMonteCarlo.jl",
    authors = "Chris Rackauckas",
    modules = [QuasiMonteCarlo],
    clean = true, doctest = false,
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/QuasiMonteCarlo/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/QuasiMonteCarlo.jl.git";
    push_preview = true)
