using Documenter, QuasiMonteCarlo

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "QuasiMonteCarlo.jl",
    authors = "Chris Rackauckas",
    modules = [QuasiMonteCarlo],
    clean = true, doctest = false, linkcheck = true,
    strict = [
        :doctest,
        :linkcheck,
        :parse_error,
        :example_block,
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
    ],
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/QuasiMonteCarlo/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/QuasiMonteCarlo.jl.git";
    push_preview = true)
