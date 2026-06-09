using Pkg

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "Core" || GROUP == "All"
    include("core.jl")
end

if (GROUP == "QA" || GROUP == "All") && isempty(VERSION.prerelease)
    Pkg.activate(joinpath(@__DIR__, "qa"))
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("qa/qa.jl")
end
