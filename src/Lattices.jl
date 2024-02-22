"""
    GridSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

A simple rectangular grid lattice. Samples `n` random samples from a grid with `dx = (ub -lb)/n``

In more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a
result, they should almost never be used for multivariate integration; their use is as
a starting point for other algorithms.
"""
Base.@kwdef @concrete struct GridSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::GridSample, T = Float64)
    samples = rand.(range.(zeros(T, d), ones(T, d); length = n + 1), Ref(n))
    randomize(mapreduce(permutedims, vcat, samples), S.R)
end

"""
    LatticeRuleSample() <: DeterministicSamplingAlgorithm

Generate a point set using a lattice rule.
"""
Base.@kwdef @concrete struct LatticeRuleSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::LatticeRuleSample, T = Float64)
    lat = LatticeRules.LatticeRule(d)
    return randomize(reduce(hcat, lat[0:(n - 1)]), S.R)
end
