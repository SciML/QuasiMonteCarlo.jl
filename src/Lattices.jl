"""
    GridSample()

A simple rectangular grid lattice.

In more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a
result, they should almost never be used for multivariate integration; their use is as
a starting point for other algorithms.
"""
Base.@kwdef @concrete struct GridSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::GridSample, T = Float64)
    n = convert(T, n)
    return randomize([(i - convert(T, 0.5)) / n for _ in 1:d, i in 1:n], S.R)
end

"""
    LatticeRuleSample()

Generate a point set using a lattice rule.
"""
Base.@kwdef @concrete struct LatticeRuleSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::LatticeRuleSample, T = Float64)
    lat = LatticeRules.LatticeRule(d)
    return randomize(reduce(hcat, lat[0:(n - 1)]), S.R)
end
