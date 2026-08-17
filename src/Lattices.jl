"""
    GridSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

A simple rectangular grid lattice. It samples from a regular grid in the unit
box and optionally applies `R`.

# Fields

- `R::RandomizationMethod = NoRand()`: Randomization applied to the grid.

In more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a
result, they should almost never be used for multivariate integration; their use is as
a starting point for other algorithms.
"""
Base.@kwdef @concrete struct GridSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::GridSample, T = Float64)
    samples = rand.(range.(zeros(T, d), ones(T, d); length = n + 1), Ref(n))
    return randomize(mapreduce(permutedims, vcat, samples), S.R)
end

"""
    LatticeRuleSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Generate a point set using a rank-1 lattice rule.

# Fields

- `R::RandomizationMethod = NoRand()`: Randomization applied to the lattice
  points.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> points = sample(8, 2, LatticeRuleSample());

julia> size(points)
(2, 8)
```
"""
Base.@kwdef @concrete struct LatticeRuleSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::LatticeRuleSample, T = Float64)
    lat = LatticeRules.LatticeRule(d)
    result = reduce(hcat, lat[0:(n - 1)])
    return randomize(T == Float64 ? result : T.(result), S.R)
end
