"""
    LatinHypercubeSample(rng::AbstractRNG = Random.TaskLocalRNG()) <: RandomSamplingAlgorithm

A Latin hypercube is a point set with the property that every one-dimensional
interval `(i / n, (i + 1) / n)` contains exactly one point. It is useful for
high-dimensional sampling because it is more uniform than independent Monte
Carlo sampling without requiring a full tensor-product grid.

# Fields

- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used to
  independently permute the strata in each dimension.

# Examples

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> points = sample(8, 2, LatinHypercubeSample(MersenneTwister(42)));

julia> size(points)
(2, 8)
```
"""
Base.@kwdef @concrete struct LatinHypercubeSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function sample(n::Integer, d::Integer, S::LatinHypercubeSample, T = Float64)
    _check_sequence(n)
    rng = S.rng
    seq = ((1:n) .- convert(T, 0.5)) / n
    # TODO: Probably a more efficient way to do this
    return reduce(vcat, [shuffle(rng, seq)' for _ in 1:d])
end
