"""
    SobolSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Samples taken from Sobol's base-2 sequence.

# Fields

- `R::RandomizationMethod = NoRand()`: Randomization applied to the Sobol
  points.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> points = sample(8, 2, SobolSample());

julia> size(points)
(2, 8)
```
"""
Base.@kwdef @concrete struct SobolSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::SobolSample, T = Float64)
    if n < 0
        throw(ArgumentError("number of samples must be non-negative"))
    end

    seq = Matrix{T}(undef, d, n)
    if n == 0
        return seq
    end

    # Use function barrier since `Sobol.SobolSeq(d)` can't be inferred
    return _sample!(seq, Sobol.SobolSeq(d), S.R)
end

function _sample!(seq::AbstractMatrix, s::Sobol.SobolSeq, R::RandomizationMethod)
    n = size(seq, 2)
    skip(s, n)
    for x in eachcol(seq)
        Sobol.next!(s, x)
    end
    return randomize(seq, R)
end
