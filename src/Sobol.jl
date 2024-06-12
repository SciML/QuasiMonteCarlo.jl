"""
    SobolSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Samples taken from Sobol's base-2 sequence.
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
    Sobol.skip!(s, n, @view(seq[:, begin]))
    for x in eachcol(seq)
        Sobol.next!(s, x)
    end
    return randomize(seq, R)
end
