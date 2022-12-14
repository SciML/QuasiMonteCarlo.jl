"""
    LatinHypercubeSample <: SamplingAlgorithm

A Latin Hypercube is a point set with the property that every one-dimensional interval `(i/n, i+1/n)` contains exactly one point. This is a good way to sample a high-dimensional space, as it is more uniform than a random sample but does not require as many points as a full net.
"""
Base.@kwdef @concrete struct LatinHypercubeSample <: SamplingAlgorithm
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function sample(n::Integer, d::Integer, S::LatinHypercubeSample, T=Float64)
    _check_sequence(n)
    rng = S.rng
    seq = ((1:n) .- convert(T, 0.5)) / n
    # TODO: Probably a more efficient way to do this
    return reduce(vcat, [shuffle(rng, seq)' for _ in 1:d])
end
