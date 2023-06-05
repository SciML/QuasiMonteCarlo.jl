"""
    SobolSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Samples taken from Sobol's base-2 sequence.
"""
Base.@kwdef @concrete struct SobolSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::SobolSample, T = Float64)
    s = Sobol.SobolSeq(zeros(T, d), ones(T, d))
    skip(s, n)
    return randomize(reduce(hcat, [next!(s) for i in 1:n]), S.R)
end
