"""
    SobolSample <: SamplingAlgorithm

Samples taken from Sobol's base-2 sequence.
"""
struct SobolSample <: SamplingAlgorithm end

function sample(n::Integer, d::Integer, ::SobolSample, T = Float64)
    s = Sobol.SobolSeq(zeros(T, d), ones(T, d))
    skip(s, n)
    return reduce(hcat, [next!(s) for i in 1:n])
end
