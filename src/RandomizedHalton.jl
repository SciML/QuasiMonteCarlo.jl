"""
    RandomizedHalton(rng::AbstractRNG = Random.GLOBAL_RNG) <: RandomSamplingAlgorithm

    A Halton-sequence randomized using the algorithm presented by A. Owen in "A randomized Halton algorithm in R" (2017).
"""
Base.@kwdef @concrete struct RandomizedHaltonSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function sample(n::Integer, d::Integer, S::RandomizedHaltonSample, T = Float64)
    _check_sequence(n)
    rng = S.rng
    bases = nextprimes(one(n), d)
    halton_seq = Matrix{T}(undef, d, n)

    for i in 1:d
        halton_seq[i, :] = randradinv(collect(1:n), S.rng, bases[i])
    end
    return halton_seq
end

function randradinv(ind::Vector{Int}, rng::AbstractRNG, b::Int = 2)
    b2r = 1 / b
    ans = ind .* 0
    res = ind

    while (1 - b2r < 1)
        dig = res .% b
        perm = shuffle(rng, collect(1:b)) .- 1
        pdig = perm[convert.(Int, dig .+ 1)]
        ans = ans .+ pdig .* b2r
        b2r = b2r / b
        res = (res .- dig) ./ b
    end
    return ans
end
