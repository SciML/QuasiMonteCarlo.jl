"""
    RandomizedHaltonSample(; rng = Random.TaskLocalRNG()) <: RandomSamplingAlgorithm

Create a randomized Halton sequence by independently permuting the radical
inverse digits in each dimension.

# Fields

- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used for
  the digit permutations.

# Examples

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> points = sample(8, 2, RandomizedHaltonSample(rng = MersenneTwister(42)));

julia> size(points)
(2, 8)
```

References:
Owen, A. (2017). *A randomized Halton algorithm in R*. https://doi.org/10.48550/arXiv.1706.02808
"""
Base.@kwdef @concrete struct RandomizedHaltonSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function sample(n::Integer, d::Integer, S::RandomizedHaltonSample, T = Float64)
    _check_sequence(n)
    bases = nextprimes(one(n), d)
    halton_seq = Matrix{T}(undef, d, n)

    ind = collect(1:n)
    for i in 1:d
        halton_seq[i, :] = randradinv(ind, S.rng, bases[i])
    end
    return halton_seq
end

function randradinv(ind::Vector{Int}, rng::AbstractRNG, b::Int = 2)
    b2r = 1 / b
    ans = ind .* 0
    res = ind

    base_ind = 1:b

    while (1 - b2r < 1)
        dig = res .% b
        perm = shuffle(rng, base_ind) .- 1
        pdig = perm[convert.(Int, dig .+ 1)]
        ans = ans .+ pdig .* b2r
        b2r = b2r / b
        res = (res .- dig) ./ b
    end
    return ans
end
