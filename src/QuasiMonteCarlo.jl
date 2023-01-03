module QuasiMonteCarlo

using Sobol, LatticeRules, Distributions, Primes, LinearAlgebra, Random
using ConcreteStructs

abstract type SamplingAlgorithm end

include("VanDerCorput.jl")
include("Faure.jl")
include("Kronecker.jl")
include("Halton.jl")
include("Sobol.jl")
include("LatinHypercube.jl")
include("Lattices.jl")
include("Section.jl")

const UB_LB_MESSAGE = "Lower bound exceeds upper bound (lb > ub)"
const ZERO_SAMPLES_MESSAGE = "Number of samples must be greater than zero"
const DIM_MISMATCH_MESSAGE = "Dimensionality of lb and ub must match"

function _check_sequence(lb, ub, n::Integer)
    @assert length(lb)==length(ub) DIM_MISMATCH_MESSAGE
    @assert all(x -> x[1] ≤ x[2], zip(lb, ub)) UB_LB_MESSAGE
    @assert n>0 ZERO_SAMPLES_MESSAGE
end

_check_sequence(n::Integer) = @assert n>0 ZERO_SAMPLES_MESSAGE

"""
    RandomSample <: SamplingAlgorithm

Randomly distributed random numbers.
"""
Base.@kwdef @concrete struct RandomSample <: SamplingAlgorithm
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function sample(n::Integer, d::Integer, S::RandomSample, T = Float64)
    _check_sequence(n)
    return rand(S.rng, T, d, n)
end

"""
    sample(n, lb, ub, sample_method)
    sample(n, d, sample_method)

Create a QMC point set, where:
- `n` is the number of points to sample.
- `lb` is the lower bound for each variable. The length determines the dimensionality.
- `ub` is the upper bound.
- `d` is the dimensionality of the point set.
- `sample_method` is the quasi-Monte Carlo sampling strategy. Note that any Distributions.jl
  distribution can be used in addition to any `SamplingAlgorithm`.
"""
function sample(n::Integer, lb::T, ub::T,
                S::SamplingAlgorithm) where {T <: Union{Base.AbstractVecOrTuple, Number}}
    _check_sequence(lb, ub, n)
    lb = float.(lb)
    ub = float.(ub)
    out = sample(n, length(lb), S, eltype(lb))
    return (ub .- lb) .* out .+ lb
end

"""
    sample(n, d, D::Distribution)

Returns a tuple containing independent random samples from distribution `D`.
"""
function sample(n::Integer, d::Integer, D::Distributions.Sampleable)
    @assert n>0 ZERO_SAMPLES_MESSAGE
    x = [[rand(D) for j in 1:d] for i in 1:n]
    return reduce(hcat, x)
end

function generate_design_matrices(n, lb, ub, sampler, num_mats = 2)
    if n <= 0
        throw(ZeroSamplesError())
    end
    @assert length(lb) == length(ub)
    d = length(lb)
    out = sample(n, repeat(lb, num_mats), repeat(ub, num_mats), sampler)
    [out[(j * d + 1):((j + 1) * d), :] for j in 0:(num_mats - 1)]
end

# See https://discourse.julialang.org/t/is-there-a-dedicated-function-computing-m-int-log-b-b-m/89776/10
function logi(b::Int, n::Int)
    m = round(Int, log(b, n))
    b^m == n || throw(ArgumentError("$n is not a power of $b"))
    return m
end

abstract type RandomizationMethod end
"""
```julia
NoRand
```

No Randomization is performed on the sampled sequence.
"""
struct NoRand <: RandomizationMethod end
randomization(x::AbstractArray, S::NoRand) = x

include("RandomizedQuasiMonteCarlo/shifting.jl")
include("RandomizedQuasiMonteCarlo/conversion.jl")
include("RandomizedQuasiMonteCarlo/scrambling_base_b.jl")

export SamplingAlgorithm,
       GridSample,
       UniformSample,
       SobolSample,
       LatinHypercubeSample,
       LatticeRuleSample,
       RandomSample,
       HaltonSample,
       VanDerCorputSample,
       GoldenSample,
       KroneckerSample,
       SectionSample,
       FaureSample,
       randomization,
       RandomizationMethod,
       NoRand,
       Shift,
       ScrambleMethod,
       OwenScramble,
       MatousekScramble,
       DigitalShift
end # module
