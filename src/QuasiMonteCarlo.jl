module QuasiMonteCarlo

using Sobol, LatticeRules, Primes, LinearAlgebra, Random
using ConcreteStructs, Accessors

abstract type SamplingAlgorithm end
abstract type RandomSamplingAlgorithm <: SamplingAlgorithm end
abstract type DeterministicSamplingAlgorithm <: SamplingAlgorithm end
abstract type RandomizationMethod end

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
    RandomSample <: RandomSamplingAlgorithm

Randomly distributed random numbers.
"""
Base.@kwdef @concrete struct RandomSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function sample(n::Integer, d::Integer, S::RandomSample, T = Float64)
    _check_sequence(n)
    return rand(S.rng, T, d, n)
end

"""
```julia
sample(n::Integer, d::Integer, S::SamplingAlgorithm, T = Float64)
sample(n::Integer,
    lb::T,
    ub::T,
    S::SamplingAlgorithm) where {T <: Union{Base.AbstractVecOrTuple, Number}}
```

Return a QMC point set where:

  - `n` is the number of points to sample.
  - `S` is the quasi-Monte Carlo sampling strategy.
    The point set is in a `d`-dimensional unit box `[0, 1]^d`.
    If the bounds are specified, the sample is transformed (translation + scaling) into a box `[lb, ub]` where:
  - `lb` is the lower bound for each variable. Its length fixes the dimensionality of the sample.
  - `ub` is the upper bound. Its dimension must match `length(lb)`.

In the first method the type of the point set is specified by `T` while in the second method the output type is infered from the bound types.
"""
function sample(n::Integer, lb::T, ub::T,
    S::D) where {T <: Union{Base.AbstractVecOrTuple, Number},
    D <: SamplingAlgorithm}
    _check_sequence(lb, ub, n)
    lb = float.(lb)
    ub = float.(ub)
    out = sample(n, length(lb), S, eltype(lb))
    return (ub .- lb) .* out .+ lb
end

# See https://discourse.julialang.org/t/is-there-a-dedicated-function-computing-m-int-log-b-b-m/89776/10
function logi(b::Integer, n::Integer)
    m = round(Int, log(b, n))
    b^m == n || throw(ArgumentError("$n is not a power of $b"))
    return m
end

include("net_utilities.jl")
include("VanDerCorput.jl")
include("Faure.jl")
include("Kronecker.jl")
include("Halton.jl")
include("Sobol.jl")
include("LatinHypercube.jl")
include("Lattices.jl")
include("Section.jl")

"""
```julia
NoRand <: RandomizationMethod
```

No randomization is performed on the sampled sequence.
"""
struct NoRand <: RandomizationMethod end
randomize(x, S::NoRand) = x

include("RandomizedQuasiMonteCarlo/shifting.jl")
include("RandomizedQuasiMonteCarlo/scrambling_base_b.jl")
include("RandomizedQuasiMonteCarlo/iterators.jl")

import Requires
@static if !isdefined(Base, :get_extension)
    function __init__()
        Requires.@require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
            include("../ext/QuasiMonteCarloDistributionsExt.jl")
        end
    end
end

export SamplingAlgorithm,
    GridSample,
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
    randomize,
    RandomizationMethod,
    NoRand,
    Shift,
    ScrambleMethod,
    OwenScramble,
    MatousekScramble,
    DigitalShift,
    DesignMatrix
end # module
