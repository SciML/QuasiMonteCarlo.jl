module QuasiMonteCarlo

using Sobol: Sobol
using LatticeRules: LatticeRules
using Primes: nextprime, nextprimes
using LinearAlgebra: Diagonal, LowerTriangular, UpperTriangular
using Random: Random, AbstractRNG, rand!, shuffle
using ConcreteStructs: @concrete
using Accessors: @set

"""
    SamplingAlgorithm

Abstract supertype for sampling strategies accepted by [`sample`](@ref).

## Extension rules

Define a concrete subtype and implement
`sample(n::Integer, d::Integer, sampler::YourSampler, T = Float64)`.
The method must return a `d`-by-`n` matrix of elements of type `T` whose columns
are points in the unit box `[0, 1]^d`. It must reject invalid sampler-specific
inputs with an informative exception and must preserve the requested number of
points and dimensions. `sample` applies user-supplied bounds itself, so extensions
should implement the unit-box method rather than a separate bounds method.
"""
abstract type SamplingAlgorithm end

"""
    RandomSamplingAlgorithm <: SamplingAlgorithm

Abstract supertype for samplers that use an RNG to generate randomized point sets.
"""
abstract type RandomSamplingAlgorithm <: SamplingAlgorithm end

"""
    DeterministicSamplingAlgorithm <: SamplingAlgorithm

Abstract supertype for deterministic low-discrepancy samplers.
"""
abstract type DeterministicSamplingAlgorithm <: SamplingAlgorithm end

"""
    RandomizationMethod

Abstract supertype for methods that randomize deterministic point sets with [`randomize`](@ref).
"""
abstract type RandomizationMethod end

const UB_LB_MESSAGE = "Lower bound exceeds upper bound (lb > ub)"
const ZERO_SAMPLES_MESSAGE = "Number of samples must be greater than zero"
const DIM_MISMATCH_MESSAGE = "Dimensionality of lb and ub must match"

function _check_sequence(lb, ub, n::Integer)
    @assert length(lb) == length(ub) DIM_MISMATCH_MESSAGE
    @assert all(x -> x[1] ≤ x[2], zip(lb, ub)) UB_LB_MESSAGE
    return @assert n > 0 ZERO_SAMPLES_MESSAGE
end

_check_sequence(n::Integer) = @assert n > 0 ZERO_SAMPLES_MESSAGE

"""
    RandomSample(; rng = Random.TaskLocalRNG()) <: RandomSamplingAlgorithm

Plain Monte Carlo sampler that draws independent uniform samples from the unit box.
"""
Base.@kwdef @concrete struct RandomSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function sample(n::Integer, d::Integer, S::RandomSample, T = Float64)
    _check_sequence(n)
    return rand(S.rng, T, d, n)
end

"""
    sample(n::Integer, d::Integer, sampler::SamplingAlgorithm, T = Float64)
    sample(n::Integer, lb, ub, sampler::SamplingAlgorithm)

Generate `n` sample points with `sampler`.

## Arguments

  - `n`: Positive number of points to generate.
  - `d`: Positive dimension of the unit box. This form returns a `d`-by-`n`
    matrix with values in `[0, 1]`.
  - `lb`: Scalar, tuple, or vector of lower bounds. Its length determines the
    dimension when it is not scalar.
  - `ub`: Scalar, tuple, or vector of upper bounds with the same shape as `lb`.
    Each lower bound must be less than or equal to its corresponding upper bound.
  - `sampler`: Concrete [`SamplingAlgorithm`](@ref) that determines the point
    set construction.

## Optional Positional Arguments

  - `T = Float64`: Element type of the unit-box result. The bounds form infers
    its output element type from `lb`.

## Returns

A matrix whose columns are the generated points. The unit-box form has size
`(d, n)`. The bounds form has size `(length(lb), n)` for collection bounds and
maps every coordinate from `[0, 1]` to its corresponding closed interval
`[lb[i], ub[i]]`.

## Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> unit_points = sample(4, 2, SobolSample());

julia> size(unit_points)
(2, 4)

julia> points = sample(4, [0.0, -1.0], [1.0, 1.0], SobolSample());

julia> all([0.0, -1.0] .<= points .<= [1.0, 1.0])
true
```

## Extension rules

To add a sampler, subtype [`SamplingAlgorithm`](@ref) and implement only the
unit-box form `sample(n, d, sampler, T)`. The implementation must return a
`d`-by-`n` matrix with elements in `[0, 1]`; the bounds form is provided by this
package and delegates to that unit-box method. Do not extend this bounds method
for a new sampler.
"""
function sample(
        n::Integer, lb::T, ub::T,
        S::D
    ) where {
        T <: Union{AbstractVector, Tuple, Number},
        D <: SamplingAlgorithm,
    }
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
include("RandomizedHalton.jl")
include("Lattices.jl")
include("Section.jl")

"""
    NoRand() <: RandomizationMethod

No randomization is performed on the sampled sequence.
"""
struct NoRand <: RandomizationMethod end

"""
    randomize(x, R::RandomizationMethod)

Apply the randomization method `R` to the sample matrix `x`.
The `NoRand()` method returns `x` unchanged.
"""
randomize(x, S::NoRand) = x

include("RandomizedQuasiMonteCarlo/shifting.jl")
include("RandomizedQuasiMonteCarlo/scrambling_base_b.jl")
include("RandomizedQuasiMonteCarlo/iterators.jl")

include("precompile.jl")

export SamplingAlgorithm,
    sample,
    generate_design_matrices,
    GridSample,
    SobolSample,
    LatinHypercubeSample,
    RandomizedHaltonSample,
    LatticeRuleSample,
    RandomSample,
    HaltonSample,
    VanDerCorputSample,
    GoldenSample,
    KroneckerSample,
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
