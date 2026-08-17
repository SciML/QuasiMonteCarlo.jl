"""
    AbstractDesignMatrix

Developer-facing interface for iterators returned by [`DesignMatrix`](@ref).

Concrete subtypes must store a `count` field and implement
`next!(iterator)`, which produces the next point-set matrix. The generic
`Base.length` and `Base.iterate` methods use those two pieces of the contract:
iteration yields exactly `count` matrices, and each call to `next!` may reuse
internal storage but must return a complete matrix before the next iteration
step mutates that storage.

Implementations should also define `Base.eltype(::Type{<:AbstractDesignMatrix})`
to describe the matrix type yielded by iteration. This interface is intended
for package developers extending the design-matrix machinery; application code
should call [`DesignMatrix`](@ref) and iterate over its result.
"""
abstract type AbstractDesignMatrix end
@public AbstractDesignMatrix

# Make an iterator so that we can do "for X in DesignMatrix(...)"
Base.length(s::AbstractDesignMatrix) = s.count
function Base.iterate(s::AbstractDesignMatrix, state = 1)
    return state > s.count ? nothing : (next!(s), state + 1)
end

"""
    OwenDesignMat{T<:Real, I<:Integer, F<:Integer} <: AbstractDesignMatrix

Create an Owen scrambling iterator for doing multiple i.i.d. [`OwenScrambling`](@ref) randomization.
One can use the common[`DesignMatrix`](@ref) interface to create the iterator.
"""
mutable struct OwenDesignMat{T <: Real, I <: Integer, F <: Integer} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    random_bits::Array{I, 3} #array of size (pad, N, d)
    bits::Array{I, 3} #array of size (pad, N, d)
    indices::Array{F, 3} #array of size (m, N, d)
    R::OwenScramble
    count::Int
end

"""
    ScrambleDesignMat{T<:Real, I<:Integer} <: AbstractDesignMatrix

Create a scrambling iterator (Digital Shift or Matousek depending on the `R` field) for doing multiple i.i.d. [`DigitalShift`](@ref) or [`MatousekScrambling`](@ref) randomization.
One can use the common[`DesignMatrix`](@ref) interface to create the iterator.
"""
mutable struct ScrambleDesignMat{T <: Real, I <: Integer} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    random_bits::Array{I, 3} #array of size (pad, N, d)
    bits::Array{I, 3} #array of size (pad, N, d)
    R::ScrambleMethod
    count::Int
end

"""
    ShiftDesignMat{T<:Real} <: AbstractDesignMatrix

Create a Shift iterator for doing multiple i.i.d [`Shift`](@ref) randomization.
One can use the common[`DesignMatrix`](@ref) interface to create the iterator.
"""
mutable struct ShiftDesignMat{T <: Real} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    R::Shift
    count::Int
end

"""
    DistributionDesignMat{T<:Real} <: AbstractDesignMatrix

Create an iterator for multiple distribution randomization. The distribution is chosen with the field `D`.
This is equivalent to using `rand!(D, X)` for some matrix `X`.
One can use the common[`DesignMatrix`](@ref) interface to create the iterator.
"""
mutable struct DistributionDesignMat{T <: Real, T2} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    D::T2 #::Distributions.Sampleable
    count::Int
end

"""
    RandomDesignMat{T<:Real} <: AbstractDesignMatrix

Create an iterator for multiple uniform randomization. This it similar to [`DistributionDesignMat`](@ref) with the field `D = Uniform()`
One can use the common[`DesignMatrix`](@ref) interface to create the iterator.
"""
mutable struct RandomDesignMat{T <: Real} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    count::Int
end

Base.eltype(::Type{OwenDesignMat{T, I, F}}) where {T, I, F} = Matrix{T}
Base.eltype(::Type{ScrambleDesignMat{T, I}}) where {T, I} = Matrix{T}
Base.eltype(::Type{ShiftDesignMat{T}}) where {T} = Matrix{T}
Base.eltype(::Type{DistributionDesignMat{T}}) where {T} = Matrix{T}
Base.eltype(::Type{RandomDesignMat{T}}) where {T} = Matrix{T} # TODO one could create a type for AbstractDistribution to include RandomDesignMat and DistributionDesignMat

"""
    DesignMatrix(n, d, sampler, num_mats, T = Float64)
    DesignMatrix(n, d, sampler, randomization, num_mats, T = Float64)

Create an iterator that produces `num_mats` point-set matrices.

The iterator is useful when a computation needs independent randomized QMC
realizations but allocating all matrices at once is undesirable. For a
deterministic sampler, the five-argument form uses the sampler's `R` field;
the six-argument form selects `randomization` explicitly. `RandomSample` also
has a five-argument form and generates independent Monte Carlo matrices.

# Arguments

- `n::Integer`: Number of points in each matrix.
- `d::Integer`: Dimension of each point. Each yielded matrix has size `(d, n)`.
- `sampler`: Deterministic or random [`SamplingAlgorithm`](@ref).
- `randomization::RandomizationMethod`: Randomization applied to a
  deterministic sampler.
- `num_mats::Integer`: Number of matrices yielded by the iterator.
- `T::DataType = Float64`: Element type of each matrix.

# Returns

- An [`AbstractDesignMatrix`](@ref) iterator whose elements are matrices of
  size `(d, n)` and element type `T`.

# Throws

- An exception from the selected sampler or randomization method if `n`, `d`,
  the randomization base, or the sampler-specific parameters are invalid.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> iterator = DesignMatrix(4, 2, SobolSample(R = Shift()), 3);

julia> length(iterator)
3

julia> size(first(iterator))
(2, 4)
```

# Extension rules

For a custom deterministic sampler, define the unit-box `sample` method from
[`SamplingAlgorithm`](@ref) and provide an `R::RandomizationMethod` field.
The built-in randomization-specific `initialize` methods are implementation
details; custom extensions should use the public `sample` and `randomize`
contracts unless they are deliberately implementing a compatible iterator.
"""
function DesignMatrix(
        N, d, S::DeterministicSamplingAlgorithm, num_mats::Integer, T::DataType = Float64
    )
    return DesignMatrix(N, d, S, S.R, num_mats, T)
end

function DesignMatrix(
        N,
        d,
        S::DeterministicSamplingAlgorithm,
        R::OwenScramble,
        num_mats::Integer,
        T::DataType = Float64
    )
    X, random_bits, bits, indices = initialize(N, d, S, R, T)
    return OwenDesignMat(X, random_bits, bits, indices, R, num_mats)
end

function DesignMatrix(
        N,
        d,
        S::DeterministicSamplingAlgorithm,
        R::ScrambleMethod,
        num_mats::Integer,
        T = Float64
    )
    X, random_bits, bits = initialize(N, d, S, R, T)
    return ScrambleDesignMat(X, random_bits, bits, R, num_mats)
end

function DesignMatrix(
        N,
        d,
        S::DeterministicSamplingAlgorithm,
        R::Shift,
        num_mats,
        T::DataType = Float64
    )
    X = initialize(N, d, S, R, T)
    return ShiftDesignMat(X, R, num_mats)
end

function DesignMatrix(N, d, D::RandomSample, num_mats, T::DataType = Float64)
    X = initialize(N, d, D, T)
    return RandomDesignMat(X, num_mats)
end

next!(DM::OwenDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.indices, DM.R)
next!(DM::ScrambleDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.R)
next!(DM::ShiftDesignMat) = randomize!(DM.X, DM.R)
next!(DM::DistributionDesignMat) = rand!(DM.D, DM.X)
next!(DM::RandomDesignMat) = rand!(DM.X)

## OwenScramble
function initialize(n, d, sampler, R::OwenScramble, T::DataType = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = permutedims(sample(n, d, no_rand_sampler, T))

    b = R.base
    bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(bits)
    random_points = similar(points)
    indices = which_permutation(bits, R.base)
    return random_points, random_bits, bits, indices
end

function scramble!(
        random_points::AbstractArray{T},
        random_bits,
        bits,
        indices,
        R::OwenScramble
    ) where {T <: Real}
    randomize_bits!(random_bits, bits, indices, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end
    return permutedims(random_points)
end

## Other scramble

function initialize(n, d, sampler, R::ScrambleMethod, T::DataType = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = permutedims(sample(n, d, no_rand_sampler, T))

    b = R.base
    bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(bits)
    random_points = similar(points)
    return random_points, random_bits, bits
end

function scramble!(
        random_points::AbstractArray{T},
        random_bits,
        bits,
        R::ScrambleMethod
    ) where {T <: Real}
    randomize_bits!(random_bits, bits, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end
    return permutedims(random_points)
end

## Shift
function initialize(n, d, sampler, R::Shift, T::DataType = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = sample(n, d, no_rand_sampler, T)
    return points
end

## Distribution
function initialize(n, d, D::RandomSample, T::DataType = Float64)
    # Generate unrandomized sequence
    X = zeros(T, d, n)
    return X
end

# generate_design_matrices

## Generic function

"""
    generate_design_matrices(n, d, sampler, num_mats, T = Float64)
    generate_design_matrices(n, lb, ub, sampler, num_mats = 2)

Generate multiple point-set matrices with `sampler`.

# Arguments

  - `n`: Positive number of points in each generated matrix.
  - `d`: Positive dimension of the unit box. This form returns matrices of size
    `(d, n)` with elements in `[0, 1]`.
  - `lb`: Collection of lower bounds. Its length determines the dimension.
  - `ub`: Collection of upper bounds with the same length as `lb`.
  - `sampler`: Concrete [`SamplingAlgorithm`](@ref) used to construct each
    point set.
  - `num_mats`: Number of matrices to generate.

# Optional positional arguments

  - `T = Float64`: Element type of the unit-box matrices.

# Returns

A vector of `num_mats` matrices. Each unit-box matrix has size `(d, n)`. The
bounds form maps every coordinate to its corresponding interval `[lb[i], ub[i]]`.

# Throws

- `AssertionError`: If the bounds have different lengths or a lower bound
  exceeds its corresponding upper bound.
- An exception from the selected sampler if `n`, `d`, or `num_mats` is invalid.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> matrices = generate_design_matrices(4, 2, RandomSample(), 3);

julia> length(matrices)
3

julia> all(size(matrix) == (2, 4) for matrix in matrices)
true
```

# Developer Notes

This function builds on the public [`sample`](@ref) contract. New sampling
algorithms should extend `sample`; they should not add methods to
`generate_design_matrices`.
"""
function generate_design_matrices(
        n, d, sampler::DeterministicSamplingAlgorithm, num_mats::Integer,
        T::DataType = Float64
    )
    return generate_design_matrices(n, d, sampler, sampler.R, num_mats, T)
end

function generate_design_matrices(
        n, d, sampler::RandomSamplingAlgorithm, num_mats::Integer,
        T::DataType = Float64
    )
    return [sample(n, d, sampler, T) for j in 1:num_mats]
end

function generate_design_matrices(
        n, lb, ub, sampler,
        num_mats = 2
    )
    if n <= 0
        throw(ZeroSamplesError())
    end
    @assert length(lb) == length(ub)

    # Generate a vector of num_mats independent "randomized" version of the QMC sequence
    out = generate_design_matrices(n, length(lb), sampler, num_mats, eltype(lb))

    # Scaling
    for j in eachindex(out)
        out[j] = (ub .- lb) .* out[j] .+ lb
    end
    return out
end

## NoRand()

"""
    generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)

`R = NoRand()` produces `num_mats` matrices each containing a different deterministic point set in `[0, 1)ᵈ`.
Note that this is an ad hoc way to produce i.i.d sequence as it creates a deterministic point in dimension `d × num_mats` and split it in `num_mats` point set of dimension `d`.
This does not have any QMC garantuees.
"""
function generate_design_matrices(
        n, d, sampler, R::NoRand, num_mats::Integer, T::DataType = Float64
    )
    out = sample(n, num_mats * d, sampler, T)
    @warn "The `generate_design_matrices(n, d, sampler, R = NoRand(), num_mats)` method does not produces true and independent QMC matrices, see [this doc warning](https://docs.sciml.ai/QuasiMonteCarlo/stable/design_matrix/) for more context.
    Prefer using randomization methods such as `R = Shift()`, `R = MatousekScrambling()`, etc., see [documentation](https://docs.sciml.ai/QuasiMonteCarlo/stable/randomization/)"
    return [out[(j * d + 1):((j + 1) * d), :] for j in 0:(num_mats - 1)]
end

function generate_design_matrices(
        n,
        d,
        sampler,
        R::RandomizationMethod,
        num_mats::Integer,
        T::DataType = Float64
    )
    return collect(DesignMatrix(n, d, sampler, R, num_mats, T))
end
