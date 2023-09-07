abstract type AbstractDesignMatrix end

# Make an iterator so that we can do "for X in DesignMatrix(...)"
Base.length(s::AbstractDesignMatrix) = s.count
function Base.iterate(s::AbstractDesignMatrix, state = 1)
    state > s.count ? nothing : (next!(s), state + 1)
end

"""
    OwenDesignMat{T<:Real, I<:Integer, F<:Integer} <: AbstractDesignMatrix
Owen scrambling iterator. 
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
Scrambling iterator used in Digital Shift and Matousek scrambling. 
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
Scrambling iterator used in shift method. 
"""
mutable struct ShiftDesignMat{T <: Real} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    R::Shift
    count::Int
end

mutable struct DistributionDesignMat{T <: Real} <: AbstractDesignMatrix
    X::Array{T, 2} #array of size (N, d)
    D::Distributions.Sampleable
    count::Int
end

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
```julia
DesignMatrix(n, d, sample_method::DeterministicSamplingAlgorithm, num_mats, T = Float64)
```
Create an iterator for doing mutliple randomization over QMC sequences where
- `num_mats` is the lenght of the iterator
- `n` is the number of points to sample.
- `d` is the dimensionality of the point set in `[0, 1)ᵈ`,
- `sample_method` is the quasi-Monte Carlo sampling strategy used to create a deterministic point set `out`.
- `T` is the `eltype` of the point sets. For some QMC methods (Faure, Sobol) this can be `Rational`
It is now compatible with all scrambling methods and shifting. One can also use it with `Distributions.Sampleable` or `RandomSample`.
"""
function DesignMatrix(N, d, S::DeterministicSamplingAlgorithm, num_mats, T = Float64)
    return DesignMatrix(N, d, S, S.R, num_mats, T)
end

function DesignMatrix(N,
    d,
    S::DeterministicSamplingAlgorithm,
    R::OwenScramble,
    num_mats,
    T = Float64)
    X, random_bits, bits, indices = initialize(N, d, S, R, T)
    return OwenDesignMat(X, random_bits, bits, indices, R, num_mats)
end

function DesignMatrix(N,
    d,
    S::DeterministicSamplingAlgorithm,
    R::ScrambleMethod,
    num_mats,
    T = Float64)
    X, random_bits, bits = initialize(N, d, S, R, T)
    return ScrambleDesignMat(X, random_bits, bits, R, num_mats)
end

function DesignMatrix(N,
    d,
    S::DeterministicSamplingAlgorithm,
    R::Shift,
    num_mats,
    T = Float64)
    X = initialize(N, d, S, R, T)
    return ShiftDesignMat(X, R, num_mats)
end

function DesignMatrix(N, d, D::Distributions.Sampleable, num_mats, T = Float64)
    X = initialize(N, d, D, T)
    return DistributionDesignMat(X, D, num_mats)
end

function DesignMatrix(N, d, D::RandomSample, num_mats, T = Float64)
    X = initialize(N, d, D, T)
    return RandomDesignMat(X, num_mats)
end

next!(DM::OwenDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.indices, DM.R)
next!(DM::ScrambleDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.R)
next!(DM::ShiftDesignMat) = randomize!(DM.X, DM.R)
next!(DM::DistributionDesignMat) = rand!(DM.D, DM.X)
next!(DM::RandomDesignMat) = rand!(DM.X)

## OwenScramble
function initialize(n, d, sampler, R::OwenScramble, T = Float64)
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

function scramble!(random_points::AbstractArray{T},
    random_bits,
    bits,
    indices,
    R::OwenScramble) where {T <: Real}
    randomize_bits!(random_bits, bits, indices, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end
    return permutedims(random_points)
end

## Other scramble

function initialize(n, d, sampler, R::ScrambleMethod, T = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = permutedims(sample(n, d, no_rand_sampler, T))

    b = R.base
    bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(bits)
    random_points = similar(points)
    return random_points, random_bits, bits
end

function scramble!(random_points::AbstractArray{T},
    random_bits,
    bits,
    R::ScrambleMethod) where {T <: Real}
    randomize_bits!(random_bits, bits, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end
    return permutedims(random_points)
end

## Shift
function initialize(n, d, sampler, R::Shift, T = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = sample(n, d, no_rand_sampler, T)
    return points
end

## Distribution
function initialize(n, d, D::Union{Distributions.Sampleable, RandomSample}, T = Float64)
    # Generate unrandomized sequence
    X = zeros(T, d, n)
    return X
end

# generate_design_matrices

## Generic function 

"""
```julia
generate_design_matrices(n, d, sample_method::DeterministicSamplingAlgorithm,
num_mats, T = Float64)
generate_design_matrices(n, d, sample_method::RandomSamplingAlgorithm,
num_mats, T = Float64)
generate_design_matrices(n, lb, ub, sample_method,
num_mats = 2)
```
Create `num_mats` matrices each containing a QMC point set, where:
- `n` is the number of points to sample.
- `d` is the dimensionality of the point set in `[0, 1)ᵈ`,
- `sample_method` is the quasi-Monte Carlo sampling strategy used to create a deterministic point set `out`.
- `T` is the `eltype` of the point sets. For some QMC methods (Faure, Sobol) this can be `Rational`
If the bound `lb` and `ub` are specified instead of `d`, the samples will be transformed into the box `[lb, ub]`.
"""
function generate_design_matrices(n, d, sampler::DeterministicSamplingAlgorithm, num_mats,
    T = Float64)
    return generate_design_matrices(n, d, sampler, sampler.R, num_mats, T)
end

function generate_design_matrices(n, d, sampler::RandomSamplingAlgorithm, num_mats,
    T = Float64)
    return [sample(n, d, sampler, T) for j in 1:num_mats]
end

function generate_design_matrices(n, lb, ub, sampler,
    num_mats = 2)
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
function generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)
    out = sample(n, num_mats * d, sampler, T)
    return [out[(j * d + 1):((j + 1) * d), :] for j in 0:(num_mats - 1)]
end

function generate_design_matrices(n,
    d,
    sampler,
    R::RandomizationMethod,
    num_mats,
    T = Float64)
    return collect(DesignMatrix(n, d, sampler, R, num_mats, T))
end