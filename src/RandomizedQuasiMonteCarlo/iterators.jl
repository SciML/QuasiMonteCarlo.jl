abstract type AbstractDesignMatrix end

# Make an iterator so that we can do "for X in design_mat(...)"
Base.length(s::AbstractDesignMatrix) = s.count
Base.iterate(s::AbstractDesignMatrix, state=1) = state > s.count ? nothing : (next!(s), state+1)



mutable struct OwenDesignMat{T<:Real, I<:Integer, F<:Integer} <: AbstractDesignMatrix
    X::Array{T,2} #array of size (N, d)
    random_bits::Array{I,3} #array of size (pad, N, d)
    bits::Array{I,3} #array of size (pad, N, d)
    indices::Array{F,3} #array of size (m, N, d)
    R::OwenScramble
    count::Int
end

mutable struct ScrambleDesignMat{T<:Real, I<:Integer} <: AbstractDesignMatrix
    X::Array{T,2} #array of size (N, d)
    random_bits::Array{I,3} #array of size (pad, N, d)
    bits::Array{I,3} #array of size (pad, N, d)
    R::ScrambleMethod
    count::Int
end

mutable struct ShiftDesignMat{T<:Real} <: AbstractDesignMatrix
    X::Array{T,2} #array of size (N, d)
    R::Shift
    count::Int
end

mutable struct DistributionDesignMat{T<:Real} <: AbstractDesignMatrix
    X::Array{T,2} #array of size (N, d)
    D::Distributions.Sampleable
    count::Int
end

Base.eltype(::Type{OwenDesignMat{T, I, F}}) where {T, I, F} = Matrix{T}
Base.eltype(::Type{ScrambleDesignMat{T, I}}) where {T, I} = Matrix{T}
Base.eltype(::Type{ShiftDesignMat{T}}) where {T} = Matrix{T}
Base.eltype(::Type{DistributionDesignMat{T}}) where {T} = Matrix{T}

function DesignMatrix(N, d, S::DeterministicSamplingAlgorithm, num_mats, T = Float64)
    return DesignMatrix(N, d, S, S.R, num_mats, T)
end

function DesignMatrix(N, d, S::DeterministicSamplingAlgorithm, R::OwenScramble, num_mats, T = Float64)
    X, random_bits, bits, indices = initialize(N, d, S, R, T)
    return OwenDesignMat(X, random_bits, bits, indices, R, num_mats)
end

function DesignMatrix(N, d, S::DeterministicSamplingAlgorithm, R::ScrambleMethod, num_mats, T = Float64)
    X, random_bits, bits = initialize(N, d, S, R, T)
    return ScrambleDesignMat(X, random_bits, bits, R, num_mats)
end

function DesignMatrix(N, d, S::DeterministicSamplingAlgorithm, R::Shift, num_mats, T = Float64)
    X = initialize(N, d, S, R, T)
    return ShiftDesignMat(X, R, num_mats)
end

function DesignMatrix(N, d, D::Distributions.Sampleable, num_mats, T = Float64)
    X = initialize(N, d, D, T)
    return DistributionDesignMat(X, D, num_mats)
end

next!(DM::OwenDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.indices, DM.R)
next!(DM::ScrambleDesignMat) = scramble!(DM.X, DM.random_bits, DM.bits, DM.R)
next!(DM::ShiftDesignMat) = randomize!(DM.X, DM.R)
next!(DM::DistributionDesignMat) = rand!(DM.D, DM.X)

# Owen 
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

function scramble!(random_points::AbstractArray{T}, random_bits, bits, indices, R::OwenScramble) where T<:Real
    randomize_bits!(random_bits, bits, indices, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end

    return permutedims(random_points)
end

# Other scramble
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

function scramble!(random_points::AbstractArray{T}, random_bits, bits, R::ScrambleMethod) where T<:Real
    randomize_bits!(random_bits, bits, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), R.base)
    end
    return permutedims(random_points)
end

# Shift
function initialize(n, d, sampler, R::Shift, T = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    points = sample(n, d, no_rand_sampler, T)
    return points
end

# Distribution
function initialize(n, d, D::Distributions.Sampleable, T = Float64)
    # Generate unrandomized sequence
    X = zeros(T, d, n)
    return X
end
