import Base.size, Base.getindex, Base.setindex, Base.axes

# Pascal triangle -- implementation needed for Faure sequences
struct PascalMatrix{I, M<:AbstractMatrix{I}} <: LinearAlgebra.AbstractTriangular{I, M}
    data::M
end

function PascalMatrix(dimension::I) where I<:Integer
    pascal = UpperTriangular(ones(I, dimension, dimension))
	for i in 2:dimension
        for j in 2:i
            pascal[j, i] = pascal[j-1, i-1] + pascal[j, i-1]
        end
    end
    return PascalMatrix(pascal)
end

Base.:size(x::PascalMatrix) = size(x.data)
Base.:getindex(x::PascalMatrix, inds...) = getindex(x.data, inds...)
Base.:setindex!(A, X::PascalMatrix, inds...) = setindex!(A, X.data, inds...)
Base.:axes(x::PascalMatrix) = Base.:axes(x.data)

"""
    FaureSample()

A Faure low-discrepancy sequence.

Faure-distributed samples cover all dimensions evenly, using the same set of points for all
variables, up to ordering.

When scrambled, randomized Faure sequences provide worst-case guarantees that variance will
be at most `exp(1) ≈ 2.718` times greater than for a purely Monte Carlo integral. However,
they are much less efficient at integrating sequences with low effective dimension.

The Faure sequence in dimension `s` forms a `(0, s)`-sequence with base `b = nextprime(s)`.

A Faure sequence must have length `k * base^s` with `k < base < 1`.

References:
Faure, H. (1982). Discrepance de suites associees a un systeme de numeration (en dimension s). *Acta Arith.*, 41, 337-351.
Owen, A. B. (1997). Monte Carlo variance of scrambled net quadrature. *SIAM Journal on Numerical Analysis*, 34(5), 1884-1910.
"""
struct FaureSample end

function Base.:^(pascal::PascalMatrix, power::Integer)
    result = similar(pascal.data)
    @inbounds for idx in eachindex(pascal.data)
        i, j = Tuple(idx)
        result[i, j] = power^(max(j-i, 0)) * pascal.data[i, j]
    end
    return result
end

function sample(n::Integer, lb, ub, ::FaureSample)
    length(lb) == length(ub) || DimensionMismatch("Lower and upper bounds do not match.")
    dimension = length(lb)
    faure = sample(n, dimension, FaureSample())
    @inbounds for (row_idx, row) in enumerate(eachrow(faure))
        @. row = (ub[row_idx] - lb[row_idx]) * row - lb
    end
end

function sample(n::Integer, dimension::Integer, ::FaureSample; skipchecks=false)
    base = nextprime(dimension)
    power = floor(Int, log(base, n))

    if !skipchecks && (n % prevpow(base, n) ≠ 0)
        n -= (n % base^power)
        throw(ArgumentError(
            "Invalid sample size: Faure sequences must be multiples of `base^power`. " *
            "Try $n or $(n+base^power) instead."
        ))
    end

    return _faure_samples(n, power, dimension)
end

@views @fastmath function _faure_samples(
    n_samples::Integer, power::Integer, dimension::Integer, ::Type{F}=Float64) where F
    base = nextprime(dimension)
    n_digits = power + 1

    # Upper triangular Pascal matrix
    pascal = PascalMatrix(n_digits)

    faure = Matrix{F}(undef, dimension, n_samples)
    @inbounds for sample_idx in 1:n_samples
        # base decomposition
        dgs = digits(sample_idx; base=base, pad=n_digits)
        faure[1, sample_idx] = (evalpoly(inv(base), dgs) / base) % base
        @simd for dim_idx in 2:dimension
            digit_vec = pascal^(dim_idx - 1) * dgs .% base
            faure[dim_idx, sample_idx] = (evalpoly(inv(base), digit_vec) / base) % base
        end
    end
    return faure
end
