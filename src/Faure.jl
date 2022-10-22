# lower-triangular Pascal matrix--needed to generate Faure points
function pascal_mat(dimension::I, base=Inf) where I<:Integer
    pascal = UpperTriangular(ones(I, dimension, dimension))
	@inbounds for i in 2:dimension
        for j in 2:i
            pascal[j, i] = (pascal[j-1, i-1] + pascal[j, i-1]) % base
        end
    end
    return pascal
end

# raise a Pascal matrix to a given power.
@fastmath function pascal_power!(result::UpperTriangular, pascal::UpperTriangular, power::Integer, base::Integer)
    @inbounds @simd for idx in eachindex(pascal)
        i, j = Tuple(idx)
        result[idx] = powermod(power, max(j-i, zero(i)), base) * pascal[idx]
    end
    return result
end

"""
    FaureSample()

A Faure low-discrepancy sequence.

Faure-distributed samples cover all dimensions evenly, using the same set of points for all
variables, up to ordering.

When scrambled, randomized Faure sequences provide worst-case guarantees that variance will
be at most `exp(1) ≈ 2.718` times greater than for a purely Monte Carlo integral. However,
they are much less efficient than the Sobol sequence at integrating functions with low
effective dimension (functions where the first few inputs dominate the evaluation).

The Faure sequence in dimension `s` forms a `(0, s)`-sequence with base `b = nextprime(s)`.

A Faure sequence must have length `k * base^s` with `k < base < 1`.

References:
Faure, H. (1982). Discrepance de suites associees a un systeme de numeration (en dimension s). *Acta Arith.*, 41, 337-351.
Owen, A. B. (1997). Monte Carlo variance of scrambled net quadrature. *SIAM Journal on Numerical Analysis*, 34(5), 1884-1910.
"""
struct FaureSample{T}
    rng::T
end

FaureSample() = FaureSample(Random.GLOBAL_RNG)


function sample(n::Integer, lb, ub, ::FaureSample)
    length(lb) == length(ub) || DimensionMismatch("Lower and upper bounds do not match.")
    dimension = length(lb)
    faure = sample(n, dimension, FaureSample())
    @inbounds for (row_idx, row) in enumerate(eachrow(faure))
        @. row = (ub[row_idx] - lb[row_idx]) * row - lb[row_idx]
    end
    return faure
end

function sample(n::Integer, dimension::Integer, sampler::FaureSample; skipchecks=false)
    base = nextprime(dimension)
    power = floor(Int, log(base, n))

    if !skipchecks && (n % prevpow(base, n) ≠ 0)
        n -= (n % base^power)
        throw(ArgumentError(
            "Invalid sample size: Faure sequences must be multiples of `base^power`. " *
            "Try $n or $(n+base^power) instead."
        ))
    end

    return _faure_samples(sampler.rng, n, power, dimension)
end

@views function _faure_samples(
    rng, n_samples::I, power::I, dimension::I, ::Type{F}=Float64
    ) where {I<:Integer, F}
    base = nextprime(dimension)
    inv_base = inv(base)
    n_digits = I(power + one(power))

    # Upper triangular Pascal matrix
    pascal = pascal_mat(n_digits, base)
    # we multiply by powers of a Pascal matrix to generate it. preallocating:
    permutation = similar(pascal)
    # preallocate space for Faure points and a vector of digits
    faure = zeros(F, dimension, n_samples)
    digit_counter = zeros(I, n_digits)
    dgs = copy(digit_counter)
    # we use index maps to randomize the points of a Faure sequence
    idx_maps = [i for (_, i) in Iterators.product(1:n_digits, 0:(base-1))]

    for dim_idx in 1:dimension
        digit_counter .= base-1
        dim_idx == 1 || pascal_power!(permutation, pascal, dim_idx-1, base)
        # each dimension of a Faure sequence is a permuted van der Corput (vdc) sequence
        vdc = faure[dim_idx, :]
        for sample_idx in eachindex(vdc)
            # increment digit counter by 1
            _base_sum!(rng, digit_counter, base, idx_maps)
            dgs .= digit_counter
            # permute later dimensions using powers of a Pascal matrix
            if dim_idx ≠ 1
                dgs .= (permutation * dgs) .% base
            end
            # randomly shuffle points
            for i in eachindex(dgs)
                dgs[i] = idx_maps[1+end-i, 1+dgs[i]]
            end
            vdc[sample_idx] = @evalpoly(inv_base, dgs..., rand(rng, F)) * inv_base
        end
    end
    return faure
end

@views function _base_sum!(
    rng, dgs::AbstractVector, base::Integer, idx_maps::AbstractArray
)
    dgs[1] = one(eltype(dgs)) + dgs[1]
    if dgs[1] == base
        dgs[1] = zero(eltype(dgs))
        shuffle!(rng, idx_maps[1, :])
        if length(dgs) ≠ 1  # ignore overflows
            _base_sum!(rng, dgs[2:end], base, idx_maps[2:end, :])
        end
    end
end
