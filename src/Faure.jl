# upper-triangular Pascal matrix--needed to generate Faure points
function pascal_mat(dimension::I, base = Inf) where {I <: Integer}
    pascal = UpperTriangular(ones(I, dimension, dimension))
    @inbounds for i in 2:dimension
        for j in 2:i
            pascal[j, i] = (pascal[j - 1, i - 1] + pascal[j, i - 1]) % base
        end
    end
    return pascal
end

# raise a Pascal matrix to a given power.
@fastmath function pascal_power!(result::UpperTriangular, pascal::UpperTriangular,
                                 power::Integer, base::Integer)
    @inbounds @simd for idx in eachindex(pascal)
        i, j = Tuple(idx)
        i ≤ j ? (result[idx] = powermod(power, j - i, base) * pascal[idx]) : result[idx]
    end
    return result
end

"""
    FaureSample(R::RandomizationMethod)

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
Faure, H. (1982). Discrepance de suites associées a un systeme de numération (en dimension s). *Acta Arith.*, 41, 337-351.
Owen, A. B. (1997). Monte Carlo variance of scrambled net quadrature. *SIAM Journal on Numerical Analysis*, 34(5), 1884-1910.
"""
Base.@kwdef struct FaureSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

@fastmath function sample(n::Integer, dimension::Integer, S::FaureSample, T = Float64;
                          skipchecks = false)
    base = nextprime(dimension)
    n_digits = ceil(Int, log(base, n))
    power = n_digits - 1

    if !skipchecks && (n % prevpow(base, n) ≠ 0)
        n -= (n % base^power)
        throw(ArgumentError("Invalid sample size: Faure sequences must be multiples of `base^power`. " *
                            "Try $n or $(n+base^power) instead."))
    end

    return randomize(_faure_samples(n, n_digits, dimension, T), S.R)
end

@fastmath @views function _faure_samples(n_samples::I, n_digits::I, dimension::I,
                                         ::Type{F} = Float64) where {I <: Integer, F}
    base = nextprime(dimension)
    inv_base = inv(base)

    # Upper triangular Pascal matrix
    pascal = pascal_mat(n_digits, base)
    # we multiply by powers of a Pascal matrix to generate it. preallocating:
    permutation = similar(pascal)
    # preallocate space for Faure points and a vector of digits
    faure = zeros(F, dimension, n_samples)
    digit_counter = zeros(I, n_digits)
    dgs = copy(digit_counter)
    for dim_idx in 1:dimension
        digit_counter .= base - 1
        dim_idx == 1 || pascal_power!(permutation, pascal, dim_idx - 1, base)
        # each dimension of a Faure sequence is a permuted van der Corput (vdc) sequence
        vdc = faure[dim_idx, :]
        for sample_idx in eachindex(vdc)
            # increment digit counter by 1
            _base_sum!(digit_counter, base)
            dgs .= digit_counter
            # permute later dimensions using powers of a Pascal matrix
            if dim_idx ≠ 1
                dgs .= (permutation * dgs) .% base
            end
            vdc[sample_idx] = (evalpoly(inv_base, dgs) + inv(2base^n_digits)) * inv_base
        end
    end
    return faure
end

# @fastmath @views function _base_sum!(dgs::AbstractVector, base::Integer)
#     dgs[1] = one(eltype(dgs)) + dgs[1]
#     if dgs[1] == base
#         dgs[1] = zero(eltype(dgs))
#         if length(dgs) ≠ 1  # ignore overflows
#             _base_sum!(dgs[2:end], base)
#         end
#     end
# end
