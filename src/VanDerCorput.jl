"""
    VanDerCorputSample(base::Integer)

The van der Corput sequence, also called the radical inverse sequence, is a one-dimensional
low-discrepancy sequence that recursively splits the unit interval into equally-sized
pieces before inserting one sample in each interval.

For example, in base 2, the sequence starts by inserting one point in each half of the unit
interval in the first pass (the first 2 samples); then one sample in each quarter, then each
eighth, and so on. This creates a well-stratified sample, so long as the number of samples
is a multiple of a power of the base.
"""
struct VanDerCorputSample{I<:Integer} <: SamplingAlgorithm
    base::I
end

function sample(n::Integer, d::Integer, S::VanDerCorputSample, T::Type=Float64)
    @assert d == 1 "Van der Corput sequence only supports 1D sampling"
    @assert n > 0 ZERO_SAMPLES_MESSAGE
    base = S.base
    log_n = log(base, n)
    t = floor(Int, log_n)
    n_digits = ceil(Int, log_n)
    λ = n ÷ (base^t)
    return _vdc(promote(λ, n_digits, base)..., T; n)
end

"""
    _vdc(λ::Integer, n_digits::Integer, base::Integer, F=Float64)
    _vdc(n::Integer, base::Integer, F=Float64)

Generate a Van der Corput sequence of length `n=λ*base^n_digits` in the given base, where
`λ` is the number of samples per subinterval, and `base^n_digits` is the number of subintervals.

Note that `λ` must be smaller than `base` for the sequence to be well-stratified.
"""
function _vdc(n::I, base::I, F=Float64) where {I <: Integer}
    return _vdc(one(n), I(log(base, n)), base, F; n)
end
@fastmath @views function _vdc(
    λ::I, n_digits::I, base::I, F=Float64; n=λ * base^n_digits
) where {I <: Integer}
    sequence = zeros(F, n)
    inv_base = convert(F, inv(base))
    dgs = fill(base - 1, n_digits)
    # offset required to place samples in center of interval
    offset = convert(F, inv(2n))
    for sample_idx in eachindex(sequence)
        # increment digit counter by 1
        _base_sum!(dgs, base)
        sequence[sample_idx] = evalpoly(inv_base, dgs) * inv_base + offset
    end
    return sequence
end

@fastmath @views function _base_sum!(dgs::AbstractVector, base::Integer)
    @label start # for recursion
    dgs[1] = one(eltype(dgs)) + dgs[1]
    if dgs[1] == base
        dgs[1] = zero(eltype(dgs))
        if length(dgs) ≠ 1  # ignore overflows
            dgs = dgs[2:end]
            @goto start  # recursive call
        end
    end
end
