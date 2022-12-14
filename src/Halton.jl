"""
    HaltonSample()

Create a Halton sequence.
"""
struct HaltonSample <: SamplingAlgorithm end

@views function sample(n::I, d::I, ::HaltonSample, T::Type=Float64) where I<:Integer
    bases = nextprimes(one(n), d)
    n_digits = ceil.(I, log.(bases, n))
    λ = n .÷ bases.^n_digits
    halton_seq = Matrix{T}(undef, d, n)
    for i in axes(halton_seq, 1)
        halton_seq[i, :] .= _vdc(λ[i], n_digits[i], bases[i], T; n)
    end
    return halton_seq
end

# @fastmath @views function vdc(n::I, base::Vector{I}, F=Float64) where I <: Integer
#     n_digits = ceil(Int, log2(n))
#     halton = zeros(F, n, length(base))
#     inv_base = inv.(base)
#     dgs = zeros(I, n_digits)
#     # offset required to place sample in center of interval
#     offset = inv(2n)
#     for sample_idx in eachindex(halton)
#         # increment digit counter by 1
#         _base_sum!(dgs, base)
#         halton[sample_idx, :] .= (evalpoly.(inv_base, dgs) + offset) * inv_base
#     end
#     return halton
# end
