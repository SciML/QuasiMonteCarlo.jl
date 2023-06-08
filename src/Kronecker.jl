"""
    KroneckerSample(generator::AbstractVector, R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

A Kronecker sequence is a point set generated using a vector and equation:
`x[i] = i * generator .% 1`

Where `i` runs from `1` through the sample size `n`. This sequence will be equidistributed
(uniform in the infinite limit) so long as the components of `generator` are linearly
independent over the field of rational numbers.

If no generator is specified, a lattice based on the generalized golden ratio is used; see
`GoldenSample` for more information.

Kronecker sequences are not recommended for use in more than 3 dimensions, as theory on them
is sparse. `LatticeRuleSample` will return rank-1 lattice rules, which behave similarly to
Kronecker sequences but have better properties.

References:
Leobacher, G., & Pillichshammer, F. (2014). *Introduction to quasi-Monte Carlo integration and applications.* Switzerland: Springer International Publishing.
https://link.springer.com/content/pdf/10.1007/978-3-319-03425-6.pdf
"""
Base.@kwdef struct KroneckerSample{V <: Union{AbstractVector, Missing}} <:
                   DeterministicSamplingAlgorithm
    generator::V = missing
    R::RandomizationMethod = NoRand()
end

function KroneckerSample(d::Integer, R = NoRand(), T = Float64)
    ratio = harmonious(d, 2eps(T))
    generator = [ratio^i for i in 1:d]
    return KroneckerSample(generator, R)
end

function sample(n::Integer, d::Integer, S::KroneckerSample{Missing}, T = Float64)
    return sample(n, d, KroneckerSample(d, S.R, T))
end

function sample(n::Integer, d::Integer, k::KroneckerSample{V},
                T) where {V <: AbstractVector}
    @assert eltype(V)==T "Sample type must match generator type."
    return randomize(sample(n, d, k), k.R)
end

function sample(n::Integer, d::Integer, k::KroneckerSample{V}) where {V <: AbstractVector}
    @assert d==length(k.generator) "Dimensions of generator and sample must match."
    return @. mod(k.generator * (1:n)', 1)
end

"""
    GoldenSample()

Generate a quasirandom Kronecker sequence using powers of the generalized golden ratio.

The harmonious, or generalized golden, ratios are defined as the solutions to the equation:
``x^d = x + 1``

Where `d` is the dimension of the sequence. The Golden sequence is then equivalent to
`Kronecker([x^-i for i in 1:d])`.

WARNING: the generalized golden sequence in more than 2 dimensions is purely experimental.
It likely has poor discrepancy in high dimensions, and should not be used without verifying
answers against a better-known quasirandom sequence. Try a rank-1 lattice rule instead.

References:
Roberts, M. (2018). The Unreasonable Effectiveness of Quasirandom Sequences.
*Extreme Learning.*
http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
"""
GoldenSample(args...; kwargs...) = KroneckerSample(args...; kwargs...)

# generate (inverse) harmonious numbers
@fastmath function harmonious(d::Integer, tol::T = 2eps(float(d))) where {T}
    y_old = one(T)
    # Approximate solution of x^(d+1) = x + 1 (nested radical)
    y = (one(T) + y_old)^inv(d + 1)
    while abs(y - y_old) > tol  # continue until precise enough
        y_old = y
        y = (one(T) + y)^inv(d + 1)
    end
    return inv(y)
end
