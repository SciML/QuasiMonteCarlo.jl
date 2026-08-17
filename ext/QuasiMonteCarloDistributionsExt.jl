module QuasiMonteCarloDistributionsExt

using QuasiMonteCarlo
import Distributions

"""
    sample(n::Integer, d::Integer, D::Distributions.Sampleable, T = eltype(D))

Draw independent samples from a `Distributions.Sampleable` distribution.

# Arguments

- `n::Integer`: Number of points to draw.
- `d::Integer`: Number of independent coordinates in each point.
- `D::Distributions.Sampleable`: Distribution sampled independently for each
  coordinate.
- `T = eltype(D)`: Element-type argument retained for compatibility with the
  general `sample` interface. The distribution controls the element type of
  the generated values.

# Returns

- A matrix with size `(d, n)`, with one sampled point per column.

# Examples

```jldoctest
julia> using Distributions, QuasiMonteCarlo

julia> points = sample(4, 2, Normal());

julia> size(points)
(2, 4)
```
"""
function QuasiMonteCarlo.sample(
        n::Integer,
        d::Integer,
        D::Distributions.Sampleable,
        T = eltype(D)
    )
    @assert n > 0 QuasiMonteCarlo.ZERO_SAMPLES_MESSAGE
    x = [[rand(D) for j in 1:d] for i in 1:n]
    return reduce(hcat, x)
end

"""
    sample(n::Integer, lb, ub, sampler::Distributions.Sampleable)

Draw a quasi-Monte Carlo point set and map it to the box defined by `lb` and
`ub` using a `Distributions.Sampleable` as the sampler argument.

# Arguments

- `n::Integer`: Number of points to draw.
- `lb`: Scalar, tuple, or vector of lower bounds.
- `ub`: Upper bounds with the same shape as `lb`.
- `sampler::Distributions.Sampleable`: Distribution used by the extension.

# Returns

- A matrix with size `(length(lb), n)` whose columns lie in the requested
  bounds.

# Throws

- `AssertionError`: If the bounds have different lengths or a lower bound
  exceeds its corresponding upper bound.
"""
function QuasiMonteCarlo.sample(
        n::Integer, lb::T, ub::T,
        S::D
    ) where {
        T <: Union{AbstractVector, Tuple, Number},
        D <: Distributions.Sampleable,
    }
    QuasiMonteCarlo._check_sequence(lb, ub, n)
    lb = float.(lb)
    ub = float.(ub)
    out = QuasiMonteCarlo.sample(n, length(lb), S, eltype(lb))
    return (ub .- lb) .* out .+ lb
end

function QuasiMonteCarlo.DesignMatrix(
        N,
        d,
        D::Distributions.Sampleable,
        num_mats,
        T = Float64
    )
    X = QuasiMonteCarlo.initialize(N, d, D, T)
    return QuasiMonteCarlo.DistributionDesignMat(X, D, num_mats)
end

function QuasiMonteCarlo.initialize(n, d, D::Distributions.Sampleable, T = Float64)
    # Generate unrandomized sequence
    X = zeros(T, d, n)
    return X
end

end
