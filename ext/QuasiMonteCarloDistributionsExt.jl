module QuasiMonteCarloDistributionsExt

using QuasiMonteCarlo
isdefined(Base, :get_extension) ? (import Distributions) : (import ..Distributions)

"""
```julia
sample(n::Integer, lb::T, ub::T, D::Distributions.Sampleable, T = eltype(D))
sample(n::Integer,
    lb::T,
    ub::T,
    D::Distributions.Sampleable) where {T <: Union{Base.AbstractVecOrTuple, Number}}
```

Return a point set from a distribution `D`:

  - `n` is the number of points to sample.
  - `D` is a `Distributions.Sampleable` from Distributions.jl.
    The point set is in a `d`-dimensional unit box `[0, 1]^d`.
    If the bounds are specified instead of just `d`, the sample is transformed (translation + scaling) into a box `[lb, ub]` where:
  - `lb` is the lower bound for each variable. Its length fixes the dimensionality of the sample.
  - `ub` is the upper bound. Its dimension must match `length(lb)`.
"""
function QuasiMonteCarlo.sample(n::Integer,
        d::Integer,
        D::Distributions.Sampleable,
        T = eltype(D))
    @assert n>0 QuasiMonteCarlo.ZERO_SAMPLES_MESSAGE
    x = [[rand(D) for j in 1:d] for i in 1:n]
    return reduce(hcat, x)
end

"""
```julia
sample(n::Integer, d::Integer, S::Distributions.Sampleable, T = Float64)
sample(n::Integer,
    lb::T,
    ub::T,
    S::Distributions.Sampleable) where {T <: Union{Base.AbstractVecOrTuple, Number}}
```

Return a QMC point set where:

  - `n` is the number of points to sample.
  - `S` is the quasi-Monte Carlo sampling strategy.
    The point set is in a `d`-dimensional unit box `[0, 1]^d`.
    If the bounds are specified, the sample is transformed (translation + scaling) into a box `[lb, ub]` where:
  - `lb` is the lower bound for each variable. Its length fixes the dimensionality of the sample.
  - `ub` is the upper bound. Its dimension must match `length(lb)`.

In the first method the type of the point set is specified by `T` while in the second method the output type is inferred from the bound types.
"""
function QuasiMonteCarlo.sample(n::Integer, lb::T, ub::T,
        S::D) where {T <: Union{Base.AbstractVecOrTuple, Number},
        D <: Distributions.Sampleable}
    QuasiMonteCarlo._check_sequence(lb, ub, n)
    lb = float.(lb)
    ub = float.(ub)
    out = QuasiMonteCarlo.sample(n, length(lb), S, eltype(lb))
    return (ub .- lb) .* out .+ lb
end

function QuasiMonteCarlo.DesignMatrix(N,
        d,
        D::Distributions.Sampleable,
        num_mats,
        T = Float64)
    X = QuasiMonteCarlo.initialize(N, d, D, T)
    return QuasiMonteCarlo.DistributionDesignMat(X, D, num_mats)
end

function QuasiMonteCarlo.initialize(n, d, D::Distributions.Sampleable, T = Float64)
    # Generate unrandomized sequence
    X = zeros(T, d, n)
    return X
end

end
