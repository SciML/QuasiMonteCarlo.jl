"""
    Shift(; rng = Random.TaskLocalRNG()) <: RandomizationMethod

Cranley-Patterson rotation.

# Fields

- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used to
  draw one shift vector for each call to [`randomize`](@ref).

# Examples

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> points = sample(8, 2, SobolSample());

julia> shifted = randomize(points, Shift(rng = MersenneTwister(42)));

julia> size(shifted) == size(points)
true
```

References: Cranley, R., & Patterson, T. N. (1976). Randomization of number theoretic methods for multiple integration. SIAM Journal on Numerical Analysis, 13(6), 904-914.
"""
Base.@kwdef @concrete struct Shift <: RandomizationMethod
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function randomize(x, R::Shift)
    y = copy(x)
    shift!(R.rng, y)
    return y
end

function randomize!(x, R::Shift)
    return shift!(R.rng, x)
end

function shift!(rng::AbstractRNG, points::AbstractMatrix{T}) where {T <: Real}
    d = size(points, 1)
    U = zeros(T, d)
    return shift!(rng, points, U)
end

function shift!(points::AbstractMatrix{T}, U::AbstractVector{T}) where {T <: Real}
    return shift!(Random.TaskLocalRNG(), points, U)
end

function shift!(
        rng::AbstractRNG, points::AbstractMatrix{T},
        U::AbstractVector{T}
    ) where {T <: Real}
    rand!(rng, U)
    for i in axes(points, 2)
        points[:, i] += U
    end
    return points[:] = frac.(points)
end

frac(y) = y - floor(y)
