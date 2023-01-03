"""
```julia
Shifting
```

Cranley-Patterson rotation aka Shifting
"""
Base.@kwdef struct Shifting <: RandomizationMethod
    rng::AbstractRNG = Random.GLOBAL
end

"""
    shift(points::AbstractArray) 
Cranley Patterson Rotation i.e. `y = (points .+ U) mod 1` where `U ∼ 𝕌([0,1]ᵈ)` and `points` is a `d×n` matrix
"""
function shift(rng::AbstractRNG, points::AbstractArray)
    y = copy(points)
    shift!(rng, y)
    return y
end

shift(points::AbstractArray) = shift(Random.default_rng(), points)
shift!(points::AbstractArray) = shift!(Random.default_rng(), points)
shift!(points::AbstractMatrix, U::AbstractVector) = shift!(Random.default_rng(), points, U)

function shift!(rng::AbstractRNG, points::AbstractMatrix)
    d = size(points, 1)
    U = zeros(d)
    shift!(rng, points, U)
end

function shift!(rng::AbstractRNG, points::AbstractMatrix, U::AbstractVector)
    rand!(rng, U)
    for i in axes(points, 2)
        points[:, i] += U
    end
    points[:] = frac.(points)
end

frac(y) = y - floor(y)