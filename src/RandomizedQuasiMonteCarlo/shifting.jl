"""
```julia
Shifting
```

Cranley-Patterson rotation aka Shifting
"""
Base.@kwdef @concrete struct Shift <: RandomizationMethod
    rng::AbstractRNG = Random.GLOBAL_RNG
end

"""
    randomize(x, R::Shift)
Cranley Patterson Rotation i.e. `y = (points .+ U) mod 1` where `U ∼ 𝕌([0,1]ᵈ)` and `points` is a `d×n` matrix
"""
function randomize(x, R::Shift)
    y = copy(x)
    shift!(R.rng, y)
    return y
end

function shift!(rng::AbstractRNG, points::AbstractMatrix)
    d = size(points, 1)
    U = zeros(d)
    shift!(rng, points, U)
end

shift!(points::AbstractMatrix, U::AbstractVector) = shift!(Random.default_rng(), points, U)

function shift!(rng::AbstractRNG, points::AbstractMatrix, U::AbstractVector)
    rand!(rng, U)
    for i in axes(points, 2)
        points[:, i] += U
    end
    points[:] = frac.(points)
end

frac(y) = y - floor(y)
