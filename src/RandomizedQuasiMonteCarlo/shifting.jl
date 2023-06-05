"""
```julia
Shifting(rng::AbstractRNG = Random.GLOBAL_RNG) <: RandomizationMethod
```

Cranley-Patterson rotation aka Shifting

References: Cranley, R., & Patterson, T. N. (1976). Randomization of number theoretic methods for multiple integration. SIAM Journal on Numerical Analysis, 13(6), 904-914.
"""
Base.@kwdef @concrete struct Shift <: RandomizationMethod
    rng::AbstractRNG = Random.GLOBAL_RNG
end

"""
    randomize(x, R::Shift)
Cranley Patterson Rotation i.e. `y = (x .+ U) mod 1` where `U ∼ 𝕌([0,1]ᵈ)` and `x` is a `d×n` matrix
"""
function randomize(x, R::Shift)
    y = copy(x)
    shift!(R.rng, y)
    return y
end

function randomize!(x, R::Shift)
    shift!(R.rng, x)
end

function shift!(rng::AbstractRNG, points::AbstractMatrix{T}) where {T <: Real}
    d = size(points, 1)
    U = zeros(T, d)
    shift!(rng, points, U)
end

function shift!(points::AbstractMatrix{T}, U::AbstractVector{T}) where {T <: Real}
    shift!(Random.default_rng(), points, U)
end

function shift!(rng::AbstractRNG, points::AbstractMatrix{T},
                U::AbstractVector{T}) where {T <: Real}
    rand!(rng, U)
    for i in axes(points, 2)
        points[:, i] += U
    end
    points[:] = frac.(points)
end

frac(y) = y - floor(y)

function generate_design_matrices(n, d, sampler, R::Shift, num_mats, T = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    out = sample(n, d, no_rand_sampler, T)

    # randomize (shift) num_mats times
    return [randomize!(out, R) for j in 1:num_mats]
end
