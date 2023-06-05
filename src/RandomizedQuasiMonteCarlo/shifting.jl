"""
```julia
Shifting <: RandomizationMethod
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

"""
```julia
generate_design_matrices(n, d, sample_method::DeterministicSamplingAlgorithm,
                         R::RandomizationMethod, num_mats, T = Float64)
```
Create `num_mats` matrices each containing a QMC point set in `[0, 1)ᵈ`, where:
- `n` is the number of points to sample.
- `d` is the dimensionality of the point set.
- `sample_method` is the quasi-Monte Carlo sampling strategy used to create a deterministic point set `out`.
- `R` is the method used to randomize `num_mats` times the point set `out`. Each randomization is i.i.d.
- `T` is the `eltype` of the point sets. For some QMC methods (Faure, Sobol) this can be `Rational`
"""
function generate_design_matrices(n, d, sampler, R::Shift, num_mats, T = Float64)
    # Generate unrandomized sequence
    no_rand_sampler = @set sampler.R = NoRand()
    out = sample(n, d, no_rand_sampler, T)

    # randomize (shift) num_mats times
    return [randomize!(out, R) for j in 1:num_mats]
end
