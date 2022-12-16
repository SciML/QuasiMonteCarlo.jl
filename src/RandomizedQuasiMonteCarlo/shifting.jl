"""
```julia
Shifting
```

Cranley-Patterson rotation aka Shifting
"""
struct Shifting <: RandMethod end

"""
```julia
DigitalShifting
```

Digital Shift
"""
struct DigitalShifting <: RandMethod end

"""
    shift(points::AbstractArray) 
Cranley Patterson Rotation i.e. `y = (points .+ U) mod 1` where `U ∼ 𝕌([0,1]ᵈ)` and `points` is a `n×d` matrix
"""
function shift(points::AbstractArray)
    y = copy(points)
    shift!(y)
    return y
end

function shift!(points::AbstractMatrix{<:Real})
    d = size(points, 2)
    s = zeros(d)
    shift!(points, s)
end

function shift!(points::AbstractMatrix{<:Real}, s::AbstractVector{<:Real})
    rand!(s)
    for i in axes(points, 1)
        points[i, :] += s
    end
    frac!(points)
end

function frac!(points::AbstractArray)
    for i in eachindex(points)
        points[i] -= floor(points[i])
    end
end

"""
    digital_shift(points::AbstractArray) 
Digital shift each corrdinate in base `b` is shifted as `yₖ = (xₖ + Uₖ) mod b` where `Uₖ ∼ 𝕌({0:b-1})`. `U` is the same for every point `points` but i.i.d along every dimensions.
"""
function digital_shift(rng::AbstractRNG, points::AbstractArray{T, N}, b::Integer;
                       M = 32) where {T, N}
    bits = unif2bits(points, b, M = M)
    random_points = copy(points)
    for s in axes(points, 3)
        digital_shift_bits!(rng, @view(bits[:, :, s]), b)
    end

    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(bits[:, i]), b)
    end

    return random_points
end

function digital_shift(points::AbstractArray, b::Integer; M = 32)
    digital_shift(Random.default_rng(), points, b; M = M)
end

function digital_shift_bits!(rng::AbstractRNG, random_bits::AbstractMatrix{<:Integer},
                             b::Integer)
    M, n = size(random_bits)
    DS = rand(rng, 0:(b - 1), M)
    for i in 1:n
        random_bits[:, i] = (random_bits[:, i] + DS) .% b
    end
end

function digital_shift_bits!(random_bits::AbstractMatrix{<:Integer}, b::Integer)
    digital_shift_bits!(Random.default_rng(), random_bits, b)
end
