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
    U = zeros(T, d)
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

"""
    digital_shift(points::AbstractArray) 
Digital shift each corrdinate in base `b` is shifted as `yₖ = (xₖ + Uₖ) mod b` where `Uₖ ∼ 𝕌({0:b-1})`. `U` is the same for every point `points` but i.i.d along every dimensions.
"""
function digital_shift(rng::AbstractRNG, points::AbstractArray, b::Integer;
                       M = 32)
    random_points = permutedims(similar(points))
    digital_shift!(rng, random_points, permutedims(points), b; M = M)
    return permutedims(random_points)
end

function digital_shift(points::AbstractArray, b::Integer; M = 32)
    digital_shift(Random.default_rng(), points, b; M = M)
end

function digital_shift!(rng::AbstractRNG, random_points::AbstractArray{T, N}, points,
                        b::Integer; M = 32) where {T, N}
    bits = unif2bits(points, b, M = M)
    for s in axes(random_points, 3)
        digital_shift_bits!(rng, @view(bits[:, :, s]), b)
    end

    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(bits[:, i]), b)
    end
end

function digital_shift!(random_points::AbstractArray, points::AbstractArray, b::Integer;
                        M = 32)
    digital_shift!(Random.default_rng(), random_points, points, b; M = M)
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
