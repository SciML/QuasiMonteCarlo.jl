"""
    GridSample(rng::AbstractRNG = Random.TaskLocalRNG()) <: RandomSamplingAlgorithm

The points of a rectangular grid in the unit box, in random order. The grid has `m` equal
cells along each axis, with `m` the smallest integer such that `m^d ≥ n`. Its points are
the `m^d` cell centers, so no coordinate is 0 or 1.

`sample(n, d, GridSample())` returns the centers of `n` distinct grid cells, drawn without
replacement and in random order. With `n = m^d` it is the whole grid. Every leading run of
points is a uniform sample of the grid, so its mean is an unbiased estimate of the grid's mean.
A narrow element type, such as `Float16` on a fine grid, can round two centers to one value.

# Fields

- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator that orders the grid.

In more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a
result, they should almost never be used for multivariate integration; their use is as
a starting point for other algorithms.

# Examples

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> sort(vec(sample(4, 1, GridSample(Xoshiro(1)))))
4-element Vector{Float64}:
 0.125
 0.375
 0.625
 0.875

julia> size(unique(sample(9, 2, GridSample()); dims = 2))
(2, 9)
```
"""
Base.@kwdef @concrete struct GridSample <: RandomSamplingAlgorithm
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function sample(n::Integer, d::Integer, S::GridSample, T = Float64)
    _check_sequence(n)
    d > 0 || throw(ArgumentError("GridSample needs at least one dimension, got d = $d"))
    cells = grid_cells(n, d)
    grid_size = big(cells)^d
    indices = distinct_indices(S.rng, grid_size <= typemax(Int) ? Int(grid_size) : grid_size, n)
    return fill_grid_points!(Matrix{T}(undef, d, n), indices, cells)
end

"""
    grid_cells(n, d)

The smallest number `m` of cells per axis such that the grid of `m^d` cells holds `n` points.

```jldoctest
julia> using QuasiMonteCarlo

julia> QuasiMonteCarlo.grid_cells(10, 2)
4
```
"""
function grid_cells(n::Integer, d::Integer)
    # The float root is near `m`; the exact checks correct it.
    cells = max(big(1), floor(BigInt, Float64(n)^(1 / d)))
    while cells^d < n
        cells += 1
    end
    while cells > 1 && (cells - 1)^d >= n
        cells -= 1
    end
    return Int(cells)
end

"""
    fill_grid_points!(points, indices, cells)

Write grid point `indices[j]` into column `j` of `points`. Coordinate `i` is the center of the
cell given by digit `i` of the index in base `cells`.

```jldoctest
julia> using QuasiMonteCarlo

julia> QuasiMonteCarlo.fill_grid_points!(zeros(2, 2), [1, 2], 2)
2×2 Matrix{Float64}:
 0.75  0.25
 0.25  0.75
```
"""
function fill_grid_points!(points::AbstractMatrix{T}, indices::AbstractVector{<:Integer}, cells::Integer) where {T}
    for (point, index) in zip(eachcol(points), indices)
        rest = index
        for axis in eachindex(point)
            rest, cell = divrem(rest, cells)
            point[axis] = grid_center(T, Int(cell), cells)
        end
    end
    return points
end

"""
    distinct_indices(rng, count, n)

`n` distinct indices of `0:(count - 1)`, drawn uniformly without replacement, in random
order: every leading run of them is a uniform sample without replacement.

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> sort(QuasiMonteCarlo.distinct_indices(Xoshiro(1), 5, 5))
5-element Vector{Int64}:
 0
 1
 2
 3
 4
```
"""
function distinct_indices(rng::AbstractRNG, count::Integer, n::Integer)
    if !(0 <= n <= count)
        throw(ArgumentError("cannot draw $n distinct indices from $count"))
    end
    # A shuffle needs `count` entries; rejection needs fewer than 2n draws when n ≤ count / 2.
    if 2n >= count
        return shuffle(rng, 0:(count - 1))[1:n]
    end
    drawn = Set{typeof(count)}()
    indices = typeof(count)[]
    while length(indices) < n
        index = rand(rng, 0:(count - 1))
        if !(index in drawn)
            push!(drawn, index)
            push!(indices, index)
        end
    end
    return indices
end

"""
    grid_center(T, cell, n)

The center `(2 cell + 1) / 2n` of cell `cell` of `n` equal cells of `[0, 1]`, as a `T`. A float
keeps it strictly inside `(0, 1)`, also when its rounding would reach 0 or 1.

```jldoctest
julia> using QuasiMonteCarlo

julia> QuasiMonteCarlo.grid_center(Float64, 3, 4)
0.875
```
"""
function grid_center(::Type{T}, cell::Integer, n::Integer) where {T <: AbstractFloat}
    return clamp(T(_exact_center(cell, n)), nextfloat(zero(T)), prevfloat(one(T)))
end
grid_center(::Type{T}, cell::Integer, n::Integer) where {T} = T(_exact_center(cell, n))
"""
    _exact_center(cell, n)

The exact center `(2 cell + 1) / 2n` as a `Rational` of widened integers, so `2n` cannot overflow.

```jldoctest
julia> using QuasiMonteCarlo

julia> QuasiMonteCarlo._exact_center(Int32(0), Int32(2)^30)
1//2147483648
```
"""
_exact_center(cell::Integer, n::Integer) = (2widen(cell) + 1) // (2widen(n))

"""
    LatticeRuleSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

Generate a point set using a rank-1 lattice rule. Each coordinate is moved up by half the
last binary digit that the output type holds, so the first point is not the origin.

# Fields

- `R::RandomizationMethod = NoRand()`: Randomization applied to the lattice
  points.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> points = sample(8, 2, LatticeRuleSample());

julia> size(points)
(2, 8)
```
"""
Base.@kwdef @concrete struct LatticeRuleSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::LatticeRuleSample, T = Float64)
    lat = _lattice_rule(d)
    result = reduce(hcat, lat[0:(n - 1)])
    return randomize(cell_midpoint.(T, result), S.R)
end

"""
    cell_midpoint(T, x)

The midpoint of the cell of `precision(T) - 1` binary digits that holds `x ∈ [0, 1)`, as a
`T`: `x` moved up by at most one last digit of `T`, strictly inside `(0, 1)`. A type other
than a float takes the `Float64` midpoint.

```jldoctest
julia> using QuasiMonteCarlo

julia> QuasiMonteCarlo.cell_midpoint(Float64, 0.0) == 2.0^-53
true
```
"""
function cell_midpoint(::Type{T}, x::AbstractFloat) where {T <: AbstractFloat}
    kept = precision(T) - 1
    # Scale in the wider type, so neither the scaling nor the conversion rounds or overflows.
    wide = precision(T) >= precision(x) ? T(x) : x
    return ldexp(T(floor(ldexp(wide, kept))) + one(T) / 2, -kept)
end
cell_midpoint(::Type{T}, x::AbstractFloat) where {T} = T(cell_midpoint(Float64, x))
