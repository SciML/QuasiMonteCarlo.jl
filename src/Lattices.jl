"""
    GridSample(R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm

A simple rectangular grid lattice. It samples from a regular grid in the unit
box and optionally applies `R`. Each coordinate is the center of one of `n` equal
cells, so no coordinate is 0 or 1.

# Fields

- `R::RandomizationMethod = NoRand()`: Randomization applied to the grid.

In more than 2 dimensions, grids have worse discrepancy than simple Monte Carlo. As a
result, they should almost never be used for multivariate integration; their use is as
a starting point for other algorithms.
"""
Base.@kwdef @concrete struct GridSample <: DeterministicSamplingAlgorithm
    R::RandomizationMethod = NoRand()
end

function sample(n::Integer, d::Integer, S::GridSample, T = Float64)
    centers = [grid_center(T, cell, n) for cell in 0:(n - 1)]
    samples = [rand(centers, n) for _ in 1:d]
    return randomize(mapreduce(permutedims, vcat, samples), S.R)
end

"""
    grid_center(T, cell, n)

The center `(2 cell + 1) / 2n` of cell `cell` of `n` equal cells of `[0, 1]`, as a `T`. A float
keeps it strictly inside `(0, 1)`, also when its rounding would reach 0 or 1.

```jldoctest
julia> QuasiMonteCarlo.grid_center(Float64, 3, 4)
0.875
```
"""
function grid_center(::Type{T}, cell::Integer, n::Integer) where {T <: AbstractFloat}
    return clamp(T((2cell + 1) // 2n), nextfloat(zero(T)), prevfloat(one(T)))
end
grid_center(::Type{T}, cell::Integer, n::Integer) where {T} = T((2cell + 1) // 2n)

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
