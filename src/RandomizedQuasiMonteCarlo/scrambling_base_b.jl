# * The scrambling codes were first inspired from Owen's `R` implementation that can be found [here](https://artowen.su.domains/code/rsobol.R). * #
#TODO the typing could probably improved, resulting in more type stable (important for QMC computation to be able to use Float32, Float64, Rational, ... for points, and Int64, Int32, Int16 etc for bits etc.)

"""
    ScrambleMethod <: RandomizationMethod

A scramble method randomizes a digital net in an integer base.

# Interface rules

Concrete subtypes must provide the fields `base::Integer`, `pad::Integer`,
and `rng::AbstractRNG`. `base` must match the base of the sequence being
scrambled, and `pad` must be at least `log(base, n)` for a point set with `n`
points. The `randomize` interface preserves the matrix size and keeps points in
the unit box.

The package provides [`DigitalShift`](@ref), [`MatousekScramble`](@ref), and
[`OwenScramble`](@ref). New implementations should add methods to the generic
`randomize`/`randomize!` interface rather than changing sampler methods.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> points = sample(8, 2, FaureSample());

julia> randomized = randomize(points, DigitalShift(base = 2, pad = 4));

julia> size(randomized) == size(points)
true
```

The scramble methods implementer are

  - `DigitalShift`.
  - `OwenScramble`: Nested Uniform Scramble which was introduced in Owen (1995).
  - `MatousekScramble`: Linear Matrix Scramble which was introduced in Matousek (1998).
"""
abstract type ScrambleMethod <: RandomizationMethod end

function randomize(x, R::ScrambleMethod)
    random_x = permutedims(copy(x))
    randomize!(random_x, permutedims(x), R)
    return permutedims(random_x)
end

function randomize!(
        random_points::AbstractMatrix{T},
        points::AbstractMatrix{T}, R::ScrambleMethod
    ) where {T <: Real}
    @assert size(points) == size(random_points)
    b = R.base
    unrandomized_bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(unrandomized_bits)
    randomize_bits!(random_bits, unrandomized_bits, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
    return
end

"""
    OwenScramble(; base::Integer, pad = 32, rng = Random.TaskLocalRNG()) <: ScrambleMethod

Nested Uniform Scramble, also known as Owen's scramble.

# Fields

- `base::Integer`: Prime base of the digital net being scrambled.
- `pad::Integer = 32`: Number of base-`base` digits scrambled in each coordinate.
- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator from which each
  `randomize` call draws one 64-bit seed.

A nested scramble gives each node of the `base`-ary digit tree a random permutation of `0:(base - 1)`.
Digit `k` of a coordinate goes through the permutation of the node that the `k - 1` digits above it select.
A hash of the seed, the dimension and these digits gives the permutation (Burley 2020), so no table is stored.
The scramble of a point thus depends only on the point and the seed:

  - `randomize` accepts any number of points.
  - The first `n` points of a scrambled sequence are the scrambled first `n` points.
    Extending a sample keeps the points already drawn.

Each node permutation is a random affine map `a ↦ h·a + g mod base`, with `h` uniform on `1:(base - 1)` and `g` uniform on `0:(base - 1)`.
For `base` 2 or 3, this is a uniform random permutation.
For larger bases, each pair of distinct digits maps to a uniform pair of distinct digits, as under a uniform permutation.
Two scrambled points thus have the same joint law as under Owen's scramble, so estimates have the same variance.
The cost is `O(pad)` hashes per coordinate.

References: Owen, A. B. (1995). Randomly permuted (t, m, s)-nets and (t, s)-sequences. In Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing: Proceedings of a conference at the University of Nevada, Las Vegas, Nevada, USA, June 23–25, 1994 (pp. 299-317). Springer New York.
Burley, B. (2020). Practical Hash-based Owen Scrambling. Journal of Computer Graphics Techniques, 9(4), 1-20.

# Examples

```jldoctest
julia> using QuasiMonteCarlo, Random

julia> points = sample(8, 2, SobolSample());

julia> first_three = randomize(points[:, 1:3], OwenScramble(base = 2, rng = Xoshiro(1)));

julia> randomize(points, OwenScramble(base = 2, rng = Xoshiro(1)))[:, 1:3] == first_three
true
```
"""
Base.@kwdef struct OwenScramble{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.TaskLocalRNG()
    function OwenScramble(base::I, pad::I, rng::AbstractRNG) where {I <: Integer}
        # A multiplier `h` that shares a factor with `base` would merge digits.
        isprime(base) || throw(ArgumentError("OwenScramble needs a prime base, got $base"))
        return new{I}(base, pad, rng)
    end
end

"""Odd 64-bit constant `2^64 / φ` that SplitMix64 adds to its state."""
const GOLDEN_GAMMA = 0x9e3779b97f4a7c15

"""
    splitmix64(word::Unsigned) -> UInt64

Return the SplitMix64 output for state `word`: a bijective avalanche of 64 bits.

# Examples

```jldoctest
julia> QuasiMonteCarlo.splitmix64(0x0)
0xe220a8397b1dcdaf
```
"""
function splitmix64(word::Unsigned)
    mixed = UInt64(word) + GOLDEN_GAMMA
    mixed = (mixed ⊻ (mixed >> 30)) * 0xbf58476d1ce4e5b9
    mixed = (mixed ⊻ (mixed >> 27)) * 0x94d049bb133111eb
    return mixed ⊻ (mixed >> 31)
end

"""
    root_node(seed::UInt64, dimension::Integer) -> UInt64

Return the hash of the digit tree's root for coordinate `dimension` under `seed`.

# Examples

```jldoctest
julia> QuasiMonteCarlo.root_node(UInt64(0), 1) != QuasiMonteCarlo.root_node(UInt64(0), 2)
true
```
"""
root_node(seed::UInt64, dimension::Integer) = splitmix64(seed ⊻ splitmix64(UInt64(dimension)))

"""
    child_node(node::UInt64, digit::Integer) -> UInt64

Return the hash of the child that `digit` selects below `node` in the digit tree.
Distinct digits give distinct children, as `splitmix64` is a bijection.

# Examples

```jldoctest
julia> QuasiMonteCarlo.child_node(0x0, 0) != QuasiMonteCarlo.child_node(0x0, 1)
true
```
"""
child_node(node::UInt64, digit::Integer) = splitmix64(node + (UInt64(digit) + 1) * GOLDEN_GAMMA)

"""
    node_shift(node::UInt64, base::Integer) -> UInt64

Return the shift in `0:(base - 1)` that `node` applies to its digit: the top 64 bits of `node * base`.

# Examples

```jldoctest
julia> QuasiMonteCarlo.node_shift(typemax(UInt64), 3)
0x0000000000000002
```
"""
node_shift(node::UInt64, base::Integer) = (widemul(node, UInt64(base)) >> 64) % UInt64

"""
    node_multiplier(node::UInt64, base::Integer) -> UInt64

Return the multiplier in `1:(base - 1)` that `node` applies to its digit, from one more hash of `node`.

# Examples

```jldoctest
julia> QuasiMonteCarlo.node_multiplier(UInt64(12345), 2)
0x0000000000000001
```
"""
node_multiplier(node::UInt64, base::Integer) = node_shift(splitmix64(node), base - 1) + 1

"""
    add_mod(digit::T, shift::T, base::T) where {T <: Unsigned} -> T

Return `(digit + shift) mod base` without overflow, for `digit, shift < base`.

# Examples

```jldoctest
julia> QuasiMonteCarlo.add_mod(typemax(UInt64) - 1, typemax(UInt64) - 1, typemax(UInt64))
0xfffffffffffffffd
```
"""
function add_mod(digit::T, shift::T, base::T) where {T <: Unsigned}
    gap = base - digit
    return shift >= gap ? shift - gap : digit + shift
end

"""
    permute_digit(node::UInt64, digit::UInt64, base::UInt64) -> UInt64

Return `h·digit + g mod base`, with multiplier `h` and shift `g` from the hash `node`.
Throw a `DomainError` if `digit` is not a base-`base` digit.

# Examples

```jldoctest
julia> sort([QuasiMonteCarlo.permute_digit(UInt64(7), UInt64(digit), UInt64(5)) for digit in 0:4])
5-element Vector{UInt64}:
 0x0000000000000000
 0x0000000000000001
 0x0000000000000002
 0x0000000000000003
 0x0000000000000004
```
"""
function permute_digit(node::UInt64, digit::UInt64, base::UInt64)
    digit < base || throw(DomainError(digit, "a base-$base digit must lie in 0:$(base - 1)"))
    multiplier = node_multiplier(node, base)
    # Every base-2 multiplier is 1, so base 2 skips the 128-bit division.
    scaled = isone(multiplier) ? digit : (widemul(multiplier, digit) % base) % UInt64
    return add_mod(scaled, node_shift(node, base), base)
end

function randomize_bits!(
        random_bits::AbstractArray{T, 3},
        origin_bits::AbstractArray{T, 3},
        R::OwenScramble
    ) where {T <: Integer}
    pad, n, d = size(origin_bits)
    seed = rand(R.rng, UInt64)
    for dimension in 1:d
        root = root_node(seed, dimension)
        for point in 1:n
            node = root
            for k in 1:pad
                digit = origin_bits[k, point, dimension]
                random_bits[k, point, dimension] = T(permute_digit(node, UInt64(digit), UInt64(R.base)))
                node = child_node(node, digit)
            end
        end
    end
    return random_bits
end

"""
    MatousekScramble(base::Integer; pad = 32, rng = Random.TaskLocalRNG()) <: ScrambleMethod

Linear Matrix Scramble, also known as Matousek's scramble.

# Fields

- `base::Integer`: Base of the digital net being scrambled.
- `pad::Integer = 32`: Number of base-`base` digits retained for each point.
- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used by
  the scramble.

`randomize(x, R::MatousekScramble)` returns a scrambled version of `x`.
The scramble method is Linear Matrix Scramble which was introduced in Matousek (1998).
`pad` is the number of bits used for each point. One needs `pad ≥ log(base, n)`.

References: Matoušek, J. (1998). On thel2-discrepancy for anchored boxes. Journal of Complexity, 14(4), 527-556.
"""
Base.@kwdef struct MatousekScramble{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.TaskLocalRNG()
end

#? Weird it should be faster than nested uniform Scramble but here it is not at all.-> look for other implementation and paper
"""
    randomize_bits!(random_bits::AbstractArray{T, 3}, origin_bits::AbstractArray{T, 3}, R::ScrambleMethod) where {T <: Integer}

In place version of a ScrambleMethod (`MatousekScramble` or `DigitalShift`) for the "bit" array.
This is faster to use this functions for multiple scramble of the same array (use `generate_design_matrices`).
"""
function randomize_bits!(
        random_bits::AbstractArray{T, 3},
        origin_bits::AbstractArray{T, 3},
        R::MatousekScramble
    ) where {T <: Integer}
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    pad, n, d = size(origin_bits)
    b = R.base
    rng = R.rng
    m = logi(b, n)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek(rng, m, b)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # matousek_M (m×m) * origin_bits (m×n) .+ matousek_C (m×1)
        @views random_bits[1:m, :, s] .= (
            matousek_M * origin_bits[1:m, :, s] .+
                matousek_C
        ) .% b
    end

    # Paste in random entries for bits after m'th one
    return if pad > m
        # random_bits[(m + 1):pad, :, :] = rand(rng, 0:(b - 1), n * d * (pad - m))
        rand!(rng, @view(random_bits[(m + 1):pad, :, :]), 0:(b - 1))
    end
end

"""
    getmatousek(rng::AbstractRNG, m::Integer, b::Integer)

Generate the Matousek linear scramble in base b for one of the d components
It produces a m x m bit matrix matousek_M and a length m bit vector matousek_C
"""
function getmatousek(rng::AbstractRNG, m::Integer, b::I) where {I <: Integer}
    matousek_M = LowerTriangular(zeros(I, m, m)) + Diagonal(rand(rng, 1:(b - 1), m)) # Mₖₖ ∼ U{1, ⋯, b-1}
    matousek_C = rand(rng, 0:(b - 1), m)
    for i in 2:m
        for j in 1:(i - 1)
            matousek_M[i, j] = rand(rng, 0:(b - 1))
        end
    end
    return matousek_M, matousek_C
end

getmatousek(m::Integer, b::Integer) = getmatousek(Random.TaskLocalRNG(), m, b)

"""
    DigitalShift(base::Integer; pad = 32, rng = Random.TaskLocalRNG()) <: ScrambleMethod

Digital shift.
`randomize(x, R::DigitalShift)` returns a scrambled version of `x`.

# Fields

- `base::Integer`: Base of the digital net being scrambled.
- `pad::Integer = 32`: Number of base-`base` digits retained for each point.
- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used by
  the shift.

The scramble method is Digital Shift.
It scrambles each coordinate in base `b` as `yₖ = (xₖ + Uₖ) mod b` where `Uₖ ∼ 𝕌({0:b-1})`.
`U` is the same for every point `points` but i.i.d. along every dimension.
"""
Base.@kwdef struct DigitalShift{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function randomize_bits!(
        random_bits::AbstractArray{T, 3},
        origin_bits::AbstractArray{T, 3},
        R::DigitalShift
    ) where {T <: Integer}
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    pad, n, d = size(origin_bits)
    b = R.base
    rng = R.rng
    m = logi(b, n)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        DS = rand(rng, 0:(b - 1), m)

        # xₖ = (aₖ + Cₖ) mod b where xₖ is the k element in base b
        # origin_bits (m×n) .+ DS (m×1)
        @views random_bits[1:m, :, s] .= (origin_bits[1:m, :, s] .+ DS) .% b
    end

    # Paste in random entries for bits after m'th one
    return if pad > m
        # random_bits[(m + 1):pad, :, :] = rand(rng, 0:(b - 1), n * d * (pad - m))
        rand!(rng, @view(random_bits[(m + 1):pad, :, :]), 0:(b - 1))
    end
end
