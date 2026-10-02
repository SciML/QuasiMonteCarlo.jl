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

"""
    OwenScramble(base::Integer; pad = 32, rng = Random.TaskLocalRNG()) <: ScrambleMethod

Nested Uniform Scramble, also known as Owen's scramble.

# Fields

- `base::Integer`: Base of the digital net being scrambled.
- `pad::Integer = 32`: Number of base-`base` digits retained for each point.
- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used by
  the scramble.

`randomize(x, R::OwenScramble)` returns a scrambled version of `x`.
The scramble method is Nested Uniform Scramble which was introduced in Owen (1995).
`pad` is the number of bits used for each point. One needs `pad ≥ log(base, n)`.

References: Owen, A. B. (1995). Randomly permuted (t, m, s)-nets and (t, s)-sequences. In Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing: Proceedings of a conference at the University of Nevada, Las Vegas, Nevada, USA, June 23–25, 1994 (pp. 299-317). Springer New York.
"""
Base.@kwdef struct OwenScramble{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.TaskLocalRNG()
end

function randomize!(
        random_points::AbstractMatrix{T},
        points::AbstractMatrix{T}, R::OwenScramble
    ) where {T <: Real}
    @assert size(points) == size(random_points)
    b = R.base
    unrandomized_bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(unrandomized_bits)
    indices = which_permutation(unrandomized_bits, b)
    randomize_bits!(random_bits, unrandomized_bits, indices, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
    return
end

"""
    randomize_bits!(random_bits::AbstractArray{T, 3}, origin_bits::AbstractArray{T, 3}, R::ScrambleMethod) where {T <: Integer}

In place version of a `OwenScramble` (Nested Uniform Scramble) for the "bit" array.
This is faster to use this functions for multiple scramble of the same array (use `generate_design_matrices`).
"""
function randomize_bits!(
        random_bits::AbstractArray{T, 3},
        origin_bits::AbstractArray{T, 3},
        indices::AbstractArray{F, 3},
        R::OwenScramble
    ) where {T <: Integer, F <: Integer}
    # in place nested uniform Scramble.
    #
    m, n, d = size(indices)
    b = R.base
    rng = R.rng
    pad = size(random_bits, 1)
    @assert m ≥ 1 "We need m ≥ 1" # m=0 causes awkward corner case below.

    for s in 1:d
        theperms = getpermset(rng, m, b)         # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            @views random_bits[k, :, s] .= (
                origin_bits[k, :, s] .+
                    theperms[k, indices[k, :, s]]
            ) .%
                b   # permutation by adding a bit modulo b
        end
    end

    # Paste in random entries for bits after m'th one
    return if pad > m
        # random_bits[(m + 1):pad, :, :] = rand(rng, 0:(b - 1), n * d * (pad - m))
        rand!(rng, @view(random_bits[(m + 1):pad, :, :]), 0:(b - 1))
    end
end

function getpermset(rng::AbstractRNG, m::Integer, b::I) where {I <: Integer}
    # Get b^(k-1) random binary permutations for k=1 ... m
    # m will ordinarily be m when there are n=b^m points
    #
    y = zeros(I, m, b^(m - 1))
    for k in 1:m
        nₖ = b^(k - 1)
        y[k, 1:nₖ] = rand(rng, 0:(b - 1), nₖ)
    end
    return y
end

getpermset(m::Integer, b::Integer) = getpermset(Random.TaskLocalRNG(), m, b)

"""
    which_permutation(bits::AbstractArray{<:Integer,3}, b)

This function is used in Nested Uniform Scramble.
It assigns for each point (in every dimension) `m` number corresponding to its location on the slices 1, 1/b, 1/b², ..., 1/bᵐ⁻¹ of the axes [0,1[.
This also can be used to verify some equidistribution prorepreties.
Here we create the `indices` array `m`, and not `pad`. Indeed `(t,m,d)-net` in base `b` are scrambled up to the `1/bᵐ` component.
Higher order components are just used i.i.d. `Uₖ ∼ 𝐔({0:b-1})` in `owen_scramble_bit!`.
"""
function which_permutation(bits::AbstractArray, b::I) where {I <: Integer}
    n, d = size(bits)[2:end]
    m = logi(b, n)

    indices = zeros(I, m, n, d)
    for j in axes(bits, 3)
        which_permutation!(@view(indices[:, :, j]), bits[:, :, j], b)
    end
    return indices
end

function which_permutation!(
        indices::AbstractMatrix{<:Integer},
        bits::AbstractMatrix{<:Integer}, b::Integer
    )
    @assert size(indices)[2:end] == size(bits)[2:end] "You need size(indices) = $(size(indices)) equal size(bits) = $(size(bits))"
    indices[1, :] .= 0 # same permutation for all observations i
    for i in axes(indices, 2)                     # Here is where we want m > 0 so the loop works ok
        for k in axes(indices, 1)[2:end]
            indices[k, i] = bits2int(@view(bits[1:(k - 1), i]), b) # index of which perms to use at bit k for each i
        end
    end
    return indices .+= 1 # array indexing starts at 1
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
    MatousekScramble(base::Integer; pad = 32, rng = Random.TaskLocalRNG()) <: ScrambleMethod

Linear Matrix Scramble, also known as Matousek's scramble.

# Fields

- `base::Integer`: Base of the digital net being scrambled.
- `pad::Integer = 32`: Number of base-`base` digits retained for each point.
- `rng::AbstractRNG = Random.TaskLocalRNG()`: Random-number generator used by
  the scramble.

`randomize(x, R::MatousekScramble)` returns a scrambled version of `x`.
The scramble method is Linear Matrix Scramble which was introduced in Matousek (1998):
in each dimension the `pad` digits `a` of a point become `M a + c` modulo `base`, where `M`
is a random lower-triangular `pad × pad` matrix with non-zero diagonal and `c` a random
digit vector. The matrix and shift are drawn from `rng` in dimension order, whatever the
number of points, so a longer sample keeps the points of a shorter one and the first `k`
dimensions of a `d`-dimensional sample are the `k`-dimensional sample's. In base 2 the
scrambled points of a digital net are a digital net (with generating matrices `M C`)
shifted by `c`. `pad` is the number of digits used for each point. A [`DigitalNetSample`](@ref)
applies the scramble to its generating matrices, which costs `O(pad²)` per column rather than
per point ([`scramble_generators`](@ref)).

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
    pad, _, d = size(origin_bits)
    for s in 1:d
        # A pad × pad matrix and shift per dimension, drawn in dimension order.
        matousek_M, matousek_C = getmatousek(R.rng, pad, R.base)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # matousek_M (pad×pad) * origin_bits (pad×n) .+ matousek_C (pad×1)
        @views random_bits[:, :, s] .= (matousek_M * origin_bits[:, :, s] .+ matousek_C) .% R.base
    end
    return random_bits
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
It scrambles each coordinate in base `b` as `yₖ = (xₖ + Uₖ) mod b` where `Uₖ ∼ 𝕌({0:b-1})`
for each of the `pad` digits `k`. `U` is the same for every point `points` but i.i.d. along
every dimension, and drawn from `rng` in dimension order, so a point's scramble depends on
the point and the first draws of `rng` alone: it does not depend on the number of points
and the first `k` dimensions of a `d`-dimensional sample are the `k`-dimensional sample's.
On the generating matrices of a [`DigitalNetSample`](@ref) it is applied to the matrices.
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
    pad, _, d = size(origin_bits)
    for s in 1:d
        # One shift of all `pad` digits per dimension, drawn in dimension order.
        digit_shift = rand(R.rng, 0:(R.base - 1), pad)
        @views random_bits[:, :, s] .= (origin_bits[:, :, s] .+ digit_shift) .% R.base
    end
    return random_bits
end

const DigitalMatrixScramble = Union{MatousekScramble, DigitalShift}

"""
    scramble_generators(generating_matrices::AbstractMatrix{U}, d::Integer, R::DigitalMatrixScramble)

The generating matrices of the first `d` dimensions of a base-2 digital net, scrambled by
`R`, and the digit shift to add to every point. In dimension `s`, with `pad` the number of
digits, the matrix `L` and the shift `c` that `randomize_bits!` draws from `R.rng`, in
dimension order, give the columns `L C` and the shift `c`. The points of the scrambled net
are the XORs of its columns selected by the digits of the point's index, XORed with `c`, so
they equal the scramble of the unscrambled points. Only the first `pad` digits of each
column and shift are kept. It costs `O(pad²)` per column, not per point.

The generating matrices are left-aligned `U` words; returns the scrambled matrices and a
vector of shifts, one per dimension.
"""
function scramble_generators(
        generating_matrices::AbstractMatrix{U}, d::Integer, R::DigitalMatrixScramble
    ) where {U <: Unsigned}
    width = 8 * sizeof(U)
    pad = Int(R.pad)
    if R.base != 2
        throw(ArgumentError("digital nets are base 2, but the scramble has base $(R.base)"))
    end
    if !(1 <= pad <= width)
        throw(ArgumentError("pad = $pad digits do not fit in a $width-bit word"))
    end
    # Dimension order, so the first `k` dimensions are the same for every `d ≥ k`.
    draws = [draw_digit_words(R.rng, R, U) for _ in 1:d]
    scrambled = stack(
        [
            [multiply_digits(masks, word) for word in generators]
                for ((masks, _), generators) in zip(draws, eachrow(@view generating_matrices[1:d, :]))
        ]; dims = 1
    )
    return scrambled, last.(draws)
end

@public scramble_generators

"""
    draw_digit_words(rng, R::DigitalMatrixScramble, U)

One dimension's scramble of `R` as left-aligned words of type `U`: the rows of the
lower-triangular matrix, as bit masks over the digits (row `k` is a mask whose parity with
a word gives digit `k` of the product), and the shift.
"""
function draw_digit_words(rng::AbstractRNG, R::MatousekScramble, ::Type{U}) where {U <: Unsigned}
    matrix, shift = getmatousek(rng, Int(R.pad), R.base)
    masks = [digits_to_word(U, row) for row in eachrow(matrix)]
    return masks, digits_to_word(U, shift)
end

function draw_digit_words(rng::AbstractRNG, R::DigitalShift, ::Type{U}) where {U <: Unsigned}
    shift = rand(rng, 0:(R.base - 1), Int(R.pad))
    identity_masks = [one(U) << (8 * sizeof(U) - digit) for digit in 1:Int(R.pad)]
    return identity_masks, digits_to_word(U, shift)
end

"""The left-aligned word of type `U` whose leading binary digits are `digits`."""
function digits_to_word(::Type{U}, digits::AbstractVector{<:Integer}) where {U <: Unsigned}
    word = zero(U)
    for (position, digit) in enumerate(digits)
        word |= U(digit) << (8 * sizeof(U) - position)
    end
    return word
end

"""Digit `k` of the product is the parity of `masks[k]` ANDed with `word`: a matrix-vector product over GF(2)."""
function multiply_digits(masks::AbstractVector{U}, word::U) where {U <: Unsigned}
    product = zero(U)
    for (position, mask) in enumerate(masks)
        product |= U(count_ones(mask & word) & 1) << (8 * sizeof(U) - position)
    end
    return product
end
