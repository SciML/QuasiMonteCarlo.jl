# * The scrambling codes were first inspired from Owen's `R` implementation that can be found [here](https://artowen.su.domains/code/rsobol.R). * #
#TODO the typing could probably improved, resulting in more type stable (important for QMC computation to be able to use Float32, Float64, Rational, ... for points, and Int64, Int32, Int16 etc for bits etc.)

"""
```julia
ScrambleMethod <: RandomizationMethod
```

A scramble method needs at least the scrambling base `b`, the number of "bits" to use `pad` (`pad=32` is the default) and a seed `rng` (`rng = Random.GLOBAL_RNG` is the default).
The scramble methods implementer are

  - `DigitalShift`.
  - `OwenScramble`: Nested Uniform Scramble which was introduced in Owen (1995).
  - `MatousekScramble`: Linear Matrix Scramble which was introduced in Matousek (1998).
"""
abstract type ScrambleMethod <: RandomizationMethod end

"""
    randomize(x, R::ScrambleMethod)

Return a scrambled version of `x`.
The scramble methods implemented are

  - `DigitalShift`.
  - `OwenScramble`: Nested Uniform Scramble which was introduced in Owen (1995).
  - `MatousekScramble`: Linear Matrix Scramble which was introduced in Matousek (1998).
"""
function randomize(x, R::ScrambleMethod)
    random_x = permutedims(copy(x))
    randomize!(random_x, permutedims(x), R)
    return permutedims(random_x)
end

"""
```julia
OwenScramble <: ScrambleMethod
```

Nested Uniform Scramble aka Owen's scramble.

`randomize(x, R::OwenScramble)` returns a scrambled version of `x`.
The scramble method is Nested Uniform Scramble which was introduced in Owen (1995).
`pad` is the number of bits used for each point. One needs `pad ≥ log(base, n)`.

References: Owen, A. B. (1995). Randomly permuted (t, m, s)-nets and (t, s)-sequences. In Monte Carlo and Quasi-Monte Carlo Methods in Scientific Computing: Proceedings of a conference at the University of Nevada, Las Vegas, Nevada, USA, June 23–25, 1994 (pp. 299-317). Springer New York.
"""
Base.@kwdef struct OwenScramble{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function randomize!(random_points::AbstractMatrix{T},
    points::AbstractMatrix{T}, R::OwenScramble) where {T <: Real}
    @assert size(points) == size(random_points)
    b = R.base
    unrandomized_bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(unrandomized_bits)
    indices = which_permutation(unrandomized_bits, b)
    randomize_bits!(random_bits, unrandomized_bits, indices, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
end

"""
    randomize_bits!(random_bits::AbstractArray{T, 3}, origin_bits::AbstractArray{T, 3}, R::ScrambleMethod) where {T <: Integer}

In place version of a `OwenScramble` (Nested Uniform Scramble) for the "bit" array.
This is faster to use this functions for multiple scramble of the same array (use `generate_design_matrices`).
"""
function randomize_bits!(random_bits::AbstractArray{T, 3},
    origin_bits::AbstractArray{T, 3},
    indices::AbstractArray{F, 3},
    R::OwenScramble) where {T <: Integer, F <: Integer}
    # in place nested uniform Scramble.
    #
    m, n, d = size(indices)
    b = R.base
    rng = R.rng
    pad = size(random_bits, 1)
    @assert m≥1 "We need m ≥ 1" # m=0 causes awkward corner case below.

    for s in 1:d
        theperms = getpermset(rng, m, b)         # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            @views random_bits[k, :, s] .= (origin_bits[k, :, s] .+
                                            theperms[k, indices[k, :, s]]) .%
                                           b   # permutation by adding a bit modulo b
        end
    end

    # Paste in random entries for bits after m'th one
    if pad > m
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

getpermset(m::Integer, b::Integer) = getpermset(Random.GLOBAL_RNG, m, b)

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

function which_permutation!(indices::AbstractMatrix{<:Integer},
    bits::AbstractMatrix{<:Integer}, b::Integer)
    @assert size(indices)[2:end]==size(bits)[2:end] "You need size(indices) = $(size(indices)) equal size(bits) = $(size(bits))"
    indices[1, :] .= 0 # same permutation for all observations i
    for i in axes(indices, 2)                     # Here is where we want m > 0 so the loop works ok
        for k in axes(indices, 1)[2:end]
            indices[k, i] = bits2int(@view(bits[1:(k - 1), i]), b) # index of which perms to use at bit k for each i
        end
    end
    indices .+= 1 # array indexing starts at 1
end

function randomize!(random_points::AbstractMatrix{T},
    points::AbstractMatrix{T}, R::ScrambleMethod) where {T <: Real}
    @assert size(points) == size(random_points)
    b = R.base
    unrandomized_bits = unif2bits(points, b, pad = R.pad)
    random_bits = similar(unrandomized_bits)
    randomize_bits!(random_bits, unrandomized_bits, R)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
end

"""
```julia
MatousekScramble <: ScrambleMethod
```

Linear Matrix Scramble aka Matousek's scramble.

`randomize(x, R::MatousekScramble)` returns a scrambled version of `x`.
The scramble method is Linear Matrix Scramble which was introduced in Matousek (1998).
`pad` is the number of bits used for each point. One needs `pad ≥ log(base, n)`.

References: Matoušek, J. (1998). On thel2-discrepancy for anchored boxes. Journal of Complexity, 14(4), 527-556.
"""
Base.@kwdef struct MatousekScramble{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.GLOBAL_RNG
end

#? Weird it should be faster than nested uniform Scramble but here it is not at all.-> look for other implementation and paper
"""
    randomize_bits!(random_bits::AbstractArray{T, 3}, origin_bits::AbstractArray{T, 3}, R::ScrambleMethod) where {T <: Integer}

In place version of a ScrambleMethod (`MatousekScramble` or `DigitalShift`) for the "bit" array.
This is faster to use this functions for multiple scramble of the same array (use `generate_design_matrices`).
"""
function randomize_bits!(random_bits::AbstractArray{T, 3},
    origin_bits::AbstractArray{T, 3},
    R::MatousekScramble) where {T <: Integer}
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    pad, n, d = size(origin_bits)
    b = R.base
    rng = R.rng
    m = logi(b, n)
    @assert m≥1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek(rng, m, b)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # matousek_M (m×m) * origin_bits (m×n) .+ matousek_C (m×1) 
        @views random_bits[1:m, :, s] .= (matousek_M * origin_bits[1:m, :, s] .+
                                          matousek_C) .% b
    end

    # Paste in random entries for bits after m'th one
    if pad > m
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
    matousek_M, matousek_C
end

getmatousek(m::Integer, b::Integer) = getmatousek(Random.GLOBAL_RNG, m, b)

"""
```julia
DigitalShift <: ScrambleMethod
```

Digital shift.
`randomize(x, R::DigitalShift)` returns a scrambled version of `x`.

The scramble method is Digital Shift.
It scrambles each coordinate in base `b` as `yₖ = (xₖ + Uₖ) mod b` where `Uₖ ∼ 𝕌({0:b-1})`.
`U` is the same for every point `points` but i.i.d. along every dimension.
"""
Base.@kwdef struct DigitalShift{I <: Integer} <: ScrambleMethod
    base::I
    pad::I = 32
    rng::AbstractRNG = Random.GLOBAL_RNG
end

function randomize_bits!(random_bits::AbstractArray{T, 3},
    origin_bits::AbstractArray{T, 3},
    R::DigitalShift) where {T <: Integer}
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    pad, n, d = size(origin_bits)
    b = R.base
    rng = R.rng
    m = logi(b, n)
    @assert m≥1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        DS = rand(rng, 0:(b - 1), m)

        # xₖ = (aₖ + Cₖ) mod b where xₖ is the k element in base b
        # origin_bits (m×n) .+ DS (m×1) 
        @views random_bits[1:m, :, s] .= (origin_bits[1:m, :, s] .+ DS) .% b
    end

    # Paste in random entries for bits after m'th one
    if pad > m
        # random_bits[(m + 1):pad, :, :] = rand(rng, 0:(b - 1), n * d * (pad - m))
        rand!(rng, @view(random_bits[(m + 1):pad, :, :]), 0:(b - 1))
    end
end
