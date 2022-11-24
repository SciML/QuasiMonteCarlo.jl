""" 
    nested_uniform_scramble(points::AbstractArray; M=32)
Return a scrambled version of the `points`. 
The scrambling method is Nested Uniform Scrambling which was introduced in Owen (1995).
`M` is the number of bits used for each points. One needs `M ≥ log(base, n)`. 
"""
function nested_uniform_scramble(rng::AbstractRNG, points::AbstractArray{T, N}, b::Integer;
                                 M = 32) where {T, N}
    random_points = similar(points)
    unrandomized_bits = unif2bits(points, b, M = M)
    indices = which_permutation(unrandomized_bits, b)
    random_bits = similar(unrandomized_bits)
    nested_uniform_scramble_bit!(rng, random_bits, unrandomized_bits, indices, b)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
    return random_points
end

function nested_uniform_scramble(points::AbstractArray, b::Integer; M = 32)
    nested_uniform_scramble(Random.default_rng(), points, b; M = M)
end

""" 
    nested_uniform_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{<:Integer, 3}, origin_bits::AbstractArray{<:Integer, 3}, indices::AbstractArray{T, 3} where {T <: Integer}, b::Integer)
In place version of Nested Uniform Scrambling (for the bit array). This is faster to use this functions for multiple scrambling of the same array.
"""
function nested_uniform_scramble_bit!(rng::AbstractRNG,
                                      random_bits::AbstractArray{<:Integer, 3},
                                      origin_bits::AbstractArray{<:Integer, 3},
                                      indices::AbstractArray{T, 3} where {T <: Integer},
                                      b::Integer)
    # in place nested uniform Scramble.
    #
    m, n, d = size(indices)
    M = size(random_bits, 3)
    @assert m≥1 "We need m ≥ 1" # m=0 causes awkward corner case below.

    for s in 1:d
        theperms = getpermset(rng, m, b)         # Permutations to apply to bits 1:m
        for k in 1:m                             # Here is where we want m > 0 so the loop works ok
            random_bits[k, :, s] = (origin_bits[k, :, s] + theperms[k, indices[k, :, s]]) .%
                                   b   # permutation by adding a bit modulo b
        end
    end
    if M > m     # Paste in random entries for bits after m'th one
        random_bits[(m + 1):M, :, :] = rand(rng, 0:(b - 1), n * d * (M - m))
    end
end

function getpermset(rng::AbstractRNG, m::Integer, b::Integer)
    # Get b^(k-1) random binary permutations for k=1 ... m
    # m will ordinarily be m when there are n=b^m points
    #
    y = zeros(Int, m, b^(m - 1))
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
It assigns for each points (in every dimensions) `m` number corresponding to its location on the slices 1, 1/b, 1/b², ..., 1/bᵐ⁻¹ of the axes [0,1[.
This also can be used to verify some equidistribution prorepreties.
Here we create the `indices` array `m`, and not `M`. Indeed `(t,m,d)-net` in base `b` are scrambled up to the `1/bᵐ` component. 
Higher order components are just used i.i.d `Uₖ ∼ 𝐔({0:b-1})` in `nested_uniform_scramble_bit!`.
"""
function which_permutation(bits::AbstractArray, b)
    n, d = size(bits)[2:end]
    m = logi(b, n)

    indices = zeros(Int, m, n, d)
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

""" 
    nested_uniform_scramble(points::AbstractArray, b::Integer; M=32)
Return a scrambled version of the `points`. 
The scrambling method is Linear Matrix Scrambling which was introduced in Matousek (1998).
`M` is the number of bits used for each points. One need `M ≥ log(base, n)`. 
"""
function linear_matrix_scramble(rng::AbstractRNG, points::AbstractArray{T, N}, b::Integer;
                                M = 32) where {T, N}
    unrandomized_bits = unif2bits(points, b, M = M)
    random_bits = similar(unrandomized_bits)
    linear_matrix_scramble_bit!(rng, random_bits, unrandomized_bits, b)
    random_points = similar(points)
    for i in CartesianIndices(random_points)
        random_points[i] = bits2unif(T, @view(random_bits[:, i]), b)
    end
    return random_points
end

function linear_matrix_scramble(points::AbstractArray, b::Integer; M = 32)
    linear_matrix_scramble(Random.default_rng(), points, b; M = M)
end

#? Weird it should be faster than nested uniform Scrambling but here it is not at all.-> look for other implementation and paper
""" 
    linear_matrix_scramble_bit!(rng::AbstractRNG, random_bits::AbstractArray{<:Integer, 3}, origin_bits::AbstractArray{<:Integer, 3}, b::Integer)
In place version of Linear Matrix Scrambling (for the bit array). This is faster to use this functions for multiple scrambling of the same array.
"""
function linear_matrix_scramble_bit!(rng::AbstractRNG,
                                     random_bits::AbstractArray{<:Integer, 3},
                                     origin_bits::AbstractArray{<:Integer, 3}, b::Integer)
    # https://statweb.stanford.edu/~owen/mc/ Chapter 17.6 around equation (17.15).
    #
    M, n, d = size(origin_bits)
    m = logi(b, n)
    @assert m≥1 "We need m ≥ 1" # m=0 causes awkward corner case below.  Caller handles that case specially.

    for s in 1:d
        # Permutations matrix and shift to apply to bits 1:m
        matousek_M, matousek_C = getmatousek(rng, m, b)

        # xₖ = (∑ₗ Mₖₗ aₗ + Cₖ) mod b where xₖ is the k element in base b
        # matousek_M (m×m) * origin_bits (m×n) .+ matousek_C (m×1) 
        random_bits[1:m, :, s] = (matousek_M * origin_bits[1:m, :, s] .+ matousek_C) .% b
    end

    # Paste in random entries for bits after m'th one
    if M > m
        random_bits[(m + 1):M, :, :] = rand(rng, Bool, n * d * (M - m))
    end
end

""" 
    getmatousek(rng::AbstractRNG, m::Integer, b::Integer)
Genereate the Matousek linear scramble in base b for one of the d components
It produces a m x m bit matrix matousek_M and a length m bit vector matousek_C
"""
function getmatousek(rng::AbstractRNG, m::Integer, b::Integer)
    matousek_M = LowerTriangular(zeros(Int, m, m)) + Diagonal(rand(rng, 1:(b - 1), m)) # Mₖₖ ∼ U{1, ⋯, b-1}
    matousek_C = rand(rng, 0:(b - 1), m)
    for i in 2:m
        for j in 1:(i - 1)
            matousek_M[i, j] = rand(rng, 0:(b - 1))
        end
    end
    matousek_M, matousek_C
end

getmatousek(m::Integer, b::Integer) = getmatousek(Random.GLOBAL_RNG, m, b)