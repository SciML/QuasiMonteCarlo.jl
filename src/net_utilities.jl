"""
    unif2bits(x<:AbstractArray, b::Integer; pad=32)

Return the b-ary decomposition of all element y of an array y = ∑ₖ yₖ/bᵏ a number yₖ∈[0,1[ -> [y₁, ⋯, yₘ]
"""
function unif2bits(x::AbstractArray, b::I; pad = 32) where I<:Integer
    bits = zeros(I, pad, size(x)...)
    unif2bits!(bits, x, b)
    return bits
end

function unif2bits!(bits::AbstractArray, x::AbstractArray, b::Integer)
    @assert all([size(bits, d + 1) == size(x, d) for d in ndims(x)]) "Bit array of size $(size(bits)) instead of (pad, $(size(x)))"
    for i in CartesianIndices(x)
        unif2bits!(@view(bits[:, i]), x[i], b)
    end
end

"""
    unif2bits(y<:Real, b::Integer; pad=32)

Return the b-ary decomposition y = ∑ₖ yₖ/bᵏ a number y∈[0,1[ -> [y₁, ⋯, yₘ]
"""
function unif2bits(y::Real, b::I; pad = 32) where I<:Integer
    bits = zeros(I, pad)
    unif2bits!(bits, y, b)
    return bits
end

# Inspired by digits!(a::AbstractVector{T}, n::Integer; base) where T<:Integer in Base at intfuncs.jl:926
function unif2bits!(bits::AbstractVector{<:Integer}, y, b::Integer; kwargs...)
    invbase = inv(b)
    for i in eachindex(bits)
        r, y = divrem(y, invbase^i)
        bits[i] = r
    end
end

#? Apparently this is not ideal to explicitly state the Type. 
#? See https://github.com/SciML/QuasiMonteCarlo.jl/issues/44#issuecomment-1328156825
#? Not sure how to do otherwise in this case though. 
"""
    bits2unif(::Type{T}, bits::AbstractVector{<:Integer},
                   b::Integer)

Convert a vector of pad "bits" in base b into a number y∈[0,1[.
"""
function bits2unif(::Type{T}, bits::AbstractVector{<:Integer},
    b::Integer) where {T <: Rational}
    # Turn sequence of bits into a point in [0,1)
    # First bits are highest order
    y = zero(T)
    for j in lastindex(bits):-1:1
        y = (y + bits[j]) // b
    end
    return y
end

function bits2unif(::Type{T}, bits::AbstractVector{<:Integer},
    b::Integer) where {T <: AbstractFloat}
    # Turn sequence of bits into a point in [0,1)
    # First bits are highest order
    y = zero(T)
    for j in lastindex(bits):-1:1
        y = (y + bits[j]) / b
    end
    return y
end

function bits2unif(bits::AbstractVector{<:Integer}, b::Integer)
    bits2unif(Float64, bits::AbstractVector{<:Integer}, b::Integer)
end

#? This seems faster than @evalpoly(b, $bit...)
#?  bi = rand(0:2,32);
#? @btime @evalpoly(3, $bi...)
#?   500.515 ns (1 allocation: 272 bytes)
#? @btime QuasiMonteCarlo.bits2int($bi, 3)
#?   13.113 ns (0 allocations: 0 bytes)  
"""
    bits2int(bit::AbstractMatrix{<:Integer}, b::Integer)

Convert a vector of pad "bits" in base b into an integer.
"""
function bits2int(bit::AbstractVector, b::I) where I<:Integer
    m = length(bit)
    y = zero(I)
    for k in m:-1:1
        y = y * b + bit[k]
    end
    return y
end
