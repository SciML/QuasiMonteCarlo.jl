"""
    unif2bits(x<:AbstractArray, b::Integer; M=32)
Return the b-adic decomposition of all element y of an array y = ∑ₖ yₖ/bᵏ a number yₖ∈[0,1[ -> [y₁, ⋯, yₘ] 
"""
function unif2bits(x::AbstractArray, b::Integer; M = 32)
    bits = zeros(Int, M, size(x)...)
    unif2bits!(bits, x, b)
    return bits
end

function unif2bits!(bits::AbstractArray, x::AbstractArray, b::Integer)
    @assert all([size(bits, d + 1) == size(x, d) for d in ndims(x)]) "Bit array of size $(size(bits)) instead of (M, $(size(x)))"
    for i in CartesianIndices(x)
        unif2bits!(@view(bits[:, i]), x[i], b)
    end
end

function unif2bits(x::AbstractArray; M = 32)
    bits = BitArray(undef, M, size(x)...)
    unif2bits!(bits, x, 2)
    return bits
end
"""
    unif2bits(y<:Real, b::Integer; M=32)
Return the b-adic decomposition y = ∑ₖ yₖ/bᵏ a number y∈[0,1[ -> [y₁, ⋯, yₘ] 
"""
function unif2bits(y::Real, b::Integer; M = 32)
    bits = zeros(Int, M)
    unif2bits!(bits, y, b)
    return bits
end
"""
    unif2bits(y; M=32)
Return the binary decomposition y = ∑ₖ yₖ/2ᵏ a number y∈[0,1[ -> [y₁, ⋯, yₘ] 
"""
function unif2bits(y; M = 32)
    bits = BitArray(undef, M)
    unif2bits!(bits, y, 2)
    return bits
end

check_zero(a::Rational; kwargs...) = a == 0
#! isequal(0.18518518518518515... -1/3^2-2/3^3 , 0) Should be true but is not! It fouls the binary expansion! Hence the isapprox with a hand tuned tolerence
check_zero(a::AbstractFloat; atol = 3e-16) = isapprox(a, 0, atol = atol)

function unif2bits!(bits::AbstractVector{<:Integer}, y, b::Integer; kwargs...)
    bits .= 0
    for j in eachindex(bits), bb in (b - 1):-1:1
        a = y - bb // b^j
        if check_zero(a; kwargs...)
            bits[j] = bb
            break # it breaks from the nested loop (see here)[https://stackoverflow.com/questions/39796234/how-to-break-out-of-nested-for-loops-in-julia]
        elseif a > 0
            bits[j] = bb
            y = a
        end
    end
end

#? Apparently this is not ideal to explicitly state the Type. 
#? See https://github.com/SciML/QuasiMonteCarlo.jl/issues/44#issuecomment-1328156825
#? Not sure how to do otherwise in this case though. 
"""
    bits2unif(::Type{T}, bits::AbstractVector{<:Integer},
                   b::Integer)
Convert a vector of M "bits" in base b into a number y∈[0,1[.
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
Convert a vector of M "bits" in base b into an integer.
"""
function bits2int(bit::AbstractVector, b::Integer)
    m = length(bit)
    y = 0
    for k in m:-1:1
        y = y * b + bit[k]
    end
    return y
end
