using QuasiMonteCarlo
using Test

using Statistics, LinearAlgebra, StatsBase, Random
using Primes, Combinatorics, Distributions, IntervalArithmetic
using HypothesisTests

# struct InertSampler <: Random.AbstractRNG end
# InertSampler(args...; kwargs...) = InertSampler()
# Random.rand(::InertSampler, ::Type{T}) where {T} = zero(T)
# Random.rand(::InertSampler) = 0
# Random.shuffle!(::InertSampler, arg::AbstractArray) = arg

@views function embiggen!(a::Vector{T}, new_len::Integer, pad::T) where {T}
    old_len = length(a)
    @assert new_len ≥ old_len "Can't embiggen something smaller"
    if new_len == old_len
        return a
    elseif new_len == old_len + 1
        push!(a, pad)
        return a
    else
        append!(a, fill(pad, new_len - old_len))
        return a
    end
end

####################
### T, pad, S NETS ###
####################

"""
    istmsnet(net::AbstractMatrix{T}; λ::I, t::I, m::I, d::I, base::I) where {I <: Integer, T <: Real}

Test if a point set `net` (`dim×n`) is a `(λ,t,m,s)`-net in base `b`.

`(λ,t,m,s)`-nets have strict equidistribution properties making them good QMC sequences.
Their definition and properties can be found in the book [Monte Carlo theory, methods, and examples](https://artowen.su.domains/mc/qmcstuff.pdf) by Art B. Owen.
See Definition 15.7 and for properties see Chapter 15 to 17.

The test is exact if the element of `net` are of type `Rational`. Indeed, in that case, one can exactly deal with points at the edge of intervals of the form [a,b)ᵈ.
The conversion `Float` to `Rational` is usually possible with usual nets, e.g., Sobol, Faure (may require `Rational{BigInt}.(net)`).
"""
function istmsnet(
        net::AbstractMatrix{T}; λ::I, t::I, m::I, s::I,
        base::I
    ) where {I <: Integer, T <: Real}
    pass = true

    @assert size(net, 2) == λ * (base^m) "Number of points must be as size(net, 2) = $(size(net, 2)) == λ * (base^m) = $(λ * (base^m))"
    @assert size(net, 1) == s "Dimension must be as size(net, 2) = $(size(net, 2)) == s = $s"
    @assert all(0 .≤ net .< 1) "All points must be in [0,1)"

    perms = multiexponents(s, m - t)
    for stepsize in perms
        intervals = mince(
            IntervalBox([interval(zero(T), one(T)) for i in 1:s]),
            NTuple{s, Int}(base .^ stepsize)
        )
        pass &= all(intervals) do intvl
            λ * base^t == count(point -> inCloseOpen(point, intvl), collect(eachcol((net))))
        end
        if !pass
            println("Errors in direction k = $stepsize")
            return pass
        end
    end
    return pass
end

"""
    in_halfopen(x, a)

Checks if the number `x` is a member of the interval `a` (close on the left and open on the right), treated as a set.
"""
function inCloseOpen(x::T, a::Interval) where {T <: Real}
    isinf(x) && return false
    return a.lo <= x < a.hi
end
inCloseOpen(X::AbstractVector, Y::IntervalBox{N, T}) where {N, T} = all(inCloseOpen.(X, Y))

rng = MersenneTwister(1776)

#1D
lb = 0.0
ub = 1.0
n = 8
d = 1
