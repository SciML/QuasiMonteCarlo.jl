using QuasiMonteCarlo
using Aqua, Test
# Aqua.test_all(QuasiMonteCarlo)

using Compat
using Statistics, LinearAlgebra, StatsBase, Random
using Primes, Combinatorics, Distributions, InvertedIndices
using HypothesisTests

# struct InertSampler <: Random.AbstractRNG end
# InertSampler(args...; kwargs...) = InertSampler()
# Random.rand(::InertSampler, ::Type{T}) where {T} = zero(T)
# Random.rand(::InertSampler) = 0
# Random.shuffle!(::InertSampler, arg::AbstractArray) = arg

@views function embiggen!(a::Vector{T}, new_len::Integer, pad::T) where {T}
    old_len = length(a)
    @assert new_len≥old_len "Can't embiggen something smaller"
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

rng = MersenneTwister(1776)

#1D
lb = 0.0
ub = 1.0
n = 8
d = 1

for point_constructor in [
    FaureSample(),
    GridSample(),
    HaltonSample(),
    KroneckerSample(),
    LatinHypercubeSample(),
    LatticeRuleSample(),
    RandomSample(),
    SobolSample(),
]
    @show point_constructor
    s = QuasiMonteCarlo.sample(n, lb, ub, point_constructor)
    @test all(all(x .<= ub) for x in eachcol(s))
    @test all(all(x .>= lb) for x in eachcol(s))
    @test isa(s, Matrix) == true
    @test size(s) == (d, n)
end

for (sampler, d) in Iterators.product([Cauchy(), Normal(0, 4)], 1:3)
    @show sampler
    s = QuasiMonteCarlo.sample(n, d, sampler)
    @test s isa Matrix
    @test size(s) == (d, n)
end

@testset "1D" begin
    # @testset "SectionSample" begin
    #     constrained_val = 1.0
    #     point_constructor = SectionSample([NaN64], RandomSample())
    #     s = QuasiMonteCarlo.sample(n, lb, ub, point_constructor)
    #     @test s isa Vector{Float64}
    #     @test length(s) == n
    #     @test all(x -> lb ≤ x ≤ ub, s)
    #     @test !all(==(constrained_val), s)
    #     @test QuasiMonteCarlo.free_dimensions(point_constructor) == [1]
    #     @test QuasiMonteCarlo.fixed_dimensions(point_constructor) == Int[]

    #     point_constructor = SectionSample([constrained_val], RandomSample())
    #     s = QuasiMonteCarlo.sample(n, lb, ub, point_constructor)
    #     @test s isa Vector{Float64} && length(s) == n && all(x -> lb ≤ x ≤ ub, s)
    #     @test all(==(constrained_val), s)
    #     @test QuasiMonteCarlo.free_dimensions(point_constructor) == Int[]
    #     @test QuasiMonteCarlo.fixed_dimensions(point_constructor) == [1]

    #     #n = 0
    #     @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
    #                                                                          point_constructor)
    # end
end

#ND
@testset "GridSample" begin
    lb = [0.0, 0.0]
    ub = [1.0, 1.0]
    n = 16
    d = 4

    s = QuasiMonteCarlo.sample(n, lb, ub, GridSample())
    s = sortslices(s; dims = 2)
    differences = diff(s; dims = 2)
    @test all(≈(differences[1]), differences)
    μ = mean(s; dims = 2)
    variance = var(s; corrected = false, dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=2 / sqrt(n)
        @test variance[i]≈1 / 12 rtol=2 / sqrt(n)
    end
end

@testset "RandomSample" begin
    lb = [0.0, 0.0]
    ub = [1.0, 1.0]
    n = 20_000
    d = length(lb)

    s = QuasiMonteCarlo.sample(n, lb, ub, RandomSample(rng))
    @test size(s) == (d, n)

    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=2 / sqrt(n)
        @test variance[i]≈1 / 12 rtol=2 / sqrt(n)
    end
    @test pvalue(SignedRankTest(eachrow(s)...)) > 0.0001
end

@testset "LHS" begin
    lb = [0.0, 0.0]
    ub = [1.0, 1.0]
    d = length(lb)
    n = 20_000

    s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample(rng))
    @test size(s) == (d, n)

    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=2 / sqrt(n)
        @test variance[i]≈1 / 12 rtol=2 / sqrt(n)
    end
    @test pvalue(SignedRankTest(eachrow(s)...)) > 0.0001

    # A LHS is a scrambled `(λ=1, t=0, m=1, s=d)`-net in base `n`
    # See Cororollary 17.1 of [Monte Carlo theory, methods, and examples](https://artowen.su.domains/mc/qmcstuff.pdf).
    @test istmsnet(s, λ = 1, t = 0, m = 1, s = d, base = n)
end

@testset "Van der Corput Sequence" begin
    lb = 0
    ub = 1
    for base in [2, 3, 4]
        n = base^4
        s = QuasiMonteCarlo.sample(n, lb, ub, VanDerCorputSample(base, NoRand()))
        @test all(diff(sort(s)) .≈ 1 / n)
        @test all(s .≥ lb)
        @test all(s .≤ ub)
        @test 1 - maximum(s) ≈ minimum(s)
        @test minimum(s) ≈ inv(2n)
        @test mean(s) ≈ 0.5
        @test var(s; corrected = false)≈1 / 12 rtol=2 / sqrt(n)
    end
end

@testset "SobolSample" begin
    lb = [0.0, 0.0, 0.0]
    ub = [1.0, 1.0, 1.0]
    d = length(lb)
    power = 6
    base = 2
    n = base^power

    s = QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
    @test s isa Matrix
    @test size(s) == (d, n)
    vdc = QuasiMonteCarlo.sample(n, 1, VanDerCorputSample(base, NoRand()))
    sort!(vdc)
    for (i, j) in combinations(1:d, 2)
        @test pvalue(SignedRankTest(s[i, :], s[j, :])) > 0.0001
    end
    @test istmsnet(s; λ = 1, t = power - d, m = power, s = d, base)
    for dim in sort.(eachrow(s))
        # all dimensions should be base-2 van der Corput
        @test dim ≈ vdc
    end
    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=2 / n
        @test variance[i]≈1 / 12 rtol=2 / n
    end
end

@testset "Faure Sample" begin
    #FaureSample()
    d = 17
    base = nextprime(d)
    power = 2
    n = 17^2
    sampler = FaureSample()
    @test_throws ArgumentError QuasiMonteCarlo.sample(d + 1, d, sampler)
    @test_throws ArgumentError QuasiMonteCarlo.sample(d^2 + 1, d, sampler)
    s = sortslices(QuasiMonteCarlo.sample(n, d, sampler); dims = 2)
    # FaureSample() generates centered boxes, unlike DiceDesign
    r = sortslices(include("rfaure.jl")'; dims = 2) .+ inv(2base^(power + 1))

    @test s isa Matrix
    @test size(s) == (d, n)
    @test mean(abs, s - r) ≤ 2 / n
    @test s ≈ r

    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=2 / n
        @test variance[i]≈1 / 12 rtol=2 / n
    end
    for (i, j) in combinations(1:d, 2)
        @test pvalue(SignedRankTest(s[i, :], s[j, :])) > 0.0001
    end
    # test 5d stratification of first 3 primes
    power = 5
    for d in (2, 3, 5)
        base = nextprime(d)
        n = base^power
        s = QuasiMonteCarlo.sample(n, d, sampler)
        @test istmsnet(s; λ = 1, t = 0, m = power, s = d, base)
    end
    # test 2d stratification of next 3 primes
    power = 2
    for d in (7, 11, 13)
        base = nextprime(d)
        λ = 2
        n = λ * base^power
        s = QuasiMonteCarlo.sample(n, d, sampler)
        @test istmsnet(s; λ, t = 0, m = power, s = d, base)  # test 3d stratification of next 3 primes
    end
end

@testset "Halton Sequence" begin
    d = 4
    lb = zeros(d)
    ub = ones(d)
    bases = nextprimes(1, d)
    n = prod(bases)^2
    s = QuasiMonteCarlo.sample(n, lb, ub, HaltonSample())
    @test isa(s, Matrix)
    @test size(s) == (d, n)
    sorted = reduce(vcat, sort.(eachslice(s; dims = 1))')
    each_dim = eachrow(sorted)

    # each 1d sequence should have base b stratification property
    # (property inherited from van der Corput)
    @test all(zip(each_dim, bases)) do (seq, base)
        theoretical_count = n ÷ base
        part = Iterators.partition(seq, theoretical_count)
        all(enumerate(part)) do (i, subseq)
            all(subseq) do x
                i - 1 ≤ base * x ≤ i
            end
        end
    end
    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=1 / sqrt(n)
        @test variance[i]≈1 / 12 rtol=1 / sqrt(n)
    end
    for (i, j) in combinations(1:d, 2)
        @test pvalue(SignedRankTest(s[i, :], s[j, :])) > 0.0001
    end
end

@testset "LatticeRuleSample" begin
    #LatticeRuleSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, LatticeRuleSample())
    @test isa(s, Matrix)
    @test size(s) == (d, n)
    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i]≈0.5 atol=3 / n
        @test variance[i]≈1 / 12 rtol=3 / n
    end
end

@testset "Kronecker" begin
    d = 2
    n = 100
    ρ = 0.7548776662466927
    s = QuasiMonteCarlo.sample(n, d, KroneckerSample([ρ, ρ^2], NoRand()))
    t = QuasiMonteCarlo.sample(n, d, GoldenSample())
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test s ≈ t
    differences = eachcol(mod.(diff(s; dims = 2), 1))
    @test all(x -> x ≈ first(differences), differences)
end

# @testset "Section Sample" begin
#     lb = [0.0, 0.0]
#     ub = [1.0, 1.0]
#     n = 5
#     d = 2
#     constrained_val = 1.0
#     point_constructor = SectionSample([NaN64, constrained_val], RandomSample())
#     s = QuasiMonteCarlo.sample(n, lb, ub, point_constructor)
#     @test all(s[2, :] .== constrained_val)
#     @test isa(s, Matrix{Float64})
#     @test size(s) == (d, n)
#     @test all(lb[1] .≤ s[1, :] .≤ ub[1])
#     @test QuasiMonteCarlo.fixed_dimensions(point_constructor) == [2]
#     @test QuasiMonteCarlo.free_dimensions(point_constructor) == [1]

#     constrained_val = 2.0
#     d = 3
#     point_constructor = SectionSample([NaN64, constrained_val, NaN64], SobolSample())
#     @test QuasiMonteCarlo.free_dimensions(point_constructor) == [1, 3]
#     @test QuasiMonteCarlo.fixed_dimensions(point_constructor) == [2]
#     lb = [-1.0, constrained_val, -1.0]
#     ub = [6.0, constrained_val, 6.0]
#     s = QuasiMonteCarlo.sample(n, lb, ub, point_constructor)

#     @test all(lb[1] .≤ s[1, :] .≤ ub[1])
#     @test all(lb[3] .≤ s[3, :] .≤ ub[3])
#     @test all(s[2, :] .== constrained_val)
#     @test isa(s, Matrix{Float64})
#     @test size(s) == (d, n)

#     #n = 0
#     @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
#                                                                          SectionSample([
#                                                                                            NaN64,
#                                                                                            constrained_val,
#                                                                                        ],
#                                                                                        RandomSample()))
# end

@testset "generate_design_matrices" begin
    d = 4
    m = 8
    lb = zeros(d)
    ub = ones(d)
    num_mat = 3
    n = 2^m
    algorithms = [
        RandomSample(),
        LatinHypercubeSample(),
        SobolSample(R = OwenScramble(base = 2, pad = m)),
        SobolSample(),
        LatticeRuleSample(R = Shift()),
        SobolSample(R = MatousekScramble(base = 2, pad = m)),
    ]
    for algorithm in algorithms
        Ms = QuasiMonteCarlo.generate_design_matrices(n, lb, ub, algorithm, num_mat)
        @test length(Ms) == num_mat
        A = Ms[1]
        B = Ms[2]
        @test isa(A, Matrix{typeof(A[1][1])}) == true
        @test size(A) == (d, n)
        @test isa(B, Matrix{typeof(B[1][1])}) == true
        @test size(B) == (d, n)
    end
end

@testset "Randomized Quasi Monte Carlo: conversion functions" begin
    b = 2
    x = (0 // 16):(1 // 16):(15 // 16)
    bits = QuasiMonteCarlo.unif2bits(x, b)
    y = [QuasiMonteCarlo.bits2unif(s, b) for s in eachcol(bits)]
    @test x == y

    b = 3
    x = (0 // 27):(1 // 27):(26 // 27)
    bits = QuasiMonteCarlo.unif2bits(x, b, pad = 8)
    y = [QuasiMonteCarlo.bits2unif(Rational, s, b) for s in eachcol(bits)]
    @test x == y
end

@testset "Randomized Quasi Monte Carlo" begin
    m = 6
    d = 5
    b = QuasiMonteCarlo.nextprime(d)
    N = b^m # Number of points
    pad = 2m

    # Unrandomized low discrepency sequence
    u_faure = QuasiMonteCarlo.sample(N, d, FaureSample())

    # Randomized version
    @test u_faure == randomize(u_faure, NoRand())
    u_nus = randomize(u_faure, OwenScramble(base = b, pad = pad))
    u_lms = randomize(u_faure, MatousekScramble(base = b, pad = pad))
    u_digital_shift = randomize(u_faure, DigitalShift(base = b, pad = pad))
    u_shift = randomize(u_faure, Shift())

    u_nus = QuasiMonteCarlo.sample(N, d, FaureSample(OwenScramble(base = b, pad = pad)))
end

@testset "Randomized Quasi Monte Carlo Rational Scrambling" begin
    m = 5
    d = 2
    b = QuasiMonteCarlo.nextprime(d)
    N = b^m # Number of points
    pad = m

    # Unrandomized low discrepency sequence
    u_sobol = Rational.(QuasiMonteCarlo.sample(N, d, SobolSample()))
    # Randomized version
    u_nus = randomize(u_sobol, OwenScramble(base = b, pad = pad))
    u_lms = randomize(u_sobol, MatousekScramble(base = b, pad = pad))
    u_digital_shift = randomize(u_sobol, DigitalShift(base = b, pad = pad))
    @test eltype(u_nus) <: Rational
    @test eltype(u_lms) <: Rational
    @test eltype(u_digital_shift) <: Rational
end

@testset "Sobol' sequence are (tₛ,m,s)-net in base 2 even after scrambling" begin
    #! Note that this is the first 2ᵐ elements of the sequence which are tested
    #! Note that the Sobol' sequence from Sobol.jl is the "controversed" version where the first element have be removed

    # Table of qualitity parameter t with respect to dimension s for Sobol' sequence.
    # See "The algebraic-geometry approach to low-discrepancy sequences" H Niederreiter, C Xing - Monte Carlo and Quasi-Monte Carlo Methods 1996, 1998
    t_sobol = [0, 0, 1, 3, 5, 8, 11, 15, 19, 23, 27, 31, 35, 40, 45, 50, 55, 60, 65, 71]

    m = 12
    base = 2
    n = base^m
    λ = 1
    v = [
        NoRand(),
        OwenScramble(base = base, pad = m),
        MatousekScramble(base = base, pad = m),
        DigitalShift(base = base, pad = m),
    ]
    pass = Array{Bool}(undef, length(v), m)
    for (s, t) in enumerate(t_sobol[1:m])
        net = Rational.(QuasiMonteCarlo.sample(n, s, SobolSample()))
        for j in eachindex(v)
            pass[j, s] = istmsnet(randomize(net, v[j]); λ, t, m, s,
                                  base = base)
        end
    end
    @test all(pass)
end

@testset "Faure sequence are (0,m,s)-net even after scrambling" begin
    # Faure sequence are exactly (t=0, s)-sequence in base b=nextprime(s)
    m = 4
    pad = m
    λ = 1
    t = 0

    pass = Array{Bool}(undef, 4, m)
    for s in 1:m
        net = Rational{BigInt}.(QuasiMonteCarlo.sample(nextprime(s)^m, s, FaureSample())) # Convert the sequence in Rational{BigInt} (needed to scramble)
        pass[1, s] = istmsnet(randomize(net, NoRand()); λ, t, m, s,
                              base = nextprime(s))
        pass[2, s] = istmsnet(randomize(net, OwenScramble(base = nextprime(s), pad = m));
                              λ, t, m, s, base = nextprime(s))
        pass[3, s] = istmsnet(randomize(net,
                                        MatousekScramble(base = nextprime(s), pad = m));
                              λ, t, m, s, base = nextprime(s))
        pass[4, s] = istmsnet(randomize(net, DigitalShift(base = nextprime(s), pad = m));
                              λ, t, m, s, base = nextprime(s))
    end
    @test all(pass)
end
