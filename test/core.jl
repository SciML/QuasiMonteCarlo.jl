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
    @assert size(net, 2) == λ * (base^m) "Number of points must be as size(net, 2) = $(size(net, 2)) == λ * (base^m) = $(λ * (base^m))"
    @assert size(net, 1) == s "Dimension must be as size(net, 2) = $(size(net, 2)) == s = $s"
    @assert all(0 .≤ net .< 1) "All points must be in [0,1)"

    # The elementary boxes are products of the 1D pieces of `mince(interval(0, 1), divisions)`, which
    # tile [0,1) as [left end, next left end), so a coordinate lies in the last piece whose left end is
    # ≤ it. Binning each point once costs O(n s log n) per direction instead of O(n × #boxes).
    for stepsize in multiexponents(s, m - t)
        divisions = base .^ stepsize
        left_ends = [inf.(mince(interval(zero(T), one(T)), k)) for k in divisions]
        counts = zeros(Int, prod(divisions))
        for point in eachcol(net)
            box, stride = 1, 1  # column-major linear index of the box holding `point`
            for (ends, x) in zip(left_ends, point)
                box += (searchsortedlast(ends, x; lt = <) - 1) * stride
                stride *= length(ends)
            end
            counts[box] += 1
        end
        if !all(==(λ * base^t), counts)
            println("Errors in direction k = $stepsize")
            return false
        end
    end
    return true
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
    # A full grid holds every cell center once, in an order the rng sets.
    for (cells, d) in ((10, 1), (8, 2), (3, 3))
        centers = Float64.((2(0:(cells - 1)) .+ 1) .// 2cells)
        points = QuasiMonteCarlo.sample(cells^d, d, GridSample(Xoshiro(1776)))
        @test size(points) == (d, cells^d)
        @test Set(collect.(eachcol(points))) == Set(collect.(Iterators.product(fill(centers, d)...)))
        @test points == QuasiMonteCarlo.sample(cells^d, d, GridSample(Xoshiro(1776)))
        @test points != QuasiMonteCarlo.sample(cells^d, d, GridSample(Xoshiro(1777)))
    end
    # Other sizes take distinct points of the smallest grid that holds them.
    for (n, d, cells) in ((5, 2, 3), (10, 3, 3), (3, 70, 2), (1, 4, 1))
        points = QuasiMonteCarlo.sample(n, d, GridSample(Xoshiro(1)))
        @test size(points) == (d, n)
        @test allunique(eachcol(points))
        @test all(in(Float64.((2(0:(cells - 1)) .+ 1) .// 2cells)), points)
    end
    @test_throws ArgumentError QuasiMonteCarlo.sample(4, 0, GridSample())
    # The float root misses the exact count at these sizes.
    @test QuasiMonteCarlo.grid_cells(Int64(10)^18, 3) === Int64(10)^6
    @test QuasiMonteCarlo.grid_cells(Int64(10)^18 + 1, 3) === Int64(10)^6 + 1
    @test QuasiMonteCarlo.grid_cells(Int64(2)^62 + 1, 2) === Int64(2)^31 + 1
    @test QuasiMonteCarlo.grid_cells(typemax(Int64), 1) === typemax(Int64)
    # Scaling a full grid keeps its exact mean and variance.
    lb, ub = [-1.0, 10.0], [3.0, 12.0]
    s = QuasiMonteCarlo.sample(64, lb, ub, GridSample(Xoshiro(1776)))
    @test vec(mean(s; dims = 2)) ≈ (lb .+ ub) ./ 2
    @test vec(var(s; corrected = false, dims = 2)) ≈ (ub .- lb) .^ 2 .* (1 - 1 / 64) ./ 12
    for T in (Float64, Float32, Float16)
        points = QuasiMonteCarlo.sample(125, 5, GridSample(), T)
        @test eltype(points) == T && all(0 .< points .< 1)
    end
end

@testset "distinct_indices draws uniformly without replacement" begin
    # Each position of the order is uniform, so every leading run of grid points is unbiased.
    # (10, 10) shuffles all indices; (27, 10) rejects repeated draws.
    for (grid_size, n) in ((10, 10), (27, 10))
        rng = Xoshiro(1776)
        draws = [QuasiMonteCarlo.distinct_indices(rng, grid_size, n) for _ in 1:20_000]
        @test all(allunique, draws) && all(draw -> all(in(0:(grid_size - 1)), draw), draws)
        for position in (1, n)
            tallies = [count(draw -> draw[position] == index, draws) for index in 0:(grid_size - 1)]
            @test pvalue(ChisqTest(tallies)) > 1.0e-4
        end
    end
    @test QuasiMonteCarlo.distinct_indices(Xoshiro(1), 5, 0) == Int[]
    @test eltype(QuasiMonteCarlo.distinct_indices(Xoshiro(1), big(2)^70, 3)) == BigInt
    @test_throws ArgumentError QuasiMonteCarlo.distinct_indices(Xoshiro(1), 3, 4)
    @test_throws ArgumentError QuasiMonteCarlo.distinct_indices(Xoshiro(1), 3, -1)
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
        @test μ[i] ≈ 0.5 atol = 2 / sqrt(n)
        @test variance[i] ≈ 1 / 12 rtol = 2 / sqrt(n)
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
        @test μ[i] ≈ 0.5 atol = 2 / sqrt(n)
        @test variance[i] ≈ 1 / 12 rtol = 2 / sqrt(n)
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
        @test var(s; corrected = false) ≈ 1 / 12 rtol = 2 / sqrt(n)
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
        @test μ[i] ≈ 0.5 atol = 2 / n
        @test variance[i] ≈ 1 / 12 rtol = 2 / n
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
    # FaureSample() generates centered boxes, unlike DiceDesign: n = base^power cells of width 1 / n.
    r = sortslices(include("shared/rfaure.jl")'; dims = 2) .+ inv(2base^power)

    @test s isa Matrix
    @test size(s) == (d, n)
    @test mean(abs, s - r) ≤ 2 / n
    @test s ≈ r

    μ = mean(s; dims = 2)
    variance = var(s; dims = 2)
    for i in eachindex(μ)
        @test μ[i] ≈ 0.5 atol = 2 / n
        @test variance[i] ≈ 1 / 12 rtol = 2 / n
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

@testset "Faure nets put each point at the center of its cell" begin
    # Float `log(5, 125)` exceeds 3, but the offset uses exact `nextpow`, so 125 points stay centered.
    for (d, power) in ((1, 3), (3, 3), (5, 2), (5, 3), (7, 2))
        base = nextprime(d)
        n = base^power
        points = QuasiMonteCarlo.sample(n, d, FaureSample())
        # Each coordinate of a full net is a permutation of the n cell centers.
        @test all(sort(row) ≈ ((0:(n - 1)) .+ 0.5) ./ n for row in eachrow(points))
    end
    # With λ < base copies of base^power, the cells have one more digit.
    points = QuasiMonteCarlo.sample(2 * 3^2, 3, FaureSample())
    @test points .* 3^3 .- 0.5 ≈ round.(points .* 3^3 .- 0.5)
end

@testset "(Randomized) Halton Sequence" begin
    d = 4
    lb = zeros(d)
    ub = ones(d)
    bases = nextprimes(1, d)
    n = prod(bases)^2
    for sample_type in [HaltonSample(), RandomizedHaltonSample()]
        s = QuasiMonteCarlo.sample(n, lb, ub, sample_type)
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
            @test μ[i] ≈ 0.5 atol = 1 / sqrt(n)
            @test variance[i] ≈ 1 / 12 rtol = 1 / sqrt(n)
        end
        for (i, j) in combinations(1:d, 2)
            @test pvalue(SignedRankTest(s[i, :], s[j, :])) > 0.0001
        end
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
        @test μ[i] ≈ 0.5 atol = 3 / n
        @test variance[i] ≈ 1 / 12 rtol = 3 / n
    end
end

@testset "deterministic samplers never return 0 or 1" begin
    samplers = [FaureSample(), HaltonSample(), KroneckerSample(), LatticeRuleSample(), SobolSample()]
    for sampler in samplers, T in (Float64, Float32), (count, dims) in ((4, 1), (16, 2), (125, 5))
        points = QuasiMonteCarlo.sample(count, dims, sampler, T)
        @test eltype(points) == T
        @test all(0 .< points .< 1)
    end
    # The lattice moves each point up by half a last digit, so its first point is not the origin.
    lattice = QuasiMonteCarlo._lattice_rule(3)
    @test QuasiMonteCarlo.sample(64, 3, LatticeRuleSample()) == reduce(hcat, lattice[0:63]) .+ 2.0^-53
    @test QuasiMonteCarlo.sample(64, 3, LatticeRuleSample(), Float32)[:, 1] == fill(2.0f0^-24, 3)
    @test QuasiMonteCarlo.cell_midpoint(Float32, prevfloat(1.0)) == prevfloat(1.0f0)
    @test QuasiMonteCarlo.cell_midpoint(Rational{Int128}, 0.0) == 1 // Int128(2)^53
    # A float wider than the exponent range of Float64 scales in its own type.
    setprecision(BigFloat, 2048) do
        points = QuasiMonteCarlo.sample(4, 1, LatticeRuleSample(), BigFloat)
        @test points == [big(2)^-2048 0.5 + big(2)^-2048 0.25 + big(2)^-2048 0.75 + big(2)^-2048]
    end
    # Cells finer than the last digit of a float still give centers inside (0, 1).
    @test QuasiMonteCarlo.grid_center(Float16, 2048, 2049) == prevfloat(Float16(1))
    @test QuasiMonteCarlo.grid_center(Float16, 0, 2^30) == nextfloat(Float16(0))
    # `2n` overflows a 32-bit integer here, so the center needs widened arithmetic.
    @test QuasiMonteCarlo.grid_center(Float64, Int32(0), Int32(2)^30) == 2.0^-31
    @test QuasiMonteCarlo.grid_center(Float64, typemax(Int) - 1, typemax(Int)) ≈ 1 - 0.5 / typemax(Int)
    @test QuasiMonteCarlo.grid_center(Rational{Int}, 1, 4) == 3 // 8
    @test all(0 .< QuasiMonteCarlo.sample(2049, 1, GridSample(), Float16) .< 1)
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

@testset "iterators" begin
    d = 4
    m = 14
    n = 2^m
    num_mat = 4

    f(x) = prod(x) * 2^length(x) # test function ∫f(x)dᵈx = 1

    algorithms = [
        RandomSample(),
        Beta(0.2, 0.8),
        # LatticeRuleSample(R = Shift()), # TODO add support for LatticeRule
        SobolSample(R = OwenScramble(base = 2, pad = m)),
        SobolSample(R = MatousekScramble(base = 2, pad = m)),
        SobolSample(R = DigitalShift(base = 2, pad = m)),
    ]
    for algorithm in algorithms
        @show algorithm
        iterator = DesignMatrix(n, d, algorithm, num_mat)
        X = QuasiMonteCarlo.next!(iterator)
        @test isa(X, Matrix{eltype(X)}) == true
        μ = [mean(f(c) for c in eachcol(X)) for X in iterator] # Check that iterator do work!
        if !isa(algorithm, Beta) # The Beta distribution is not Uniform hence we do not expect the μ to be close to 1. We could compute the expected results but it would involve a β function dependency.
            @test μ ≈ ones(num_mat) atol = 5.0e-2 # the results for different randomization should all be close to 1. Arbitrarily we allow 5% error.
        end
    end
end

@testset "Types of output and intermediate arrays" begin
    # Scrambling methods use intermediate array for the bits scrambling. It is in general unnecessary to have this large array stored as Int64.
    #TODO test other randomization methods (here just scrambling) and QMC sequence (here just Sobol)
    #TODO here we test output and constructions of intermediate array, we could test that all operation are type stables
    d = 4
    m = 12
    n = 2^m

    for type in [Float32, Float64]
        for Int_type in [Int32, Int64]
            b, pad = Int_type(2), Int_type(32)
            scramblings = [
                OwenScramble(base = b, pad = pad)
                MatousekScramble(base = b, pad = pad)
                DigitalShift(base = b, pad = pad)
            ]
            for scrambling in scramblings
                output,
                    other_arrays... = QuasiMonteCarlo.initialize(
                    n,
                    d,
                    SobolSample(),
                    scrambling,
                    type
                )
                @test output isa Matrix{type}
                @test all(
                    [
                        other_array isa Array{Int_type, 3}
                            for other_array in other_arrays
                    ]
                )

                iterator = DesignMatrix(n, d, SobolSample(), scrambling, 1, type)
                scrambled_output = QuasiMonteCarlo.next!(iterator)
                @test scrambled_output isa Matrix{type}
            end
        end
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

    # Unrandomized low discrepancy sequence
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

    # Unrandomized low discrepancy sequence
    u_sobol = Rational.(QuasiMonteCarlo.sample(N, d, SobolSample()))
    # Randomized version
    u_nus = randomize(u_sobol, OwenScramble(base = b, pad = pad))
    u_lms = randomize(u_sobol, MatousekScramble(base = b, pad = pad))
    u_digital_shift = randomize(u_sobol, DigitalShift(base = b, pad = pad))
    @test eltype(u_nus) <: Rational
    @test eltype(u_lms) <: Rational
    @test eltype(u_digital_shift) <: Rational
end

@testset "istmsnet rejects point sets that are not nets" begin
    # Stratified along the first coordinate only: of the 2×2 boxes, the lower two hold 2 points each.
    x = (0:3) .// 4
    stratified_in_x = permutedims([x zero(x)])
    @test istmsnet(stratified_in_x[1:1, :]; λ = 1, t = 0, m = 2, s = 1, base = 2)
    @test !istmsnet(stratified_in_x; λ = 1, t = 0, m = 2, s = 2, base = 2)
    @test !istmsnet(stratified_in_x; λ = 1, t = 1, m = 2, s = 2, base = 2)
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
            pass[j, s] = istmsnet(
                randomize(net, v[j]); λ, t, m, s,
                base = base
            )
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
        pass[1, s] = istmsnet(
            randomize(net, NoRand()); λ, t, m, s,
            base = nextprime(s)
        )
        pass[2, s] = istmsnet(
            randomize(net, OwenScramble(base = nextprime(s), pad = m));
            λ, t, m, s, base = nextprime(s)
        )
        pass[3, s] = istmsnet(
            randomize(
                net,
                MatousekScramble(base = nextprime(s), pad = m)
            );
            λ, t, m, s, base = nextprime(s)
        )
        pass[4, s] = istmsnet(
            randomize(net, DigitalShift(base = nextprime(s), pad = m));
            λ, t, m, s, base = nextprime(s)
        )
    end
    @test all(pass)
end

"""The 32 leading digits of each coordinate of `points` scrambled by `Scramble` seeded by `seed`."""
function scrambled_words(Scramble, points, seed)
    scrambled = randomize(points, Scramble(base = 2, pad = 32, rng = Xoshiro(seed)))
    return UInt32.(floor.(scrambled .* 2.0^32))
end

"""The sorted numbers of points in each box of the grid that splits dimension `s` into `2^splits[s]` cells."""
function box_counts(words, splits)
    counts = Dict{Vector{UInt32}, Int}()
    for point in eachcol(words)
        box = UInt32[word >> (32 - split) for (word, split) in zip(point, splits)]
        counts[box] = get(counts, box, 0) + 1
    end
    return sort!(collect(values(counts)))
end

@testset "$Scramble scrambles a digital net into a shifted digital net of all pad digits" for Scramble in (MatousekScramble, DigitalShift)
    d, m = 4, 6
    n = 2^m
    sobol = QuasiMonteCarlo.sample(n, d, SobolSample())
    plain = UInt32.(floor.(sobol .* 2.0^32))
    words = scrambled_words(Scramble, sobol, 3)
    @test words != plain

    # Affine over GF(2): a ⊕ b ⊕ (any point) is again a point, so the points of each
    # dimension are a shifted subspace of dimension m.
    for dimension in axes(words, 1)
        points = Set(words[dimension, :])
        anchor = words[dimension, begin]
        @test length(points) == n
        @test all((a ⊻ b ⊻ anchor) in points for a in points for b in points)
    end

    # A longer sample keeps the points of a shorter one, and a wider one the dimensions of a
    # narrower one: the scramble does not depend on the number of points or dimensions.
    longer = QuasiMonteCarlo.sample(2n, d, SobolSample())
    @test scrambled_words(Scramble, longer, 3)[:, 1:37] == scrambled_words(Scramble, longer[:, 1:37], 3)
    @test scrambled_words(Scramble, sobol, 3)[1:2, :] == scrambled_words(Scramble, sobol[1:2, :], 3)

    # The leading digits of a scrambled point depend on the leading digits alone.
    cells = rand(Xoshiro(1), 0:(2^10 - 1), 1, 200)
    coarse = (cells .+ rand(Xoshiro(2), 0:(2^22 - 1), 1, 200) ./ 2^22) ./ 2^10
    fine = (cells .+ rand(Xoshiro(3), 0:(2^22 - 1), 1, 200) ./ 2^22) ./ 2^10
    @test floor.(Int, randomize(coarse, Scramble(base = 2, pad = 32, rng = Xoshiro(4))) .* 2^10) ==
        floor.(Int, randomize(fine, Scramble(base = 2, pad = 32, rng = Xoshiro(4))) .* 2^10)
    @test randomize(coarse, Scramble(base = 2, pad = 32, rng = Xoshiro(4))) != randomize(fine, Scramble(base = 2, pad = 32, rng = Xoshiro(4)))

    # Box counts are those of the unscrambled net: one point per elementary interval in the
    # first two dimensions, and the same counts as the unscrambled net in all four.
    for first_digits in 0:m
        @test box_counts(words[1:2, :], (first_digits, m - first_digits)) == ones(Int, n)
    end
    for splits in Iterators.product(fill(0:m, d)...)
        sum(splits) == m || continue
        @test box_counts(words, splits) == box_counts(plain, splits)
    end

    # Each scrambled point is uniform: its mean over 4000 seeds is within 4 standard errors of ½.
    draws = [scrambled_words(Scramble, fill(0.3, 1, 1), seed)[1, 1] / 2.0^32 for seed in 1:4000]
    @test abs(mean(draws) - 0.5) < 4 * sqrt(1 / 12 / length(draws))
end

@testset "MatousekScramble and DigitalShift do not depend on the thread" begin
    for Scramble in (MatousekScramble, DigitalShift)
        sobol = QuasiMonteCarlo.sample(64, 4, SobolSample())
        serial = [scrambled_words(Scramble, sobol, seed) for seed in 1:8]
        tasks = [Threads.@spawn(scrambled_words(Scramble, sobol, seed)) for seed in 1:8]
        @test fetch.(tasks) == serial
    end
end
