using QuasiMonteCarlo, Distributions, StatsBase, Random
using Test

# For testing randomized QMC sequences by using the deterministic version

struct InertSampler <: Random.AbstractRNG end
InertSampler(args...; kwargs...) = InertSampler()
Random.rand(::InertSampler, ::Type{T}) where T = zero(T)
Random.shuffle!(::InertSampler, arg::AbstractArray) = arg


#1D
lb = 0.0
ub = 5.0
n = 5
d = 1
QuasiMonteCarlo.sample(n, lb, ub, GridSample(0.1))
QuasiMonteCarlo.sample(n, lb, ub, UniformSample())
QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())
QuasiMonteCarlo.sample(n, lb, ub, LatticeRuleSample())
QuasiMonteCarlo.sample(5, d, Cauchy())
QuasiMonteCarlo.sample(5, d, Normal(0, 4))

@testset "1D" begin
    @testset "LowDiscrepancySample" begin
        s = QuasiMonteCarlo.sample(n, 0.0, 1.0, LowDiscrepancySample(2, false))
        @test isa(s, Vector{Float64})
        @test size(s) == (n,)
        @test s ≈ [0.5, 0.25, 0.75, 0.125, 0.625]

        s = QuasiMonteCarlo.sample(n, 0, 1, LowDiscrepancySample(2, false))
        @test isa(s, Vector{Float64})
        @test size(s) == (n,)
        @test s ≈ [0.5, 0.25, 0.75, 0.125, 0.625]

        s = QuasiMonteCarlo.sample(n, zero(Float32), one(Float32),
                                   LowDiscrepancySample(2, false))
        @test isa(s, Vector{Float32})
        @test size(s) == (n,)
        @test s≈[0.5, 0.25, 0.75, 0.125, 0.625] rtol=1e-7
    end

    @testset "KroneckerSample" begin
        s = QuasiMonteCarlo.sample(n, lb, ub, KroneckerSample(sqrt(2), 0))
        @test s isa Vector{Float64} && length(s) == n && all(x -> lb ≤ x ≤ ub, s)
    end

    @testset "GoldenSample" begin
        s = QuasiMonteCarlo.sample(n, lb, ub, GoldenSample())
        @test s isa Vector{Float64} && length(s) == n && all(x -> lb ≤ x ≤ ub, s)
    end

    @testset "SectionSample" begin
        constrained_val = 1.0
        sampler = SectionSample([NaN64], UniformSample())
        s = QuasiMonteCarlo.sample(n, lb, ub, sampler)
        @test s isa Vector{Float64} && length(s) == n && all(x -> lb ≤ x ≤ ub, s)
        @test !all(==(constrained_val), s)
        @test QuasiMonteCarlo.free_dimensions(sampler) == [1]
        @test QuasiMonteCarlo.fixed_dimensions(sampler) == Int[]

        sampler = SectionSample([constrained_val], UniformSample())
        s = QuasiMonteCarlo.sample(n, lb, ub, sampler)
        @test s isa Vector{Float64} && length(s) == n && all(x -> lb ≤ x ≤ ub, s)
        @test all(==(constrained_val), s)
        @test QuasiMonteCarlo.free_dimensions(sampler) == Int[]
        @test QuasiMonteCarlo.fixed_dimensions(sampler) == [1]
    end
end

#ND
lb = [0.0, 0.0]
ub = [1.0, 1.0]
n = 5
d = 2

@testset "GridSample" begin
    #GridSample{T}
    s = QuasiMonteCarlo.sample(n, lb, ub, GridSample([0.1, 0.5]))
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "UniformSample" begin
    #UniformSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, UniformSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "SobolSample" begin
    #SobolSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "LHS" begin
    #LHS
    s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample(threading = true))
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "Faure Sample" begin
    #FaureSample()
    d = 17
    n = 17^2
    rng = InertSampler()
    @test_throws ArgumentError QuasiMonteCarlo.sample(d+1, d, FaureSample())
    @test_throws ArgumentError QuasiMonteCarlo.sample(d^2+1, d, FaureSample())
    s = sortslices(QuasiMonteCarlo.sample(n, d, FaureSample(rng)); dims=2)
    s == sortslices(QuasiMonteCarlo.sample(n, zeros(d), ones(d), FaureSample(rng)); dims=2)
    r = sortslices(include("rfaure.jl")'; dims=2)
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test mean(abs2, s - r) ≤ eps(Float32)
    @test s ≈ r

    # check RQMC stratification properties
    s = QuasiMonteCarlo.sample(n, d, FaureSample(MersenneTwister(0)))
    fallsin(x, args...) = all(@. (args-1) < x ≤ args)
    @test all(1:d) do dim_idx  # for every dimension, check 1d stratification properties
        all(isone, [count(x->fallsin(n*x, i), s[dim_idx, :]) for i in 1:n])
    end

    all(1:d) do dim_idx  # for every dimension, check 2d stratification properties
        all(isone,
            [count(x->fallsin(d*x, i, j), s[dim_idx, :]) for i in 1:d for j in 1:d]
        )
    end
end

@testset "LatticeRuleSample" begin
    #LatticeRuleSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, LatticeRuleSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "LDS" begin
    #LDS
    s = QuasiMonteCarlo.sample(n, lb, ub, LowDiscrepancySample([2, 3], false))
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test s[1, :] ≈ [0.5, 0.25, 0.75, 0.125, 0.625]
    @test s[2, :] ≈ [1 / 3, 2 / 3, 1 / 9, 4 / 9, 7 / 9]

    s = QuasiMonteCarlo.sample(n, Int.(lb), Int.(ub), LowDiscrepancySample([2, 3], false))
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test s[1, :] ≈ [0.5, 0.25, 0.75, 0.125, 0.625]
    @test s[2, :] ≈ [1 / 3, 2 / 3, 1 / 9, 4 / 9, 7 / 9]

    s = QuasiMonteCarlo.sample(n, zeros(Float32, 2), ones(Float32, 2),
                               LowDiscrepancySample([2, 3], false))
    @test isa(s, Matrix{Float32})
    @test size(s) == (d, n)
    @test s[1, :]≈[0.5, 0.25, 0.75, 0.125, 0.625] rtol=1e-7
    @test s[2, :]≈[1 / 3, 2 / 3, 1 / 9, 4 / 9, 7 / 9] rtol=1e-7

    testsample = []
    for i in 1:1000000
        push!(testsample,
              mean(QuasiMonteCarlo.sample(n, lb, ub, LowDiscrepancySample([2, 3], true))))
    end
    @test round(mean(testsample), sigdigits = 1)≈0.5 rtol=1e-7
end

@testset "Distribution 1" begin
    #Distribution 1
    s = QuasiMonteCarlo.sample(n, d, Cauchy())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "Distribution 2" begin
    #Distribution 2
    s = QuasiMonteCarlo.sample(n, d, Normal(3, 5))
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "GoldenSample" begin
    s = QuasiMonteCarlo.sample(n, lb, ub, GoldenSample())
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    for id in 1:d
        @test all(lb[id] .≤ s[id, :] .≤ ub[id])
    end
end

@testset "Kronecker" begin
    s = QuasiMonteCarlo.sample(n, lb, ub, KroneckerSample([sqrt(2), 3.1415], [0, 0]))
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
end

@testset "Section Sample" begin
    lb = [0.0, 0.0]
    ub = [1.0, 1.0]
    n = 5
    d = 2
    constrained_val = 1.0
    sampler = SectionSample([NaN64, constrained_val], UniformSample())
    s = QuasiMonteCarlo.sample(n, lb, ub, sampler)
    @test all(s[2, :] .== constrained_val)
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test all(lb[1] .≤ s[1, :] .≤ ub[1])
    @test QuasiMonteCarlo.fixed_dimensions(sampler) == [2]
    @test QuasiMonteCarlo.free_dimensions(sampler) == [1]

    constrained_val = 2.0
    d = 3
    sampler = SectionSample([NaN64, constrained_val, NaN64], SobolSample())
    @test QuasiMonteCarlo.free_dimensions(sampler) == [1, 3]
    @test QuasiMonteCarlo.fixed_dimensions(sampler) == [2]
    lb = [-1.0, constrained_val, -1.0]
    ub = [6.0, constrained_val, 6.0]
    s = QuasiMonteCarlo.sample(n, lb, ub, sampler)

    @test all(lb[1] .≤ s[1, :] .≤ ub[1])
    @test all(lb[3] .≤ s[3, :] .≤ ub[3])
    @test all(s[2, :] .== constrained_val)
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
end

@testset "generate_design_matrices" begin
    num_mat = 3
    Ms = QuasiMonteCarlo.generate_design_matrices(n, lb, ub, UniformSample(), num_mat)
    @test length(Ms) == num_mat
    A = Ms[1]
    B = Ms[2]
    @test isa(A, Matrix{typeof(A[1][1])}) == true
    @test size(A) == (d, n)
    @test isa(B, Matrix{typeof(B[1][1])}) == true
    @test size(B) == (d, n)
end
