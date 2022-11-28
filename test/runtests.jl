using QuasiMonteCarlo, Random
using IntervalArithmetic, Primes, Combinatorics, Distributions, StatsBase
using Test

# For testing randomized QMC sequences by using the deterministic version

# struct InertSampler <: Random.AbstractRNG end
# InertSampler(args...; kwargs...) = InertSampler()
# Random.rand(::InertSampler, ::Type{T}) where {T} = zero(T)
# Random.rand(::InertSampler) = 0
# Random.shuffle!(::InertSampler, arg::AbstractArray) = arg

function Base.resize!(a::Vector{T}, nl::Integer, pad::T) where {T}
    l = length(a)
    if nl > l
        Base._growend!(a, nl - l)
        a[(l + 1):end] .= pad
    elseif nl != l
        if nl < 0
            throw(ArgumentError("new length must be ≥ 0"))
        end
        Base._deleteend!(a, l - nl)
    end
    return a
end

#1D
lb = 0.0
ub = 5.0
n = 5
d = 1

for sampler in [GridSample(0.1), UniformSample(), SobolSample(), LatinHypercubeSample(), LatticeRuleSample(), Cauchy(), Normal(0, 4), GoldenSample()]
    @show sampler
    A = QuasiMonteCarlo.sample(n, lb, ub, sampler)
    @test all(all(x .<= ub) for x in eachcol(A))
    @test all(all(x .>= lb) for x in eachcol(A))
end

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

        #n = 0
        @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, 0.0, 1.0,
                                                                             LowDiscrepancySample(2,
                                                                                                  false))
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

        #n = 0
        @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                             sampler)
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

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         GridSample([
                                                                                        0.1,
                                                                                        0.5,
                                                                                    ]))
end

@testset "UniformSample" begin
    #UniformSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, UniformSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         UniformSample())
end

@testset "SobolSample" begin
    #SobolSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         SobolSample())
end

@testset "LHS" begin
    #LHS
    s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample(threading = true))
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         LatinHypercubeSample())
end

@testset "Faure Sample" begin
    #FaureSample()
    d = 17
    base = nextprime(d)
    power = 2
    n = 17^2
    @test_throws ArgumentError QuasiMonteCarlo.sample(d + 1, d, FaureSample())
    @test_throws ArgumentError QuasiMonteCarlo.sample(d^2 + 1, d, FaureSample())
    s = sortslices(QuasiMonteCarlo.sample(n, d, FaureSample()); dims = 2)
    # FaureSample() generates centered boxes, unlike DiceDesign
    r = sortslices(include("rfaure.jl")'; dims = 2) .+ inv(2base^(power + 1))

    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test mean(abs2, s - r) ≤ inv(2n)
    @test s ≈ r

    # check RQMC stratification properties
    # Deterministic Faure
    function test_tms(d, n, power)
        pass = true
        base = nextprime(d)
        s = QuasiMonteCarlo.sample(n, d, FaureSample())
        parts = resize!.(partitions(power, d), d, 0)
        perms = Iterators.map(parts) do x
            multiset_permutations(x, d)
        end
        for stepsize in Iterators.flatten(perms)
            intervals = mince(IntervalBox([interval(0, 1) for i in 1:d]),
                              Tuple(base .^ stepsize))
            pass = pass && all(intervals) do intvl
                count(point -> point ∈ intvl, eachslice(s; dims = 2)) == 1
            end
            if !pass
                println("Errors in dimension $d, interval $stepsize, sample size $n")
                return pass
            end
        end
        return pass
    end
    power = 5
    for d in (3, 5, 7)
        @test test_tms(d, d^power, power)  # test 5d stratification of first 3 primes
    end
    power = 3
    for d in (11, 13, 17)
        @test test_tms(d, d^power, power)  # test 3d stratification of next 3 primes
    end

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         FaureSample())
end

@testset "LatticeRuleSample" begin
    #LatticeRuleSample()
    s = QuasiMonteCarlo.sample(n, lb, ub, LatticeRuleSample())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         LatticeRuleSample())
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

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         LowDiscrepancySample([
                                                                                                  2,
                                                                                                  3,
                                                                                              ],
                                                                                              false))
end

@testset "Distribution 1" begin
    #Distribution 1
    s = QuasiMonteCarlo.sample(n, d, Cauchy())
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, d, Cauchy())
end

@testset "Distribution 2" begin
    #Distribution 2
    s = QuasiMonteCarlo.sample(n, d, Normal(3, 5))
    @test isa(s, Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, d, Normal(3, 5))
end

@testset "Kronecker" begin
    ρ = 0.7548776662466927
    s = QuasiMonteCarlo.sample(n, 2, KroneckerSample([ρ, ρ^2], [0, 0]))
    t = QuasiMonteCarlo.sample(n, 2, KroneckerSample([ρ, ρ^2]))
    s ≈ t
    t .= QuasiMonteCarlo.sample(n, 2, GoldenSample())
    @test isa(s, Matrix{Float64})
    @test size(s) == (d, n)
    @test s ≈ t
    @test all(x -> (mod(x, 1) ≈ s[1, 2] - s[1, 1]), diff(s[1, :]))

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, 2,
                                                                         GoldenSample())
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, 2,
                                                                         KroneckerSample([
                                                                                             ρ,
                                                                                             ρ^2,
                                                                                         ],
                                                                                         [
                                                                                             0,
                                                                                             0,
                                                                                         ]))
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

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.sample(0, lb, ub,
                                                                         SectionSample([
                                                                                           NaN64,
                                                                                           constrained_val,
                                                                                       ],
                                                                                       UniformSample()))
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

    #n = 0
    @test_throws QuasiMonteCarlo.ZeroSamplesError QuasiMonteCarlo.generate_design_matrices(0,
                                                                                           lb,
                                                                                           ub,
                                                                                           UniformSample(),
                                                                                           num_mat)
end
