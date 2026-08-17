using QuasiMonteCarlo
using Test

struct ConstantSampler <: QuasiMonteCarlo.RandomSamplingAlgorithm end

struct GenericDeterministicSampler <: QuasiMonteCarlo.DeterministicSamplingAlgorithm
    R::RandomizationMethod
end

GenericDeterministicSampler() = GenericDeterministicSampler(Shift())

hasdoc(mod::Module, name::Symbol) = haskey(Base.Docs.meta(mod), Base.Docs.Binding(mod, name))

function QuasiMonteCarlo.sample(
        n::Integer, d::Integer, ::ConstantSampler, T = Float64
    )
    n > 0 || throw(ArgumentError("n must be positive"))
    d > 0 || throw(ArgumentError("d must be positive"))
    return fill(convert(T, 0.5), d, n)
end

@testset "Public sampling APIs" begin
    @test isdefined(QuasiMonteCarlo, :sample)
    @static if VERSION >= v"1.11"
        @test Base.ispublic(QuasiMonteCarlo, :sample)
    else
        @test Base.isexported(QuasiMonteCarlo, :sample)
    end
    @test hasdoc(QuasiMonteCarlo, :sample)

    @test isdefined(QuasiMonteCarlo, :generate_design_matrices)
    @test Base.isexported(QuasiMonteCarlo, :generate_design_matrices)
    @test hasdoc(QuasiMonteCarlo, :generate_design_matrices)

    unit_points = sample(3, 2, ConstantSampler(), Float32)
    @test size(unit_points) == (2, 3)
    @test eltype(unit_points) == Float32
    @test all(0 .<= unit_points .<= 1)

    bounded_points = sample(3, [-2.0, 4.0], [2.0, 6.0], ConstantSampler())
    @test bounded_points == [-0.0 -0.0 -0.0; 5.0 5.0 5.0]

    matrices = generate_design_matrices(3, 2, RandomSample(), 2, Float32)
    @test length(matrices) == 2
    @test all(size(matrix) == (2, 3) for matrix in matrices)
    @test all(eltype(matrix) == Float32 for matrix in matrices)

    bounded_matrices = generate_design_matrices(
        3, [-2.0, 4.0], [2.0, 6.0], RandomSample(), 2
    )
    @test all(all([-2.0, 4.0] .<= matrix .<= [2.0, 6.0]) for matrix in bounded_matrices)
end

@testset "Developer interface metadata" begin
    for name in (
            :RandomSamplingAlgorithm, :DeterministicSamplingAlgorithm,
            :AbstractDesignMatrix,
        )
        @test isdefined(QuasiMonteCarlo, name)
        @test hasdoc(QuasiMonteCarlo, name)
        @static if VERSION >= v"1.11"
            @test Base.ispublic(QuasiMonteCarlo, name)
        end
    end
end

function QuasiMonteCarlo.sample(
        n::Integer, d::Integer, ::GenericDeterministicSampler, T = Float64
    )
    return fill(convert(T, 0.25), d, n)
end

@testset "Generic sampler contracts" begin
    sampler = GenericDeterministicSampler()

    unit_points = sample(3, 2, sampler, Float32)
    @test size(unit_points) == (2, 3)
    @test eltype(unit_points) == Float32
    @test all(0 .<= unit_points .<= 1)

    bounded_points = sample(3, [-2.0, 4.0], [2.0, 6.0], sampler)
    @test size(bounded_points) == (2, 3)
    @test all([-2.0, 4.0] .<= bounded_points .<= [2.0, 6.0])

    matrices = generate_design_matrices(3, 2, sampler, 2, Float32)
    @test length(matrices) == 2
    @test all(size(matrix) == (2, 3) for matrix in matrices)
    @test all(eltype(matrix) == Float32 for matrix in matrices)

    iterator = DesignMatrix(3, 2, sampler, 2, Float32)
    @test length(iterator) == 2
    @test all(size(matrix) == (2, 3) for matrix in iterator)
end

@testset "Developer sequence validation API" begin
    @test isdefined(QuasiMonteCarlo, :_check_sequence)
    @test !Base.isexported(QuasiMonteCarlo, :_check_sequence)
    @static if VERSION >= v"1.11"
        @test Base.ispublic(QuasiMonteCarlo, :_check_sequence)
    end
    @test hasdoc(QuasiMonteCarlo, :_check_sequence)

    @test isnothing(QuasiMonteCarlo._check_sequence(1))
    @test isnothing(QuasiMonteCarlo._check_sequence([0.0, -1.0], [1.0, 1.0], 8))
    @test_throws AssertionError QuasiMonteCarlo._check_sequence(0)
    @test_throws AssertionError QuasiMonteCarlo._check_sequence([0.0], [1.0, 2.0], 8)
    @test_throws AssertionError QuasiMonteCarlo._check_sequence([0.0, 2.0], [1.0, 1.0], 8)
end
