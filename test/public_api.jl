using QuasiMonteCarlo
using Test

struct ConstantSampler <: SamplingAlgorithm end

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
