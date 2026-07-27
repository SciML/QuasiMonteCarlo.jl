using QuasiMonteCarlo
using Test

struct ConstantSampler <: SamplingAlgorithm end

function QuasiMonteCarlo.sample(
        n::Integer, d::Integer, ::ConstantSampler, T = Float64
    )
    n > 0 || throw(ArgumentError("n must be positive"))
    d > 0 || throw(ArgumentError("d must be positive"))
    return fill(convert(T, 0.5), d, n)
end

@testset "Public sampling API" begin
    @test isdefined(QuasiMonteCarlo, :sample)
    @test Base.isexported(QuasiMonteCarlo, :sample)
    @test Base.Docs.doc(Base.Docs.Binding(QuasiMonteCarlo, :sample)) !== nothing

    unit_points = sample(3, 2, ConstantSampler(), Float32)
    @test size(unit_points) == (2, 3)
    @test eltype(unit_points) == Float32
    @test all(0 .<= unit_points .<= 1)

    bounded_points = sample(3, [-2.0, 4.0], [2.0, 6.0], ConstantSampler())
    @test bounded_points == [-0.0 -0.0 -0.0; 5.0 5.0 5.0]
end
