using QuasiMonteCarlo, Distributions, Test
using Test

#1D
lb = 0.0
ub = 5.0
n = 5
d = 1
QuasiMonteCarlo.sample(n,lb,ub,GridSample(0.1))
QuasiMonteCarlo.sample(n,lb,ub,UniformSample())
QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
QuasiMonteCarlo.sample(20,lb,ub,LowDiscrepancySample(10))
QuasiMonteCarlo.sample(5,d,Cauchy())
QuasiMonteCarlo.sample(5,d,Normal(0,4))

#ND
lb = [0.1,-0.5]
ub = [1.0,20.0]
n = 5
d = 2

@testset "GridSample" begin
    #GridSample{T}
    s = QuasiMonteCarlo.sample(n,lb,ub,GridSample([0.1,0.5]))
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "UniformSample" begin
    #UniformSample()
    s = QuasiMonteCarlo.sample(n,lb,ub,UniformSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "SobolSample" begin
    #SobolSample()
    s = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "LHS" begin
    #LHS
    s = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "LatticeRuleSample" begin
    #LatticeRuleSample()
    s = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "LDS" begin
    #LDS
    s = QuasiMonteCarlo.sample(n,lb,ub,LowDiscrepancySample([10,3]))
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "Distribution 1" begin
    #Distribution 1
    s = QuasiMonteCarlo.sample(n,d,Cauchy())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "Distribution 2" begin
    #Distribution 2
    s = QuasiMonteCarlo.sample(n,d,Normal(3,5))
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (d, n)
end

@testset "generate_design_matrices" begin
    num_mat = 3
    Ms = QuasiMonteCarlo.generate_design_matrices(n,lb,ub,UniformSample(), num_mat)
    @test length(Ms) == num_mat
    A = Ms[1]
    B = Ms[2]
    @test isa(A, Matrix{typeof(A[1][1])}) == true
    @test size(A) == (d, n)
    @test isa(B, Matrix{typeof(B[1][1])}) == true
    @test size(B) == (d, n)
end
