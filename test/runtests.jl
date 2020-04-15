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
    @test size(s) == (n, d)
end

@testset "UniformSample" begin
    #UniformSample()
    s = QuasiMonteCarlo.sample(n,lb,ub,UniformSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end

@testset "SobolSample" begin
    #SobolSample()
    s = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end

@testset "LHS" begin
    #LHS
    s = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end

@testset "LDS" begin
    #LDS
    s = QuasiMonteCarlo.sample(n,lb,ub,LowDiscrepancySample([10,3]))
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end

@testset "Distribution 1" begin
    #Distribution 1
    s = QuasiMonteCarlo.sample(n,d,Cauchy())
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end

@testset "Distribution 2" begin
    #Distribution 2
    s = QuasiMonteCarlo.sample(n,d,Normal(3,5))
    @test isa(s,Matrix{typeof(s[1][1])}) == true
    @test size(s) == (n, d)
end
