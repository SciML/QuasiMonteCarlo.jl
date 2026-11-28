using QuasiMonteCarlo, Test

@testset "Samplers reject a scramble in another base" begin
    @test_throws ArgumentError QuasiMonteCarlo.sample(8, 2, SobolSample(R = OwenScramble(base = 3)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(27, 3, FaureSample(R = OwenScramble(base = 2)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(27, 3, FaureSample(R = MatousekScramble(base = 5)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(9, 1, VanDerCorputSample(base = 3, R = DigitalShift(base = 2)))
    @test_throws ArgumentError DesignMatrix(8, 2, SobolSample(), OwenScramble(base = 3), 2)
    @test_throws ArgumentError DesignMatrix(8, 2, SobolSample(R = MatousekScramble(base = 3)), 2)
    @test size(QuasiMonteCarlo.sample(27, 3, FaureSample(R = OwenScramble(base = 3)))) == (3, 27)
    # `_vdc` returns a vector, which every scramble rejects.
    @test_broken size(QuasiMonteCarlo.sample(9, 1, VanDerCorputSample(base = 3, R = DigitalShift(base = 3)))) == (1, 9)
end

@testset "HaltonSample rejects every one-base scramble" begin
    for scramble in (OwenScramble(base = 2), MatousekScramble(base = 2), DigitalShift(base = 2))
        @test_throws ArgumentError QuasiMonteCarlo.sample(16, 2, HaltonSample(R = scramble))
        @test_throws ArgumentError DesignMatrix(16, 2, HaltonSample(R = scramble), 2)
    end
    @test size(QuasiMonteCarlo.sample(16, 2, HaltonSample(R = Shift()))) == (2, 16)
end
