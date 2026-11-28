using QuasiMonteCarlo, Test

@testset "Samplers reject a scramble in another base" begin
    @test_throws ArgumentError QuasiMonteCarlo.sample(8, 2, SobolSample(R = OwenScramble(base = 3)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(27, 3, FaureSample(R = OwenScramble(base = 2)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(27, 3, FaureSample(R = MatousekScramble(base = 5)))
    @test_throws ArgumentError QuasiMonteCarlo.sample(9, 1, VanDerCorputSample(base = 3, R = DigitalShift(base = 2)))
    @test_throws ArgumentError DesignMatrix(8, 2, SobolSample(), OwenScramble(base = 3), 2)
    @test_throws ArgumentError DesignMatrix(8, 2, SobolSample(R = MatousekScramble(base = 3)), 2)
    @test size(QuasiMonteCarlo.sample(27, 3, FaureSample(R = OwenScramble(base = 3)))) == (3, 27)
end

@testset "HaltonSample rejects a one-base scramble in more than one dimension" begin
    for scramble in (OwenScramble(base = 2), MatousekScramble(base = 2), DigitalShift(base = 2))
        @test_throws ArgumentError QuasiMonteCarlo.sample(16, 2, HaltonSample(R = scramble))
        @test_throws ArgumentError DesignMatrix(16, 2, HaltonSample(R = scramble), 2)
    end
    @test size(QuasiMonteCarlo.sample(16, 2, HaltonSample(R = Shift()))) == (2, 16)
    # One coordinate is a base-2 sequence, so only a base-2 scramble fits.
    @test size(QuasiMonteCarlo.sample(16, 1, HaltonSample(R = OwenScramble(base = 2)))) == (1, 16)
    @test all(size(matrix) == (1, 16) for matrix in DesignMatrix(16, 1, HaltonSample(R = OwenScramble(base = 2)), 2))
    @test_throws ArgumentError QuasiMonteCarlo.sample(16, 1, HaltonSample(R = OwenScramble(base = 3)))
    @test_throws ArgumentError DesignMatrix(16, 1, HaltonSample(R = OwenScramble(base = 3)), 2)
end
