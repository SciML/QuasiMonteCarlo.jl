using QuasiMonteCarlo, Random, Statistics, Test
using Primes: prevprime

"""
    shared_digits(x::AbstractVector{<:Integer}, y::AbstractVector{<:Integer}) -> Int

Return the number of leading digits that `x` and `y` share.
"""
shared_digits(x::AbstractVector{<:Integer}, y::AbstractVector{<:Integer}) = something(findfirst(x .!= y), length(x) + 1) - 1

"""
    scramble_digits(digit_columns::AbstractMatrix{<:Integer}, base::Integer, seed::Integer) -> Matrix

Scramble each column of the `pad × n` matrix `digit_columns` with one `OwenScramble` seeded by `seed`.
"""
function scramble_digits(digit_columns::AbstractMatrix{<:Integer}, base::Integer, seed::Integer)
    origin = reshape(digit_columns, size(digit_columns, 1), :, 1)
    scramble = OwenScramble(base = base, pad = oftype(base, size(digit_columns, 1)), rng = Xoshiro(seed))
    return QuasiMonteCarlo.randomize_bits!(similar(origin), origin, scramble)[:, :, 1]
end

@testset "OwenScramble: hash reference values" begin
    # First two outputs of the SplitMix64 reference generator seeded with 0.
    @test QuasiMonteCarlo.splitmix64(0x00) == 0xe220a8397b1dcdaf
    @test QuasiMonteCarlo.splitmix64(QuasiMonteCarlo.GOLDEN_GAMMA) == 0x6e789e6aa1b965f4
    @test QuasiMonteCarlo.node_shift(UInt64(0), 3) == 0
    @test QuasiMonteCarlo.node_shift(typemax(UInt64), 3) == 2
    # The gap between sibling hashes must depend on the parent; a hash linear in its key fails this.
    parents = rand(Xoshiro(1), UInt64, 1000)
    @test allunique(QuasiMonteCarlo.child_node(parent, 1) - QuasiMonteCarlo.child_node(parent, 0) for parent in parents)
end

@testset "OwenScramble: nested bijection of each digit grid" begin
    # Checked on digits, as `unif2bits` does not round-trip every base-3 fraction.
    for (base, pad) in ((2, 10), (3, 6), (5, 4))
        grid = stack(reverse(digits(k; base, pad)) for k in 0:(base^pad - 1)) # sorted columns
        scrambled = scramble_digits(grid, base, 3)
        @test sortslices(scrambled, dims = 2) == grid
        @test all(
            shared_digits(view(grid, :, i), view(grid, :, j)) == shared_digits(view(scrambled, :, i), view(scrambled, :, j))
                for i in 1:100, j in axes(grid, 2)
        )
        # Unlike a digital shift, the last digit's shift varies with the digits above it.
        @test !allequal(mod.(scrambled[end, :] .- grid[end, :], base))
        @test scrambled != scramble_digits(grid, base, 4)
    end
end

@testset "OwenScramble: digit sums do not overflow" begin
    # 126 plus a shift of 2 or more exceeds typemax(Int8).
    top_digits = fill(Int8(126), 4, 50)
    @test all(0 .<= scramble_digits(top_digits, Int8(127), 1) .< 127)
    # Digits whose scaled values straddle `scaled + root shift == base` must keep distinct images.
    huge_base = prevprime(typemax(UInt64))
    root = QuasiMonteCarlo.root_node(rand(Xoshiro(1), UInt64), 1)
    root_shift = QuasiMonteCarlo.node_shift(root, huge_base)
    inverse = invmod(big(QuasiMonteCarlo.node_multiplier(root, huge_base)), huge_base)
    straddle = [UInt64(mod(big(scaled) * inverse, huge_base)) for scaled in (huge_base - root_shift - 1, huge_base - root_shift)]
    @test root_shift >= 2
    @test scramble_digits(permutedims(straddle), huge_base, 1) == [huge_base - 1 0]
    @test QuasiMonteCarlo.add_mod(huge_base - 1, huge_base - 1, huge_base) == huge_base - 2
    @test QuasiMonteCarlo.add_mod(UInt8(3), UInt8(4), UInt8(5)) == 2
    @test QuasiMonteCarlo.add_mod(UInt8(1), UInt8(2), UInt8(5)) == 3
    @test_throws DomainError(UInt64(5), "a base-5 digit must lie in 0:4") QuasiMonteCarlo.permute_digit(UInt64(1), UInt64(5), UInt64(5))
end

@testset "OwenScramble: node permutations are affine with a uniform multiplier" begin
    @test_throws ArgumentError OwenScramble(base = 4)
    @test_throws ArgumentError OwenScramble(base = 1)
    seeds = 1:6000
    # Base 3 has 3! = 3 · 2 affine maps, so the root permutation must take all six values equally often.
    permutations = [vec(scramble_digits([0 1 2], 3, seed)) for seed in seeds]
    @test all(abs(count(==(permutation), permutations) / length(seeds) - 1 / 6) < 4 * sqrt(5 / 36 / length(seeds)) for permutation in unique(permutations))
    @test length(unique(permutations)) == 6
    # Base 5: the gap between the images of digits 0 and 1 is the multiplier, uniform on 1:4.
    gaps = [mod(-(scramble_digits([0 1], 5, seed)...), 5) for seed in seeds]
    @test all(abs(count(==(gap), gaps) / length(seeds) - 1 / 4) < 4 * sqrt(3 / 16 / length(seeds)) for gap in 1:4)
end

@testset "OwenScramble: sibling nodes draw independent permutations" begin
    # The two points share digit 1 and differ at digit 2, so digit 3 goes through two sibling nodes.
    pair = [0 0; 0 1; 0 0]
    seeds = 1:4000
    for base in (2, 3)
        scrambles = [scramble_digits(pair, base, seed) for seed in seeds]
        @test all(scrambled[1, 1] == scrambled[1, 2] for scrambled in scrambles)
        agreement = mean(scrambled[3, 1] == scrambled[3, 2] for scrambled in scrambles)
        @test abs(agreement - 1 / base) < 4 * sqrt((1 / base) * (1 - 1 / base) / length(seeds))
    end
end

@testset "OwenScramble: a child's permutation is independent of its parent's" begin
    # Digit 1 of these points gives the root's multiplier; digit 2 of point `child` is that child's shift.
    points = [0 1 2; 0 0 0]
    seeds = 1:6000
    scrambles = [scramble_digits(points, 3, seed) for seed in seeds]
    multipliers = [mod(scrambled[1, 2] - scrambled[1, 1], 3) for scrambled in scrambles]
    for child in 1:3
        cells = [(multiplier, scrambled[2, child]) for (multiplier, scrambled) in zip(multipliers, scrambles)]
        @test all(abs(count(==((multiplier, shift)), cells) / length(seeds) - 1 / 6) < 4 * sqrt(5 / 36 / length(seeds)) for multiplier in 1:2, shift in 0:2)
    end
end

@testset "OwenScramble: any point count, and extension keeps drawn points" begin
    points = QuasiMonteCarlo.sample(1000, 3, SobolSample())
    scramble(seed::Integer) = OwenScramble(base = 2, pad = 32, rng = Xoshiro(seed))
    long = randomize(points, scramble(1))
    @test randomize(points[:, 1:37], scramble(1)) == long[:, 1:37]
    @test all(0 .<= long .< 1)
    @test OwenScramble(base = 2).pad == 32
    @test_throws AssertionError QuasiMonteCarlo.randomize!(zeros(3, 10), points, scramble(1))

    # Each scrambled point is uniform: over 4000 seeds, its mean and variance lie within 4 standard errors of ½ and ¹⁄₁₂.
    # The variance check also fails if `randomize` returns its input unchanged.
    draws = [randomize(points[:, 1:3], scramble(seed))[2, 3] for seed in 1:4000]
    @test abs(mean(draws) - 1 / 2) < 4 * sqrt(1 / 12 / length(draws))
    @test abs(var(draws) - 1 / 12) < 4 * sqrt(1 / 180 / length(draws))

    # A design-matrix iterator draws a new seed for each matrix.
    first_matrix, second_matrix = DesignMatrix(1000, 3, SobolSample(R = OwenScramble(base = 2)), 2)
    @test size(first_matrix) == (3, 1000)
    @test first_matrix != second_matrix
end

# Supposition 0.3.5 fails to load on 32-bit Julia: its `choice!` passes a `UInt32` to a `UInt64`-only method.
if Sys.WORD_SIZE == 64
    include("shared/owen_properties.jl")
else
    @test_skip "Supposition loads on 32-bit Julia"
end
