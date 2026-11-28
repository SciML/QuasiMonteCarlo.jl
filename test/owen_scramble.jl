using QuasiMonteCarlo, Random, Statistics, Test
using Supposition: Data, @check

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
    @test QuasiMonteCarlo.splitmix64(0x0) == 0xe220a8397b1dcdaf
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
    # Digits on both sides of `digit + root shift == base` must keep distinct images.
    huge_base = typemax(UInt64)
    root = QuasiMonteCarlo.root_node(rand(Xoshiro(1), UInt64), 1)
    root_shift = QuasiMonteCarlo.node_shift(root, huge_base)
    @test root_shift >= 2
    @test allunique(scramble_digits(UInt64[huge_base - root_shift huge_base - root_shift + 1], huge_base, 1))
    @test QuasiMonteCarlo.add_mod(huge_base - 1, huge_base - 1, huge_base) == huge_base - 2
    @test QuasiMonteCarlo.add_mod(UInt8(3), UInt8(4), UInt8(5)) == 2
    @test QuasiMonteCarlo.add_mod(UInt8(1), UInt8(2), UInt8(5)) == 3
    @test_throws DomainError QuasiMonteCarlo.add_mod(UInt8(5), UInt8(0), UInt8(5))
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

@testset "OwenScramble: any point count, and extension keeps drawn points" begin
    points = QuasiMonteCarlo.sample(1000, 3, SobolSample())
    scramble(seed::Integer) = OwenScramble(base = 2, pad = 32, rng = Xoshiro(seed))
    long = randomize(points, scramble(1))
    @test randomize(points[:, 1:37], scramble(1)) == long[:, 1:37]
    @test all(0 .<= long .< 1)

    # Each scrambled point is uniform: its mean over 4000 seeds is within 4 standard errors of ½.
    draws = [randomize(points[:, 1:2], scramble(seed))[1, 2] for seed in 1:4000]
    @test abs(mean(draws) - 0.5) < 4 * sqrt(1 / 12 / length(draws))

    # A design-matrix iterator draws a new seed for each matrix.
    first_matrix, second_matrix = DesignMatrix(1000, 3, SobolSample(R = OwenScramble(base = 2)), 2)
    @test size(first_matrix) == (3, 1000)
    @test first_matrix != second_matrix
end

@testset "OwenScramble: property: scramble keeps shared leading digits" begin
    raw_digits = Data.Vectors(Data.Integers{UInt8}(); min_size = 12, max_size = 12)
    @check function keeps_shared_digits(
            base = Data.SampledFrom((2, 3, 5, 7)), pad = Data.Integers(1, 12), shared = Data.Integers(0, 12),
            x_raw = raw_digits, y_raw = raw_digits, seed = Data.Integers{UInt32}()
        )
        x = Int.(x_raw[1:pad] .% base)
        y = vcat(x[1:min(shared, pad)], Int.(y_raw[(min(shared, pad) + 1):pad] .% base))
        scrambled = scramble_digits(hcat(x, y), base, seed)
        return all(0 .<= scrambled .< base) && shared_digits(x, y) == shared_digits(scrambled[:, 1], scrambled[:, 2])
    end
end

@testset "OwenScramble: property: a prefix of the points scrambles to a prefix" begin
    points = QuasiMonteCarlo.sample(300, 3, SobolSample())
    @check function prefix_scrambles_to_prefix(count = Data.Integers(1, 300), seed = Data.Integers{UInt32}())
        scramble() = OwenScramble(base = 2, pad = 32, rng = Xoshiro(seed))
        return randomize(points[:, 1:count], scramble()) == randomize(points, scramble())[:, 1:count]
    end
end
