using Supposition: Data, @check

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
