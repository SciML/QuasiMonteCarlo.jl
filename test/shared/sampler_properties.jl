using QuasiMonteCarlo, Test, Primes
using Supposition: Data, @check

# Samplers that accept every size `n ≥ 1` in every dimension.
const ANY_SIZE_SAMPLERS = (
    GridSample(), SobolSample(), LatticeRuleSample(), HaltonSample(), GoldenSample(),
    KroneckerSample(), RandomSample(), LatinHypercubeSample(), RandomizedHaltonSample(),
)
const DETERMINISTIC_SAMPLERS = filter(s -> s isa QuasiMonteCarlo.DeterministicSamplingAlgorithm, ANY_SIZE_SAMPLERS)
const ALL_SAMPLERS = (ANY_SIZE_SAMPLERS..., FaureSample())

inside_open_box(points) = all(x -> 0 < x < 1, points)
# Uniform draws from 1:2000 rarely give the smallest sizes, where edge cases live.
const SIZES = Data.Integers(1, 8) | Data.Integers(1, 2000)
cell_centers(n) = ((0:(n - 1)) .+ 0.5) ./ n
"""Whether `run()` throws an `E`. Any other exception propagates, so the property reports its type."""
function throws(::Type{E}, run) where {E <: Exception}
    try
        run()
    catch err
        err isa E && return true
        rethrow()
    end
    return false
end

@testset "samplers: property: valid sizes give n points inside the open unit box" begin
    @check max_examples = 300 function any_size_inside(
            sampler = Data.SampledFrom(DETERMINISTIC_SAMPLERS), n = SIZES, d = Data.Integers(1, 8)
        )
        points = QuasiMonteCarlo.sample(n, d, sampler)
        return size(points) == (d, n) && inside_open_box(points)
    end
    @check max_examples = 300 function van_der_corput_inside(base = Data.Integers(2, 16), n = SIZES)
        points = QuasiMonteCarlo.sample(n, 1, VanDerCorputSample(; base))
        return length(points) == n && inside_open_box(points)
    end
    # Faure accepts `n = copies * base^power` with `copies < base`.
    @check max_examples = 300 function faure_inside(
            d = Data.Integers(1, 12), copies = Data.Integers(1, 12), power = Data.Integers(0, 2)
        )
        base = nextprime(d)
        n = min(copies, base - 1) * base^power
        points = QuasiMonteCarlo.sample(n, d, FaureSample())
        return size(points) == (d, n) && inside_open_box(points)
    end
end

# Exact powers are rare among random sizes, so these properties draw `n = base^power` directly.
@testset "radical-inverse samplers: property: a full net puts one point at each cell center" begin
    @check max_examples = 300 function van_der_corput_centers(base = Data.Integers(2, 10), power = Data.Integers(0, 4))
        n = base^power
        return sort(QuasiMonteCarlo.sample(n, 1, VanDerCorputSample(; base))) ≈ cell_centers(n)
    end
    @check max_examples = 300 function halton_centers(
            d = Data.Integers(1, 6), coordinate = Data.Integers(1, 6), power = Data.Integers(0, 4)
        )
        coordinate = min(coordinate, d)
        n = prime(coordinate)^power
        return sort(QuasiMonteCarlo.sample(n, d, HaltonSample())[coordinate, :]) ≈ cell_centers(n)
    end
    @check max_examples = 300 function faure_centers(d = Data.Integers(1, 7), power = Data.Integers(0, 5))
        n = nextprime(d)^power
        points = QuasiMonteCarlo.sample(n, d, FaureSample())
        return all(sort(row) ≈ cell_centers(n) for row in eachrow(points))
    end
end

@testset "radical-inverse samplers: property: the index digit count covers the largest index" begin
    covers(digits, n, base) = digits >= 1 && n - 1 < big(base)^digits && (digits == 1 || big(base)^(digits - 1) <= n - 1)
    @check function digits_near_powers(base = Data.Integers(2, 36), power = Data.Integers(0, 12), step = Data.Integers(-1, 1))
        n = max(1, base^power + step)
        return covers(QuasiMonteCarlo._index_digits(n, base), n, base)
    end
    @check function digits_anywhere(base = Data.Integers(2, 36), n = Data.Integers(1, typemax(Int)))
        return covers(QuasiMonteCarlo._index_digits(n, base), n, base)
    end
end

@testset "samplers: property: invalid sizes and bounds throw argument errors" begin
    @check max_examples = 300 function nonpositive_sizes_throw(sampler = Data.SampledFrom(ALL_SAMPLERS), n = Data.Integers(-1000, 0))
        return throws(ArgumentError, () -> QuasiMonteCarlo.sample(n, zeros(2), ones(2), sampler))
    end
    @check max_examples = 300 function inverted_bounds_throw(
            sampler = Data.SampledFrom(ALL_SAMPLERS), d = Data.Integers(1, 5), flipped = Data.Integers(1, 5),
            gap = Data.Floats{Float64}(; nans = false, minimum = 1.0e-9, maximum = 1.0e9)
        )
        lb, ub = zeros(d), ones(d)
        ub[min(flipped, d)] = -gap
        return throws(ArgumentError, () -> QuasiMonteCarlo.sample(8, lb, ub, sampler))
    end
    @check max_examples = 300 function mismatched_bounds_throw(
            sampler = Data.SampledFrom(ALL_SAMPLERS), d = Data.Integers(1, 5), extra = Data.Integers(1, 3)
        )
        return throws(DimensionMismatch, () -> QuasiMonteCarlo.sample(8, zeros(d), ones(d + extra), sampler))
    end
end

@testset "samplers: property: scaled points stay inside the closed box" begin
    coordinate = Data.Floats{Float64}(; nans = false, infs = false, minimum = -1.0e6, maximum = 1.0e6)
    width = Data.Floats{Float64}(; nans = false, infs = false, minimum = 0.0, maximum = 1.0e6)
    @check max_examples = 300 function scaled_inside(
            sampler = Data.SampledFrom(ANY_SIZE_SAMPLERS), n = Data.Integers(1, 64),
            box = Data.Vectors(Data.Pairs(coordinate, width); min_size = 1, max_size = 4)
        )
        lb = first.(box)
        ub = lb .+ last.(box)
        points = QuasiMonteCarlo.sample(n, lb, ub, sampler)
        return all(lb .<= points .<= ub)
    end
end
