using PrecompileTools: @compile_workload, @setup_workload

@setup_workload begin
    # Minimal precompilation workload to reduce time-to-first-sample
    # Focus on the most common use cases with Float64 and Vector{Float64} bounds

    @compile_workload begin
        # Common sample sizes and dimensions for precompilation
        n = 100
        lb = [0.0, 0.0]
        ub = [1.0, 1.0]
        d = 2

        # Precompile the main sampling methods with bounds (most common usage)
        # SobolSample - most used quasi-random sampler
        sample(n, lb, ub, SobolSample())

        # LatinHypercubeSample - very common for sensitivity analysis
        sample(n, lb, ub, LatinHypercubeSample())

        # HaltonSample - another common quasi-random sequence
        sample(n, lb, ub, HaltonSample())

        # GridSample - simple grid sampling
        sample(n, lb, ub, GridSample())

        # RandomSample - basic uniform random sampling
        sample(n, lb, ub, RandomSample())

        # LatticeRuleSample - lattice-based sampling. LatticeRules bounds the point
        # count by `typemax(UInt32) + 1`, which promotes to UInt32 and wraps to 0
        # where `Int` is 32 bits, so every lattice rule throws there.
        if Sys.WORD_SIZE > 32
            sample(n, lb, ub, LatticeRuleSample())
        end

        # GoldenSample (Kronecker) - golden ratio based sampling
        sample(n, lb, ub, GoldenSample())

        # KroneckerSample with explicit dimension
        sample(n, lb, ub, KroneckerSample(d))

        # FaureSample - requires n to be a power of base (base = nextprime(d) = 2 for d=2)
        # Use n=64 which is 2^6
        sample(64, lb, ub, FaureSample())

        # Unit box sampling (without bounds) - also common
        # These are important for users who work directly with unit hypercube
        sample(n, d, SobolSample())
        sample(n, d, LatinHypercubeSample())
        sample(n, d, HaltonSample())
        sample(n, d, GridSample())
        sample(n, d, RandomSample())
        if Sys.WORD_SIZE > 32
            sample(n, d, LatticeRuleSample())
        end
        sample(n, d, GoldenSample())
        sample(n, d, KroneckerSample(d))
        sample(64, d, FaureSample())
        # VanDerCorputSample is 1D only
        sample(n, 1, VanDerCorputSample(base = 2))
        # RandomizedHaltonSample has significant compilation time
        sample(n, d, RandomizedHaltonSample())

        # Precompile with Shift randomization (common randomization method)
        sample(n, lb, ub, SobolSample(R = Shift()))
        sample(n, lb, ub, HaltonSample(R = Shift()))
    end
end
