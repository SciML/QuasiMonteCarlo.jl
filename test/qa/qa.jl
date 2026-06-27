using SciMLTesting, QuasiMonteCarlo, Test

run_qa(
    QuasiMonteCarlo;
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractVecOrTuple,  # Base (internal)
                :GLOBAL_RNG,          # Random (non-public)
                :default_rng,         # Random (non-public)
                :skip!,               # Sobol (non-public)
            ),
        ),
    ),
)
