using SciMLTesting, QuasiMonteCarlo

# Extension modules only exist once their trigger package is loaded, and
# ExplicitImports can only scan modules that exist, so the weakdeps have to be
# loaded here for QuasiMonteCarloDistributionsExt to be covered by the checks.
using Distributions

# QuasiMonteCarlo (own package): an extension is part of the package, but it is
# loaded as its own top-level module, so it does not share a `Base.moduleroot`
# with QuasiMonteCarlo and ExplicitImports' `allow_internal_accesses` escape
# hatch does not apply. These are QuasiMonteCarlo's own internals with no public
# spelling, used by QuasiMonteCarloDistributionsExt to implement `sample` and
# `DesignMatrix` for `Distributions.Sampleable`.
ei_kwargs = (;
    all_qualified_accesses_are_public = (;
        ignore = (
            :DistributionDesignMat, :ZERO_SAMPLES_MESSAGE,
            :initialize,
        ),
    ),
)

run_qa(QuasiMonteCarlo; ei_kwargs)
