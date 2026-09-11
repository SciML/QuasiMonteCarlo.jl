using QuasiMonteCarlo, BenchmarkTools

const SUITE = BenchmarkGroup()

lb = zeros(4)
ub = ones(4)
n = 10_000

# =============================================================================
# Sequence generation
# =============================================================================

SUITE["sample"] = BenchmarkGroup()

SUITE["sample"]["sobol"] = @benchmarkable QuasiMonteCarlo.sample(
    $n, $lb, $ub, SobolSample()
)
SUITE["sample"]["halton"] = @benchmarkable QuasiMonteCarlo.sample(
    $n, $lb, $ub, HaltonSample()
)
SUITE["sample"]["latinhypercube"] = @benchmarkable QuasiMonteCarlo.sample(
    $n, $lb, $ub, LatinHypercubeSample()
)
SUITE["sample"]["faure"] = @benchmarkable QuasiMonteCarlo.sample(
    9375, $lb, $ub, FaureSample()
)
SUITE["sample"]["kronecker"] = @benchmarkable QuasiMonteCarlo.sample(
    $n, $lb, $ub, KroneckerSample()
)
SUITE["sample"]["lattice_rule"] = @benchmarkable QuasiMonteCarlo.sample(
    $n, $lb, $ub, LatticeRuleSample()
)
