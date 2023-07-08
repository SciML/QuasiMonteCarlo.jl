# [Design Matrices](@id DesignMatrices)

## API

It is often convenient to generate multiple independent sequences, for error estimation (uncertainty quantification).
The resulting sequences can be stored in what is often called a design matrix.
In this package, this is achieved with the `generate_design_matrices(n, d, ::DeterministicSamplingAlgorithm), ::RandomizationMethod, num_mats)` function. `num_mats` is the number of independent realizations. The resulting design matrix is a vector of matrix of length `num_mats`.

```@docs
QuasiMonteCarlo.generate_design_matrices
```

!!! warning
    
    The method `generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)` is an ad hoc way to produce a Design Matrix. Indeed, it creates a deterministic point set in dimension `d × num_mats` and splits it into `num_mats` point set of dimension `d`. The resulting sequences have no QMC guarantees.
    This seems to have been proposed in Section 5.1 of [*Saltelli, A. et al. (2010)*](https://www.sciencedirect.com/science/article/pii/S0010465509003087)[^1] to do uncertainty quantification.

[^1]: Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., & Tarantola, S. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Computer physics communications, 181(2), 259-270.
## Example

```@example 2
using QuasiMonteCarlo, Random
Random.seed!(1234)
m = 4
d = 3
b = QuasiMonteCarlo.nextprime(d)
N = b^m # Number of points
pad = m # Can also choose something as `2m` to get "better" randomization

# 5 independent Randomized Faure sequences
QuasiMonteCarlo.generate_design_matrices(N,
    d,
    FaureSample(R = OwenScramble(base = b, pad = m)),
    5)
```
