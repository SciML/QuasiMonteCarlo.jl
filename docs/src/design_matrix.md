# [Design Matrices](@id DesignMatrices)

## API

It is often convenient to generate multiple independent sequence, for error estimation (uncertainty quantification).
The resulting sequences can be stored in what is often called a design matrix.
In this package, this is achieved with the `generate_design_matrices(n, d, ::DeterministicSamplingAlgorithm), ::RandomizationMethod, num_mats)` function. `num_mats` is the number of independent realization. The resulting design matrix is a vector of matrix of length `num_mats`.

```@docs
QuasiMonteCarlo.generate_design_matrices
```

!!! warning
    
    The method `generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)` is an ad hoc way to produce a Design Matrix. Indeed, it creates a deterministic point set in dimension `d × num_mats` and splits it into `num_mats` point set of dimension `d`. The resulting sequences have no QMC guarantees.
    This seems to have been proposed in Section 5.1 of [*Saltelli, A. et al. (2010)*](https://d1wqtxts1xzle7.cloudfront.net/76482087/PUBLISHED_PAPER-libre.pdf?1639660270=&response-content-disposition=inline%3B+filename%3DVariance_based_sensitivity_analysis_of_m.pdf&Expires=1688379297&Signature=d5EIx-lAwTgdGaZwbquxELAzzZESzN3hfKE-XOlw1Zn7MkPgJ1m%7EoHhM2q0QJa4rOteYa%7E6eyJFVCcmSBWWvQBUHmik8OgopOQmpGGgUDRheg7tI1i1asKkxJ6mK4NAAX7smWW0sob8rknpbnqN8zG2XiIxDlJHj2NcWYDy1Xo0Gl2gkyHelLZHYhYJYlCWsOGzTKJrZEX5wFVdCT47C%7EpM5o5d7fb9zqjc6etMM%7EKuJVvyRilYPX5rO7cNEAPJgBXm-ODCBN8mUBHnWDcZtDqDAethmf%7E06H8za9vaYtX679v2E1ezE1FV9oA2GViViLXk94WNths7gHl5nmxw2Kw__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA)[^1] to do uncertainty quantification.

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
