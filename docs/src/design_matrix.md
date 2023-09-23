# [Design Matrices](@id DesignMatrices)

## API

It is often convenient to generate multiple independent sequences, for error estimation (uncertainty quantification).
The resulting sequences can be stored in what is often called a design matrix.
In this package, this is achieved with the `generate_design_matrices(n, d, ::DeterministicSamplingAlgorithm), ::RandomizationMethod, num_mats)` function. `num_mats` is the number of independent realizations. The resulting design matrix is a vector of matrix of length `num_mats`.

```@docs
QuasiMonteCarlo.generate_design_matrices
```

Instead of generating `num_mats` matrices, it is possible (and more memory efficient) to randomize the same matrix multiple times and perform an operation after each randomization. This is possible using iterators.
To build an iterator use the `DesignMatrix` function.

```@docs
QuasiMonteCarlo.DesignMatrix
```

!!! warning
    
    The method `generate_design_matrices(n, d, sampler, R::NoRand, num_mats, T = Float64)` is an ad hoc way to produce a Design Matrix. Indeed, it creates a deterministic point set in dimension `d × num_mats` and splits it into `num_mats` point set of dimension `d`. The resulting sequences have no QMC guarantees.
    This seems to have been proposed in Section 5.1 of [*Saltelli, A. et al. (2010)*](https://www.sciencedirect.com/science/article/pii/S0010465509003087)[^1] to do uncertainty quantification.
    See [this discussion](https://github.com/SciML/QuasiMonteCarlo.jl/pull/79#discussion_r1222946133) for a visual proof.

[^1]: Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., & Tarantola, S. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Computer physics communications, 181(2), 259-270.
## Example

```@example 2
using QuasiMonteCarlo, Random, StatsBase
Random.seed!(1234)
m = 4
d = 3
b = QuasiMonteCarlo.nextprime(d)
N = b^m # Number of points
pad = 2m # Can also choose something as `2m` to get "better" randomization
num_mats = 5

f(x) = prod(x) * 2^length(x) # test function ∫f(x)dᵈx = 1

# Randomize over num_mats = 5 independent Randomized Faure sequences
iterator = DesignMatrix(N, d, FaureSample(R = OwenScramble(base = b, pad = pad)), num_mats)

μ = [mean(f(c) for c in eachcol(X)) for X in iterator]
```

Using `std(μ)` then gives you the estimated variance of your RQMC prediction.

```@example 2
# Or using `generate_design_matrices`. Note that this is less memory efficient since it allocate space for 5 large big matrices.
μ = [mean(f(c) for c in eachcol(X)) for X in QuasiMonteCarlo.generate_design_matrices(N,
    d,
    FaureSample(R = OwenScramble(base = b, pad = pad)),
    num_mats)]
```
