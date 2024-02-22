# QuasiMonteCarlo.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/QuasiMonteCarlo/stable/)

[![codecov](https://codecov.io/gh/SciML/QuasiMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/QuasiMonteCarlo.jl)
[![Build Status](https://github.com/SciML/QuasiMonteCarlo.jl/workflows/CI/badge.svg)](https://github.com/SciML/QuasiMonteCarlo.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

This is a lightweight package for generating Quasi-Monte Carlo (QMC) samples
using various different methods.

## Tutorials and Documentation

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/QuasiMonteCarlo/stable/). Use the
[in-development documentation](https://docs.sciml.ai/QuasiMonteCarlo/dev/) for the version of
the documentation, which contains the unreleased features.

## Example

```julia
using QuasiMonteCarlo, Distributions
lb = [0.1, -0.5]
ub = [1.0, 20.0]
n = 5
d = 2

s = QuasiMonteCarlo.sample(n, lb, ub, GridSample())
s = QuasiMonteCarlo.sample(n, lb, ub, Uniform())
s = QuasiMonteCarlo.sample(n, lb, ub, SobolSample())
s = QuasiMonteCarlo.sample(n, lb, ub, LatinHypercubeSample())
s = QuasiMonteCarlo.sample(n, lb, ub, LatticeRuleSample())
s = QuasiMonteCarlo.sample(n, lb, ub, HaltonSample())
```

The output `s` is a matrix, so one can use things like `@uview` from
[UnsafeArrays.jl](https://github.com/oschulz/UnsafeArrays.jl) for a stack-allocated
view of the `i`th point:

```julia
using UnsafeArrays
@uview s[:, i]
```

## API

Everything has the same interface:

```julia
A = QuasiMonteCarlo.sample(n, lb, ub, sample_method, output_type = Float64)
```

or to generate points directly in the unit box $[0,1]^d$

```julia
A = QuasiMonteCarlo.sample(n, d, sample_method, output_type = Float64) # = QuasiMonteCarlo.sample(n,zeros(d),ones(d),sample_method)
```

where:

  - `n` is the number of points to sample.
  - `lb` is the lower bound for each variable. The length determines the dimensionality.
  - `ub` is the upper bound.
  - `d` is the dimension of the unit box.
  - `sample_method` is the quasi-Monte Carlo sampling strategy.
  - `output_type` controls the output type, `Float64`, `Float32`, `Rational` (for exact digital net representation), etc. This feature does not yet work with every QMC sequence.

Additionally, there is a helper function for generating design matrices:

```julia
k = 2
As = QuasiMonteCarlo.generate_design_matrices(n,
    lb,
    ub,
    sample_method,
    k,
    output_type = Float64)
```

which returns `As` which is an array of `k` design matrices `A[i]` that are
all sampled from the same low-discrepancy sequence.

## Available Sampling Methods

Sampling methods `SamplingAlgorithm` are divided into two subtypes

  - `DeterministicSamplingAlgorithm`
    
      + `GridSample` for samples on a regular grid.
      + `SobolSample` for the Sobol sequence.
      + `FaureSample` for the Faure sequence.
      + `LatticeRuleSample` for a randomly-shifted rank-1 lattice rule.
      + `HaltonSample` for the Halton sequence.
      + `GoldenSample` for a Golden Ratio sequence.
      + `KroneckerSample(alpha, s0)` for a Kronecker sequence, where alpha is a length-`d` vector of irrational numbers (often `sqrt(d)`) and `s0` is a length-`d` seed vector (often `0`).

  - `RandomSamplingAlgorithm`
    
      + `UniformSample` for uniformly distributed random numbers.
      + `LatinHypercubeSample` for a Latin Hypercube.
      + Additionally, any `Distribution` can be used, and it will be sampled from.
    
    <!-- - `SectionSample(x0, sampler)` where `sampler` is any sampler above and `x0` is a vector of either `NaN` for a free dimension or some scalar for a constrained dimension. Not currently supported. -->

## Adding a new sampling method

Adding a new sampling method is a two-step process:

 1. Add a new SamplingAlgorithm type.
 2. Overload the sample function with the new type.

All sampling methods are expected to return a matrix with dimension `d` by `n`, where `d` is the dimension of the sample space and `n` is the number of samples.

**Example**

```julia
struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end

function sample(n, lb, ub, ::NewAmazingSamplingAlgorithm)
    if lb isa Number
        ...
        return x
    else
        ...
        return reduce(hcat, x)
    end
end
```

## Randomization of QMC sequences

Most of the previous methods are deterministic, i.e. `sample(n, d, Sampler()::DeterministicSamplingAlgorithm)` always produces the same sequence $X = (X_1, \dots, X_n)$.
There are two ways to obtain a randomized sequence:

  - Either directly use `QuasiMonteCarlo.sample(n, d, DeterministicSamplingAlgorithm(R = RandomizationMethod()))` or `sample(n, lb, up, DeterministicSamplingAlgorithm(R = RandomizationMethod()))`.
  - Or, given $n$ points $d$-dimensional points, all in $[0,1]^d$ one can do `randomize(X, ::RandomizationMethod())` where $X$ is a $d\times n$-matrix.

The currently available randomization methods are:

  - Scrambling methods `ScramblingMethods(b, pad, rng)` where `b` is the base used to scramble and `pad` the number of bits in `b`-ary decomposition.
    `pad` is generally chosen as $\gtrsim \log_b(n)$.
    The implemented `ScramblingMethods` are
    
      + `DigitalShift`
      + `MatousekScramble` a.k.a. Linear Matrix Scramble.
      + `OwenScramble` a.k.a. Nested Uniform Scramble is the most understood theoretically, but is more costly to operate.

  - `Shift(rng)` a.k.a. Cranley Patterson Rotation.

For numerous independent randomization, use `generate_design_matrices(n, d, ::DeterministicSamplingAlgorithm), ::RandomizationMethod, num_mats)` where `num_mats` is the number of independent `X` generated.

### Randomization Example

Randomization of a Faure sequence with various methods.

```julia
using QuasiMonteCarlo
m = 4
d = 3
b = QuasiMonteCarlo.nextprime(d)
N = b^m # Number of points
pad = m # Length of the b-ary decomposition number = sum(y[k]/b^k for k in 1:pad)

# Unrandomized (deterministic) low discrepancy sequence
x_faure = QuasiMonteCarlo.sample(N, d, FaureSample())

# Randomized version
x_nus = randomize(x_faure, OwenScramble(base = b, pad = pad)) # equivalent to sample(N, d, FaureSample(R = OwenScramble(base = b, pad = pad)))
x_lms = randomize(x_faure, MatousekScramble(base = b, pad = pad))
x_digital_shift = randomize(x_faure, DigitalShift(base = b, pad = pad))
x_shift = randomize(x_faure, Shift())
x_uniform = rand(d, N) # plain i.i.d. uniform
```

```julia
using Plots
# Setting I like for plotting
default(fontfamily = "Computer Modern",
    linewidth = 1,
    label = nothing,
    grid = true,
    framestyle = :default)
```

Plot the resulting sequences along dimensions `1` and `3`.

```julia
d1 = 1
d2 = 3 # you can try every combination of dimensions (d1, d2)
sequences = [x_uniform, x_faure, x_shift, x_digital_shift, x_lms, x_nus]
names = [
    "Uniform",
    "Faure (deterministic)",
    "Shift",
    "Digital Shift",
    "Matousek Scramble",
    "Owen Scramble"
]
p = [plot(thickness_scaling = 1.5, aspect_ratio = :equal) for i in sequences]
for (i, x) in enumerate(sequences)
    scatter!(p[i], x[d1, :], x[d2, :], ms = 0.9, c = 1, grid = false)
    title!(names[i])
    xlims!(p[i], (0, 1))
    ylims!(p[i], (0, 1))
    yticks!(p[i], [0, 1])
    xticks!(p[i], [0, 1])
    hline!(p[i], range(0, 1, step = 1 / 4), c = :gray, alpha = 0.2)
    vline!(p[i], range(0, 1, step = 1 / 4), c = :gray, alpha = 0.2)
    hline!(p[i], range(0, 1, step = 1 / 2), c = :gray, alpha = 0.8)
    vline!(p[i], range(0, 1, step = 1 / 2), c = :gray, alpha = 0.8)
end
plot(p..., size = (1500, 900))
```

![Different randomization methods of the same initial set of points](img/various_randomization.svg)
