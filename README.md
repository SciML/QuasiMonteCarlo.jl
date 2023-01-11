# QuasiMonteCarlo.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/QuasiMonteCarlo/stable/)

[![codecov](https://codecov.io/gh/SciML/QuasiMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SciML/QuasiMonteCarlo.jl)
[![Build Status](https://github.com/SciML/QuasiMonteCarlo.jl/workflows/CI/badge.svg)](https://github.com/SciML/QuasiMonteCarlo.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
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
lb = [0.1,-0.5]
ub = [1.0,20.0]
n = 5
d = 2

s = QuasiMonteCarlo.sample(n,lb,ub,GridSample([0.1,0.5]))
s = QuasiMonteCarlo.sample(n,lb,ub,UniformSample())
s = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
s = QuasiMonteCarlo.sample(n,lb,ub,HaltonSample([10,3], false))
```

The output `s` is a matrix, so one can use things like `@uview` from
[UnsafeArrays.jl](https://github.com/oschulz/UnsafeArrays.jl) for a stack-allocated
view of the `i`th point:

```julia
using UnsafeArrays
@uview s[:,i]
```

## API

Everything has the same interface:

```julia
A = QuasiMonteCarlo.sample(n,lb,ub,sample_method)
```

where:

- `n` is the number of points to sample.
- `lb` is the lower bound for each variable. The length determines the dimensionality.
- `ub` is the upper bound.
- `sample_method` is the quasi-Monte Carlo sampling strategy.

Additionally, there is a helper function for generating design matrices:

```julia
k=2
As = QuasiMonteCarlo.generate_design_matrices(n,lb,ub,sample_method,k)
```

which returns `As` which is an array of `k` design matrices `A[i]` that are
all sampled from the same low-discrepancy sequence.

## Available Sampling Methods

* `GridSample(dx)` where the grid is given by `lb:dx[i]:ub` in the ith direction.
* `UniformSample` for uniformly distributed random numbers.
* `SobolSample` for the Sobol sequence.
* `LatinHypercubeSample` for a Latin Hypercube.
* `LatticeRuleSample` for a randomly-shifted rank-1 lattice rule.
* `HaltonSample(base)` where `base[i]` is the base in the ith direction.
* `GoldenSample` for a Golden Ratio sequence.
* `KroneckerSample(alpha, s0)` for a Kronecker sequence, where alpha is an length-d vector of irrational numbers (often sqrt(d)) and s0 is a length-d seed vector (often 0).
* `SectionSample(x0, sampler)` where `sampler` is any sampler above and `x0` is a vector of either `NaN` for a free dimension or some scalar for a constrained dimension.
* Additionally, any `Distribution` can be used, and it will be sampled from.

## Adding a new sampling method

Adding a new sampling method is a two-step process:

1. Add a new SamplingAlgorithm type.
2. Overload the sample function with the new type.

All sampling methods are expected to return a matrix with dimension `d` by `n`, where `d` is the dimension of the sample space and `n` is the number of samples.

**Example**

```julia
struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end

function sample(n,lb,ub,::NewAmazingSamplingAlgorithm)
    if lb isa Number
        ...
        return x
    else
        ...
        return reduce(hcat, x)
    end
end
```

## Randomization

Note that this feature is currently experimental and is thus subject to interface changes in
non-breaking (minor) releases.

Given a matrix `x` of size `d×n` and `xᵢₛ∈[0,1]ᵈ` one obtain a randomized version `y` using one the following methods
* `owen_scramble(x, b; pad = pad)` where `b` is the base used to scramble and `pad` the number of bits in base `b` used to represent digits.
* `matousek_scramble(x, base; pad = pad)`.
* `digital_shift(x, base; pad = pad)`.
* `shift(x)`.

All these functions guarantee that the resulting array will have its components uniformly distributed `yᵢₛ∼𝐔([0,1]ᵈ)` (but not independent).

### Example 

Randomization of a Faure sequence with various methods.

```julia
    m = 4
    d = 3
    b = QuasiMonteCarlo.nextprime(d)
    N = b^m # Number of points
    pad = m

    # Unrandomized low discrepency sequence
    x_faure = QuasiMonteCarlo.sample(N, d, FaureSample())

    # Randomized version
    x_nus = randomize(x_faure, OwenScramble(base = b, pad = pad))
    x_lms = randomize(x_faure, MatousekScramble(base = b, pad = pad))
    x_digital_shift = randomize(x_faure, DigitalShift(base = b, pad = pad))
    x_shift = randomize(x_faure, Shift())
    x_uniform = rand(d, N) # plain i.i.d. uniform
```

```julia
using Plots
# Setting I like for plotting
default(fontfamily="Computer Modern", linewidth=1, label=nothing, grid=true, framestyle=:default)
```

Plot the resulting sequences along dimensions `1` and `3`.

```julia
begin
    d1 = 1
    d2 = 3
    sequences = [x_uniform, x_faure, x_shift, x_digital_shift, x_lms, x_nus]
    names = ["Uniform", "Faure (unrandomized)", "Shift", "Digital Shift", "Linear Matrix Scrambling", "Nested Uniform Scrambling"]
    p = [plot(thickness_scaling=1.5, aspect_ratio=:equal) for i in sequences]
    for (i, x) in enumerate(sequences)
        scatter!(p[i], x[d1, :], x[d2, :], ms=0.9, c=1, grid=false)
        title!(names[i])
        xlims!(p[i], (0, 1))
        ylims!(p[i], (0, 1))
        yticks!(p[i], [0, 1])
        xticks!(p[i], [0, 1])
        hline!(p[i], range(0, 1, step=1 / 4), c=:gray, alpha=0.2)
        vline!(p[i], range(0, 1, step=1 / 4), c=:gray, alpha=0.2)
        hline!(p[i], range(0, 1, step=1 / 2), c=:gray, alpha=0.8)
        vline!(p[i], range(0, 1, step=1 / 2), c=:gray, alpha=0.8)
    end
    plot(p..., size=(1500, 900))
end
```

![Different randomize methods of the same initial set of points](img/various_randomization.svg)

Faure nets and scrambled versions of Faure nets are digital $(t,d,m)$-net ([see this nice book by A. Owen](https://artowen.su.domains/mc/qmcstuff.pdf)). It basically means that they have strong equipartition properties.
Here we can (visually) verify that with Nested Uniform Scrambling (it also works with Linear Matrix Scrambling and Digital Shift).
You must see one point per rectangle of volume $1/b^m$.

```julia
begin
    d1 = 1 
    d2 = 3
    x = x_nus
    p = [plot(thickness_scaling=1.5, aspect_ratio=:equal) for i in 0:m]
    for i in 0:m
        j = m - i
        xᵢ = range(0, 1, step=1 / b^(i))
        xⱼ = range(0, 1, step=1 / b^(j))
        scatter!(p[i+1], x[d1, :], x[d2, :], ms=2, c=1, grid=false)
        xlims!(p[i+1], (0, 1.01))
        ylims!(p[i+1], (0, 1.01))
        yticks!(p[i+1], [0, 1])
        xticks!(p[i+1], [0, 1])
        hline!(p[i+1], xᵢ, c=:gray, alpha=0.2)
        vline!(p[i+1], xⱼ, c=:gray, alpha=0.2)
    end
    plot(p..., size=(1500, 900))
end
```

![All the different elementary rectangle contain only one points](img/equidistribution.svg)
