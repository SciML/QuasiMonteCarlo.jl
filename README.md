# QuasiMonteCarlo.jl

[![Build Status](https://travis-ci.org/JuliaDiffEq/QuasiMonteCarlo.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/QuasiMonteCarlo.jl)

This is a lightweight package for generating Quasi Monte Carlo (QMC) samples
using various different methods.

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
s = QuasiMonteCarlo.sample(n,lb,ub,LowDiscrepancySample([10,3]))
```

## API

Everything has the same interface:

```julia
QuasiMonteCarlo.sample(n,lb,ub,sample_method)
```

where:

- `n` is the number of points to sample.
- `lb` is the lower bound for each variable. The length determines the dimensionality.
- `ub` is the upper bound.
- `sample_method` is the quasi-Monte Carlo sampling strategy.

## Available Sampling Methods

* `GridSample(dx)` where the grid is given by `lb:dx[i]:ub` in the ith direction.
* `UniformSample` for uniformly distributed random numbers.
* `SobolSample` for the Sobol sequence.
* `LatinHypercubeSample` for a Latin Hypercube.
* `LowDiscrepancySample(base)` where `base[i]` is the base in the ith direction.
* Additionally, any `Distribution` can be used, and it will be sampled from.

## Adding a new sampling method

Adding a new sampling method is a two step process:

1. Add a new SamplingAlgorithm type
2. Overload the sample function with the new type.

**Example**

```
struct NewAmazingSamplingAlgorithm{OPTIONAL} <: SamplingAlgorithm end

function sample(n,lb,ub,::NewAmazingSamplingAlgorithm)
    if lb is  Number
        ...
        return x
    else
        ...
        return Tuple.(x)
    end
end
```
