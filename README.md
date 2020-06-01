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
s = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LowDiscrepancySample([10,3]))
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

Additionally there is a helper function for generating design matrices:

```julia
As = QuasiMonteCarlo.generate_design_matrices(n,lb,ub,sample_method,k=2)
```

which returns `As` which is an array of `k` design matrices `A[i]` which are
all sampled from the same low discrepancy sequence.

## Available Sampling Methods

* `GridSample(dx)` where the grid is given by `lb:dx[i]:ub` in the ith direction.
* `UniformSample` for uniformly distributed random numbers.
* `SobolSample` for the Sobol sequence.
* `LatinHypercubeSample` for a Latin Hypercube.
* `LatticeRuleSample` for a randomly-shifted rank-1 lattice rule.
* `LowDiscrepancySample(base)` where `base[i]` is the base in the ith direction.
* Additionally, any `Distribution` can be used, and it will be sampled from.

## Adding a new sampling method

Adding a new sampling method is a two step process:

1. Add a new SamplingAlgorithm type
2. Overload the sample function with the new type.

**Example**

```julia
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
