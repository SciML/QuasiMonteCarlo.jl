# Sampler APIs

## Sample

```@docs
QuasiMonteCarlo.sample
```

## Samplers

Samplers are divided into two subtypes

```@docs
SamplingAlgorithm
```

The extension contracts for `RandomSamplingAlgorithm` and
`DeterministicSamplingAlgorithm` are documented on the [Developer API](@ref
DeveloperAPI) page.

### Deterministic Sampling Algorithm

All `DeterministicSamplingAlgorithm` have `NoRand()` as their default `RandomizationMethod`, see [Randomization methods](@ref Randomization) and [Design Matrices section](@ref DesignMatrices) for more information on randomization.

```@docs
GridSample
```

```@docs
SobolSample
```

!!! warning
    
    The QuasiMonteCarlo.jl package relies on the [Sobol.jl](https://github.com/JuliaMath/Sobol.jl) package to sample Sobol nets. The choice, there is to NOT start the sequence at `0`. This is debatable, see this [issue](https://github.com/JuliaMath/Sobol.jl/issues/31#issuecomment-1528136486) and ref therein for more context.

```@docs
VanDerCorputSample
FaureSample
LatticeRuleSample
HaltonSample
GoldenSample
KroneckerSample
```

### Random Sampling Algorithm

```@docs
RandomSample
LatinHypercubeSample
RandomizedHaltonSample
```
