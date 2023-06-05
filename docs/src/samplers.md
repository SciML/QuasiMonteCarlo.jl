# Sampler APIs

## Sample

```@docs
QuasiMonteCarlo.sample
```

## Samplers

Samplers are divided into two subtypes

```julia
abstract type SamplingAlgorithm end
abstract type RandomSamplingAlgorithm <: SamplingAlgorithm end
abstract type DeterministicSamplingAlgorithm <: SamplingAlgorithm end
```

### Deterministic Sampling Algorithm

All `DeterministicSamplingAlgorithm` have `NoRand()` as their default `RandomizationMethod`, see [Randomization methods](@ref Randomization) section.

```@docs
GridSample
SobolSample
FaureSample
LatticeRuleSample
HaltonSample
GoldenSample
KroneckerSample
```

### Random Sampling Algorithm

```@docs
LatinHypercubeSample
```
