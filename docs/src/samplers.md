# Sampler APIs

## Sample

```@docs
sample
```

## Samplers

Samplers are divided into two subtypes

```julia
abstract type SamplingAlgorithm end
abstract type RandomSamplingAlgorithm <: SamplingAlgorithm end
abstract type DeterministicSamplingAlgorithm <: SamplingAlgorithm end
```

### Deterministic Sampling Algorithm

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
UniformSample
LatinHypercubeSample
```