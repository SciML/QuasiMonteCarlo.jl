# [Developer API](@id DeveloperAPI)

!!! warning "Developer API"
    The interfaces on this page are public and versioned for package developers
    building integrations with QuasiMonteCarlo.jl. They are not exported and are
    not intended as the ordinary user-facing sampling API. Use [`sample`](@ref)
    for application code.

## Sequence validation

The shared sequence validator lets downstream sampler implementations apply the
same bounds and sample-count requirements as QuasiMonteCarlo.jl.

```@docs
QuasiMonteCarlo._check_sequence
```

## Extension interfaces

The following interfaces are public for package developers. They describe the
minimum contracts used by the generic sampling and design-matrix functions.
They are not exported, so application code should use the exported
[`sample`](@ref), [`randomize`](@ref), [`DesignMatrix`](@ref), and
[`generate_design_matrices`](@ref) functions instead of calling implementation
hooks directly.

```@docs
QuasiMonteCarlo.RandomSamplingAlgorithm
QuasiMonteCarlo.DeterministicSamplingAlgorithm
QuasiMonteCarlo.AbstractDesignMatrix
```

### Generic contract example

An extension only needs to implement the unit-box `sample` method. The bounds
and design-matrix methods are supplied by QuasiMonteCarlo.jl:

```julia
struct CenterSampler <: QuasiMonteCarlo.SamplingAlgorithm end

function QuasiMonteCarlo.sample(n::Integer, d::Integer,
        ::CenterSampler, T = Float64)
    return fill(convert(T, 0.5), d, n)
end

points = QuasiMonteCarlo.sample(4, [-1.0, 2.0], [1.0, 4.0], CenterSampler())
```

The implementation must return a `d`-by-`n` matrix in the unit box. It must
not add a competing bounds method, because the generic bounds method validates
and scales the unit-box result.
