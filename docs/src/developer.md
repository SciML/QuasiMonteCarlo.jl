# Developer API

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
