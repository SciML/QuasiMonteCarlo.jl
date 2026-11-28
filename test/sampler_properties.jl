using Test

# Supposition 0.3.5 fails to load on 32-bit Julia: its `choice!` passes a `UInt32` to a `UInt64`-only method.
if Sys.WORD_SIZE == 64
    include("shared/sampler_properties.jl")
else
    @test_skip "Supposition loads on 32-bit Julia"
end
