# LatticeRules 0.0.1: on 32-bit Julia, `1` is Int32 so `typemax(UInt32) + 1`
# wraps to `0x0`, and every `LatticeRule32` construction fails the `n ≤ …`
# check (PieterjanRobbe/LatticeRules.jl#5). Redefine the broken methods until
# a fixed LatticeRules is released.
if Sys.WORD_SIZE == 32
    @eval LatticeRules begin
        LatticeRule32(z::Vector{UInt32}, s::Integer) =
            LatticeRule32(z, s, Int64(typemax(UInt32)) + one(Int64))
        function LatticeRule32(z::Vector{UInt32}, s::Integer, n::Integer)
            s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
            s ≤ length(z) || throw(
                ArgumentError(
                    "number of dimensions s must be less than or equal to the length of the generating vector z",
                ),
            )
            n > 0 || throw(ArgumentError("maximum number of points n must be larger than 0"))
            n ≤ Int64(typemax(UInt32)) + one(Int64) || throw(
                ArgumentError(
                    "maximum number of points n must be less than or equal to 2^32, consider implementing a LatticeRule64 type",
                ),
            )
            return LatticeRule32{s}(view(z, 1:s), n)
        end
    end
end
