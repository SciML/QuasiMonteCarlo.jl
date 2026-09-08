# Prefer LatticeRules' 3-arg constructor with an Int64 point bound so we never
# hit LatticeRules 0.0.1's broken `typemax(UInt32) + 1` path on 32-bit Julia
# (Int32 wrap → ArgumentError). Avoids mutating LatticeRules via `eval`, which
# breaks incremental compilation / precompile. Safe once LatticeRules ≥ 0.0.2
# as well (same numeric bound).
function _lattice_rule(d::Integer)
    d > 0 || throw(ArgumentError("number of dimensions d must be larger than 0"))
    nmax = Int64(typemax(UInt32)) + one(Int64)
    if d ≤ 250
        return LatticeRules.LatticeRule32(LatticeRules.CKN_250_20, d, Int64(2)^20)
    else
        return LatticeRules.LatticeRule32(LatticeRules.K_3600_32, d, nmax)
    end
end
