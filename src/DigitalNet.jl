"""
    DigitalNetSample(generating_matrices::AbstractMatrix{<:Unsigned}; R::RandomizationMethod = NoRand()) <: DeterministicSamplingAlgorithm
    DigitalNetSample(file; R::RandomizationMethod = NoRand())

A base-2 digital net or sequence defined by its generating matrices, such as a net
constructed by [LatNet Builder](https://github.com/umontreal-simul/latnetbuilder).

Row `j` of `generating_matrices` holds the generating matrix `C_j` of dimension `j`:
entry `[j, k]` is column `k` of `C_j` as a left-aligned integer, whose most significant
bit is the first row of the column. Point `i` (counting from `0`) XORs, in each
dimension, the columns selected by the binary digits of `i`, least significant digit
first. The points come in natural order, so the first `2^m` points form a net whenever
the leading `m` columns do.

`sample(n, d, DigitalNetSample(generating_matrices))` requires
`d ≤ size(generating_matrices, 1)` and `n ≤ 2^size(generating_matrices, 2)`. Each
coordinate keeps as many leading bits as the output type `T` can represent exactly.

The second form reads the generating matrices from `file`, a path or an `IO`, in the
`dnet` format of [LDData](https://github.com/QMCSoftware/LDData) or in LatNet
Builder's `-O net` output format, which is the same without the base line.

# Fields

- `generating_matrices::AbstractMatrix{<:Unsigned}`: Left-aligned generating-matrix
  columns, one row per dimension and one column per bit of the point index.
- `R::RandomizationMethod = NoRand()`: Randomization applied to the digital net.
  Scrambles must use `base = 2`.

# Examples

```jldoctest
julia> using QuasiMonteCarlo

julia> generating_matrices = UInt32[0x80000000 0x40000000; 0x80000000 0xc0000000];

julia> sample(4, 2, DigitalNetSample(generating_matrices))
2×4 Matrix{Float64}:
 0.0  0.5  0.25  0.75
 0.0  0.5  0.75  0.25
```

References:
Dick, J., & Pillichshammer, F. (2010). *Digital Nets and Sequences: Discrepancy Theory and Quasi-Monte Carlo Integration.* Cambridge University Press.
L'Ecuyer, P., Marion, P., Godin, M., & Puchhammer, F. (2022). A tool for custom construction of QMC and RQMC point sets. In *Monte Carlo and Quasi-Monte Carlo Methods: MCQMC 2020* (pp. 51-70). Springer.
"""
@concrete struct DigitalNetSample <: DeterministicSamplingAlgorithm
    generating_matrices::AbstractMatrix{<:Unsigned}
    R::RandomizationMethod
end

function DigitalNetSample(
        generating_matrices::AbstractMatrix{<:Unsigned};
        R::RandomizationMethod = NoRand()
    )
    return DigitalNetSample(generating_matrices, R)
end

function DigitalNetSample(file::Union{AbstractString, IO}; R::RandomizationMethod = NoRand())
    lines = collect(eachline(file))
    values = filter(!isempty, [strip(first(split(line, '#'))) for line in lines])
    # LDData `dnet` files state the base before s, k and r; LatNet's `-O net` omits it.
    header = 3 + (!isempty(lines) && startswith(first(lines), r"#\s*dnet"))
    if header == 4 && parse(Int, first(values)) != 2
        throw(ArgumentError("only base-2 digital nets are supported, but the file has base $(first(values))"))
    end
    s, k, r = parse.(Int, values[(header - 2):header])
    rows = values[(header + 1):end]
    if r > 64
        throw(ArgumentError("columns with $r binary digits do not fit in 64 bits"))
    end
    if length(rows) != s
        throw(ArgumentError("the header gives $s dimensions, but the file has $(length(rows)) generating matrices"))
    end
    generating_matrices = stack(row -> parse.(UInt64, split(row)) .<< (64 - r), rows; dims = 1)
    if size(generating_matrices, 2) != k
        throw(ArgumentError("the header gives $k columns, but the generating matrices have $(size(generating_matrices, 2))"))
    end
    return DigitalNetSample(generating_matrices, R)
end

function sample(n::Integer, d::Integer, S::DigitalNetSample, T::Type = Float64)
    C = S.generating_matrices
    if n < 0
        throw(ArgumentError("number of samples must be non-negative"))
    end
    if d > size(C, 1)
        throw(ArgumentError("requested $d dimensions, but the digital net has $(size(C, 1)) generating matrices"))
    end
    if n > exp2(size(C, 2))
        throw(
            ArgumentError(
                "requested $n points, but generating matrices with $(size(C, 2)) columns " *
                    "define at most 2^$(size(C, 2)) points"
            )
        )
    end

    digits = zeros(eltype(C), d, n)
    @views for i in 1:(n - 1)
        # `i & (i - 1)` clears the lowest set bit of `i`, so that point precedes point `i`.
        digits[:, i + 1] .= digits[:, (i & (i - 1)) + 1] .⊻ C[1:d, trailing_zeros(i) + 1]
    end
    return randomize(_digits2unif.(T, digits), S.R)
end

"""
    _digits2unif(T, x::Unsigned)

The left-aligned integer `x` as a fraction in [0, 1). A float `T` keeps the leading bits
it holds exactly, which keeps the result below 1; other types keep 53 bits.
"""
function _digits2unif(::Type{T}, x::Unsigned) where {T <: AbstractFloat}
    shift = max(0, 8 * sizeof(x) - precision(T))
    return ldexp(T(x >> shift), shift - 8 * sizeof(x))
end
_digits2unif(::Type{T}, x::Unsigned) where {T} = T(_digits2unif(Float64, x))
