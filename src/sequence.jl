abstract type AbstractPartialSequence end
abstract type AbstractTruncatedSequence <: AbstractPartialSequence end

struct FullLengthSequence <: AbstractPartialSequence end
struct TruncatedSequence <: AbstractTruncatedSequence end
struct NTerminalTruncatedSequence <: AbstractTruncatedSequence end
struct CTerminalTruncatedSequence <: AbstractTruncatedSequence end
struct NorCTerminalTruncatedSequence <: AbstractTruncatedSequence end
struct NandCTerminalTruncatedSequence <: AbstractTruncatedSequence end

partialcode(::FullLengthSequence) = "00"
partialcode(::NTerminalTruncatedSequence) = "10"
partialcode(::CTerminalTruncatedSequence) = "01"
partialcode(::NandCTerminalTruncatedSequence) = "11"

function Base.filter!(records::Vector{FASTA.Record}, p::AbstractPartialSequence)
    pl = Regex("PL="*partialcode(p))
    return filter!(rec->occursin(pl, description(rec)), records)
end

function Base.filter!(records::Vector{FASTA.Record}, p::NorCTerminalTruncatedSequence)
    pl_n = Regex("PL="*partialcode(NTerminalTruncatedSequence()))
    pl_c = Regex("PL="*partialcode(CTerminalTruncatedSequence()))
    return filter!(records) do rec
        d = description(rec)
        occursin(pl_n, d) || occursin(pl_c, d)
    end
end

function Base.filter!(records::Vector{FASTA.Record}, p::TruncatedSequence)
    pl_n = Regex("PL="*partialcode(NTerminalTruncatedSequence()))
    pl_c = Regex("PL="*partialcode(CTerminalTruncatedSequence()))
    pl_nc = Regex("PL="*partialcode(NandCTerminalTruncatedSequence()))
    return filter!(records) do rec
        d = description(rec)
        occursin(pl_n, d) || occursin(pl_c, d) || occursin(pl_nc, d)
    end
end

Base.filter(records::Vector{FASTA.Record}, p::AbstractPartialSequence) = filter!(deepcopy(records), p)

function Base.map(f::Function, records::Vector{FASTA.Record}, p::AbstractPartialSequence)
    pl = Regex("PL="*CATHBase.partialcode(p))
    xs = map(records) do rec
        if occursin(pl, description(rec))
            f(rec)
        else
            missing
        end
    end
    filter!(!ismissing, xs)
    return convert(Vector{nonmissingtype(eltype(xs))}, xs)
end

function Base.map(f::Function, records::Vector{FASTA.Record}, p::NorCTerminalTruncatedSequence)
    pl_n = Regex("PL="*partialcode(NTerminalTruncatedSequence()))
    pl_c = Regex("PL="*partialcode(CTerminalTruncatedSequence()))
    xs = map(records) do rec
        if occursin(pl_n, description(rec)) || occursin(pl_c, description(rec))
            f(rec)
        else
            missing
        end
    end
    filter!(!ismissing, xs)
    return convert(Vector{nonmissingtype(eltype(xs))}, xs)
end

function biome_counts(fasta_file::AbstractString)
    xs = zeros(Int, 13)
    readstream(fasta_file) do io
        for line in eachline(io)
            if startswith(line, r"^>")
                if (m = match(r"BIOMES=(\d{13})", line)) isa RegexMatch
                    xs += parse.(Int, collect(m.captures[1]))
                end
            end
        end
    end
    return xs
end
