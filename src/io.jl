abstract type AbstractFile end

abstract type AbstractMatchFile <: AbstractFile end

"""
A `.domtbl` file output by `hmmscan`.
"""
struct DOMTBLFile <: AbstractMatchFile
    path::String
end

"""
A `.crh` file output by `cath-resolve-hits`.
"""
struct CRHFile <: AbstractMatchFile
    path::String
end

"""
A `.mda` file output by `assign_cath_superfamilies.py` after processing a `.crh` file.
"""
struct MDAFile <: AbstractMatchFile
    path::String
end

"""
    pathof(fileobj::AbstractFile)

Return the path of `fileobj`.
"""
Base.pathof(x::AbstractFile) = x.path

"""
    cut(f, path; [dlm])
    cut(path, columns...; [dlm])

Extract columns from lines in a delimited file `path` using a function `f` or column numbers
`columns`.
"""
function cut(f::Function, path::AbstractString; dlm=isspace, keepempty=false)
    readstream(path) do io
        return map(eachline(io)) do line
            f(string.(split(line, dlm, keepempty=keepempty)))
        end
    end
end

function cut(path::AbstractString, columns::Integer...; dlm=isspace)
    cut(path; dlm=dlm) do line
        map(i->getindex(line, i), columns)
    end
end

cut(path::AbstractString, column::Integer; dlm=isspace) = cut(x->x[column], path; dlm=dlm)
cut(path::AbstractString, columns::AbstractRange; dlm=isspace) = cut(path, columns...; dlm=dlm)

function grep(rx::Regex, path::AbstractString)
    matches = String[]
    readstream(path) do io
        for line in eachline(io)
            if occursin(rx, line)
                push!(matches, line)
            end
        end
    end
    return matches
end

grep(pattern::AbstractString, file::AbstractString) = grep(Regex(pattern), file)

function readcodec(path::AbstractString)
    return if endswith(path, ".gz")
        GzipDecompressorStream
    # Additional required codecs should be added here
    else
        NoopStream
    end
end

function writecodec(path::AbstractString)
    return if endswith(path, ".gz")
        GzipCompressorStream
    # Additional required codecs should be added here
    else
        NoopStream
    end
end

function readstream(f::Function, path::AbstractString, args...; kwargs...)
    open(f, readcodec(path), path, args...; kwargs...)
end

function writestream(f::Function, path::AbstractString, args...; kwargs...)
    open(f, writecodec(path), path, "w", args...; kwargs...)
end

"""
    readfasta([f,] fasta)

Read a FASTA file `fasta`. If a function `f` is supplied, only records for which the
function evaluates to `true` are returned.
"""
function readfasta(f::Function, fasta::AbstractString)
    readstream(fasta) do io
        return readfasta(f, io)
    end
end

readfasta(f::Function, io::IO) = readfasta(f, FASTA.Reader(io))

function readfasta(f::Function, reader::FASTA.Reader)
    records = FASTA.Record[]
    record = FASTA.Record()
    while !eof(reader)
        read!(reader, record)
        if f(record)
            push!(records, copy(record))
        end
    end
    return records
end

readfasta(path::AbstractString) = readfasta(_->true, path)
readfasta(reader::FASTA.Reader) = readfasta(_->true, reader)

"""
    writefasta(path, records)

Write a vector of FASTA records `records` to FASTA file `path`.
"""
function writefasta(path::AbstractString, records::AbstractVector{FASTA.Record})
    writestream(path) do io
        return writefasta(io, records)
    end
end

writefasta(io::IO, records::AbstractVector{FASTA.Record}) = writefasta(FASTA.Writer(io), records)

function writefasta(writer::FASTA.Writer, records::AbstractVector{FASTA.Record})
    for record in records
        write(writer, record)
    end
    return
end

"""
    parse_sf_tsv(path)

Parse a superfamily `.tsv` file `path` to a `Vector{FASTA.Record}`, where `path` takes the
form:

    sequence_id    domain_id    mda    sequence
"""
parse_sf_tsv(path::AbstractString) = map(x->FASTA.Record(x...), cut(l->(l[2], l[4]), path))
