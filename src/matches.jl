abstract type AbstractMatch end

"""
    ModelMatch <: AbstractMatch

A type representing a hidden markov model match.
"""
struct ModelMatch <: AbstractMatch
    identifier::String
    model::Union{String, Missing}
    boundary::Union{DomainBoundary, Missing}
end

ModelMatch(id::AbstractString) = ModelMatch(id, missing, missing)
ModelMatch(id::AbstractString, model::AbstractString) = ModelMatch(id, model, missing)
ModelMatch(id::AbstractString, model::AbstractString, start::Integer, stop::Integer) =
    ModelMatch(id, model, DomainBoundary(start, stop))
ModelMatch(id::AbstractString, model::AbstractString, boundary::AbstractString) =
    ModelMatch(id, model, DomainBoundary(boundary))

identifier(x::ModelMatch) = x.identifier
model(x::ModelMatch) = x.model
boundary(x::ModelMatch) = x.boundary
start(x::ModelMatch) = start(boundary(x))
stop(x::ModelMatch) = stop(boundary(x))
Base.length(x::ModelMatch) = length(boundary(x))
Base.getindex(seq::BioSequence, m::ModelMatch) = getindex(seq, boundary(m))
Base.hash(x::ModelMatch, h::UInt) = hash(:ModelMatch, hash(identifier(x), hash(model(x), hash(boundary(x), h))))
Base.:(==)(x::ModelMatch, y::ModelMatch) = hash(x) == hash(y)

const ModelMatches = AbstractVector{ModelMatch}

"""
    matches(file)

Return the `hmmscan` matches in `file`. `file` can be either a `.domtbl` file output by
`hmmscan` or a `.crh` file output by `cath-resolve-hits`.
"""
function BioSequences.matches(path::AbstractString; kwargs...)
    if endswith(path, ".domtbl")
        return matches(DOMTBLFile(path); kwargs...)
    elseif endswith(path, ".crh")
        return matches(CRHFile(path))
    end
end

# function BioSequences.matches(file::DOMTBLFile; model=false)
#     matches = Set{ModelMatch}()
#     readstream(pathof(file)) do io
#         for line in eachline(io)
#             if !startswith(line, "#")
#                 l = split(line, limit=5)
#                 id = l[1]
#                 match = model ? ModelMatch(id, l[4]) : ModelMatch(id)
#                 if match ∉ matches
#                     push!(matches, match)
#                 end
#             end
#         end
#     end
#     return collect(matches)
# end

function BioSequences.matches(file::DOMTBLFile) #; model=false)
    # matches = Set{ModelMatch}()
    matches = ModelMatch[]
    readstream(pathof(file)) do io
        for line in eachline(io)
            if !startswith(line, "#")
                l = split(line)
                id = l[1]
                model = l[4]
                start = parse(Int, l[18])
                stop = parse(Int, l[19])
                match = ModelMatch(id, model, start, stop)
                # if match ∉ matches
                push!(matches, match)
                # end
            end
        end
    end
    # return collect(matches)
    return matches
end

function BioSequences.matches(file::CRHFile)
    matches = ModelMatch[]
    readstream(pathof(file)) do io
        for line in eachline(io)
            if !startswith(line, "#")
                l = split(line)
                id = l[1]
                model = l[2]
                boundaries = l[5]
                push!(matches, ModelMatch(id, model, boundaries))
            end
        end
    end
    return matches
end

function BioSequences.matches(file::MDAFile)
    matches = ModelMatch[]
    readstream(pathof(file)) do io
        for line in eachline(io)
            if !startswith(line, "#")
                l = split(line)
                sf = l[2]
                identifier = l[3]
                boundary = l[7]
                push!(matches, ModelMatch(identifier, sf, DomainBoundary(boundary)))
            end
        end
    end
    return matches
end

function dictify(matches::ModelMatches)
    d = Dict{String,Vector{ModelMatch}}()
    for m in matches
        id = identifier(m)
        if haskey(d, id)
            push!(d[id], m)
        else
            d[id] = [m]
        end
    end
    return d
end

"""
    sequences(inputfasta, matches[, outputfasta])

Filter sequences in `inputfasta` if their identifier is in `matches`. Optionally write the
filtered sequences in `outputfasta`.
"""
function sequences(inputfasta::AbstractString, matches::ModelMatches)
    matches = Set(identifier.(matches))
    seqs = FASTA.Record[]
    readstream(inputfasta) do input
        reader = FASTA.Reader(input)
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            if identifier(record) ∈ matches
                push!(seqs, copy(record))
            end
        end
    end
    return seqs
end

function sequences(inputfasta::AbstractString, matches::ModelMatches,
                   outputfasta::AbstractString)
    matches = Set(identifier.(matches))
    readstream(inputfasta) do input
        writestream(outputfasta) do output
            reader = FASTA.Reader(input)
            writer = FASTA.Writer(output)
            record = FASTA.Record()
            while !eof(reader)
                read!(reader, record)
                if identifier(record) ∈ matches
                    write(writer, record)
                end
            end
        end
    end
end

function domains(inputfasta::AbstractString, matches::ModelMatches)
    d = dictify(matches)
    doms = FASTA.Record[]
    readstream(inputfasta) do input
        reader = FASTA.Reader(input)
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            id = identifier(record)

            # parse UniProt identifier
            if (m = match(r"(sp|tr)\|(.+)\|.+", id)) isa RegexMatch
                id = m.captures[2]
            end

            if id ∈ keys(d)
                for m in d[id]
                    domain_sequence = getindex(sequence(record), m)
                    push!(doms, FASTA.Record(
                        id*"/"*string(boundary(m)), # add domain boundary to identifier
                        description(record),
                        domain_sequence
                        )
                    )
                end
            end
        end
    end
    return doms
end

function domains(inputfasta::AbstractString, matches::ModelMatches,
                 outputfasta::AbstractString)
    d = dictify(matches)
    readstream(inputfasta) do input
        writestream(outputfasta) do output
            reader = FASTA.Reader(input)
            writer = FASTA.Writer(output)
            record = FASTA.Record()
            while !eof(reader)
                read!(reader, record)
                id = identifier(record)

                # parse UniProt identifier
                if (m = match(r"(sp|tr)\|(.+)\|.+", id)) isa RegexMatch
                    id = m.captures[2]
                end

                if id ∈ keys(d)
                    for m in d[id]
                        domain_sequence = getindex(sequence(record), m)
                        write(writer, FASTA.Record(
                            id*"/"*string(boundary(m)), # add domain boundary to identifier
                            description(record),
                            domain_sequence
                            )
                        )
                    end
                end
            end
        end
    end
end
