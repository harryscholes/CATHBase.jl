struct MultiDomainArchitecture
    identifier::String
    domains::Vector{Domain}
end

identifier(x::MultiDomainArchitecture) = x.identifier
domains(x::MultiDomainArchitecture) = x.domains
mda(x::MultiDomainArchitecture) = superfamily.(domains(x))
Base.hash(x::MultiDomainArchitecture, h::UInt) = hash(:MultiDomainArchitecture, hash(identifier(x), hash(domains(x), h)))
Base.:(==)(x::MultiDomainArchitecture, y::MultiDomainArchitecture) = hash(x) == hash(y)

@inline function optionalsort!(mda::Vector{Domain})
    if length(mda) > 1 && !issorted(mda, by=centerofmass)
        sort!(mda, by=centerofmass)
    end
    return mda
end

function mda(mda_file::MDAFile)
    mdas = MultiDomainArchitecture[]
    buffer = sizehint!(Domain[], 1000)

    readstream(pathof(mda_file)) do io
        # parse first line
        # this step is required so that `curr_id` is set to the identifier in the first line
        line = readline(io)
        while startswith(line, "#")
            line = readline(io)
        end
        l = split(line, limit=8)
        curr_id = l[3]; sf = l[2]; boundary = l[7]
        push!(buffer, Domain(sf, boundary))

        # parse remaining lines
        for line in eachline(io)
            startswith(line, "#") && continue
            l = split(line, limit=8)
            identifier = l[3]; sf = l[2]; boundary = l[7]
            if identifier != curr_id
                optionalsort!(buffer)
                push!(mdas, MultiDomainArchitecture(curr_id, copy(buffer)))
                empty!(buffer)
                curr_id = identifier
            end

            push!(buffer, Domain(sf, boundary))
        end

        # sort!(buffer, by=centerofmass)
        optionalsort!(buffer)
        push!(mdas, MultiDomainArchitecture(curr_id, buffer))
    end

    return mdas
end

# old and slow implementation
# function mda(mda_file::MDAFile)
#     d = DefaultDict{String,Vector{Domain}}([])
#     readstream(pathof(mda_file)) do io
#         for line in eachline(io)
#             startswith(line, "#") && continue
#             l = split(line, limit=8)
#             identifier = l[3]
#             sf = l[2]
#             boundary = l[7]
#             push!(d[identifier], Domain(sf, boundary))
#         end
#     end
#     mdas = MultiDomainArchitecture[]
#     for (identifier, domains) in d
#         sort!(domains, by=centerofmass)
#         push!(mdas, MultiDomainArchitecture(identifier, domains))
#     end
#     return mdas
# end
