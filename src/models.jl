function modelsof(superfamily::SuperFamily)
    models = String[]
    for line = eachline(joinpath(@__DIR__, "..", "data", "model_to_family_map.csv"))
        l = split(line, ",")
        id = replace(l[1], "\""=>"")
        sf = replace(l[2], "\""=>"")
        if sf == identifier(superfamily)
            push!(models, id)
        end
    end
    return models
end

function map_model_to_superfamily()
    d = Dict{String,SuperFamily}()
    for line = eachline(joinpath(@__DIR__, "..", "data", "model_to_family_map.csv"))
        l = split(line, ",")
        id = replace(l[1], "\""=>"")
        sf = replace(l[2], "\""=>"")
        d[id] = SuperFamily(sf)
    end
    return d
end
