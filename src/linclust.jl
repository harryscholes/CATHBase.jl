function linclust(fasta_filepath::AbstractString;
                  output_dir::AbstractString=dirname(fasta_filepath),
                  clust_file::AbstractString=replace(basename(fasta_filepath), r".fa$"=>".clust.tsv"),
                  c::AbstractFloat=0.9, min_seq_id::AbstractFloat=0.9, cov_mode::Integer=0,
                  cluster_mode::Integer=0, kmer_per_seq::Integer=20, threads::Integer=1)
    fasta_file = basename(fasta_filepath)
    clust_filepath = joinpath(output_dir, clust_file)
    isfile(fasta_filepath) || throw(ArgumentError("FASTA file does not exist: "*fasta_filepath))
    isdir(output_dir) || throw(ArgumentError("Directory does not exist: "*output_dir))
    isfile(clust_filepath) && throw(ArgumentError("Output file already exists: "*clust_filepath))
    temp_dir = mktempdir()
    try
        run(`mmseqs createdb $fasta_filepath $temp_dir/$fasta_file.db`)
        run(```
            mmseqs linclust
                -c $c
                --min-seq-id $min_seq_id
                --cov-mode $cov_mode
                --cluster-mode $cluster_mode
                --kmer-per-seq $kmer_per_seq
                --threads $threads
                $temp_dir/$fasta_file.db $temp_dir/$fasta_file.clust $temp_dir/tmp
            ```)
        run(```
            mmseqs createtsv $temp_dir/$fasta_file.db $temp_dir/$fasta_file.db
                $temp_dir/$fasta_file.clust $temp_dir/$clust_file
            ```)
        cp("$temp_dir/$clust_file", clust_filepath)
    catch e
        throw(e)
    finally
        rm(temp_dir, recursive=true)
    end
end

struct LinclustFile <: AbstractFile
    path::String
end

function Base.read(f::Function, x::LinclustFile)
    d = DefaultDict{String,Vector{String}}([])
    readstream(pathof(x)) do io
        for line in eachline(io)
            l = split(line)
            if f(l)
                push!(d[l[1]], l[2])
            end
        end
    end
    return Dict(d)
end

Base.read(x::LinclustFile) = read(l->true, x)

ismgyid(x) = startswith(x, "MGY")
isuniprotid(x) = !ismgyid(x)

ismgyids(xs) = all(ismgyid.(xs))
isuniprotids(xs) = all(isuniprotid.(xs))
ismixedids(xs) = !(ismgyids(xs) || isuniprotids(xs))

fraction_mgy(xs) = count(ismgyid.(xs)) / length(xs)
