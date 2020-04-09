struct KmerIterator{T}
    seq::T
    k::Int
    start::Int
    n::Int
end

function KmerIterator(seq, k::Integer)
    KmerIterator(seq, k, 1, length(seq))
end

Base.length(iter::KmerIterator) = iter.n - iter.k + 1
Base.eltype(iter::KmerIterator{T}) where T = T

function Base.iterate(iter::KmerIterator, state=(iter.start,))
    start = state[1]
    stop = start+iter.k-1
    stop > iter.n && return nothing
    return (iter.seq[start:stop], (start+1,))
end

function BioSequences.minhash(seq::LongAminoAcidSeq, k::Integer, s::Integer)
    heap = BinaryMaxHeap{UInt64}()
    sizehint!(heap.valtree, s)
    kmerhashes = BioSequences.kmerminhash!(seq, k, s, heap)
    length(kmerhashes) < s && error("failed to generate enough hashes")
    return MinHashSketch(kmerhashes, k)
end

function BioSequences.kmerminhash!(seq::LongAminoAcidSeq, k::Integer, s::Integer, heap::BinaryMaxHeap{UInt64})
    for kmer in KmerIterator(seq, k)
        h = hash(kmer)
        if length(heap) < s
            push!(heap, h)
        elseif h < top(heap) && h ∉ heap.valtree
            pop!(heap)
            push!(heap, h)
        end
    end
    return heap.valtree
end

const VecOrSet{T} = Union{AbstractVector{T},AbstractSet{T}}

function jaccard(a::VecOrSet{T}, b::VecOrSet{T}) where T
    length(a) == length(b) == 0 && return 1.
    return length(a ∩ b) / length(a ∪ b)
end

jaccard(a::MinHashSketch, b::MinHashSketch) = jaccard(a.sketch, b.sketch)

szymkiewicz_simpson(a::VecOrSet{T}, b::VecOrSet{T}) where T = length(a ∩ b) / minimum(length.((a, b)))
