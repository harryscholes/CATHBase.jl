module CATHBase

export
    SuperFamily,
    @sf_str,
    FunctionalFamily,
    @ff_str,
    DomainSection,
    DomainBoundary,
    mass,
    centerofmass,
    ModelMatch,
    identifier,
    model,
    boundary,
    iscontinuous,
    start,
    stop,
    Domain,
    DOMTBLFile,
    CRHFile,
    MDAFile,
    pathof,
    matches,
    sequences,
    domains,
    modelsof,
    map_model_to_superfamily,
    cut,
    grep,
    readcodec,
    writecodec,
    readstream,
    writestream,
    readfasta,
    writefasta,
    parse_sf_tsv,
    linclust,
    LinclustFile,
    ismgyid,
    isuniprotid,
    ismgyids,
    isuniprotids,
    ismixedids,
    fraction_mgy,
    MultiDomainArchitecture,
    mda,
    FullLengthSequence,
    TruncatedSequence,
    NTerminalTruncatedSequence,
    CTerminalTruncatedSequence,
    NorCTerminalTruncatedSequence,
    NandCTerminalTruncatedSequence,
    biome_counts,
    minhash,
    jaccard,
    szymkiewicz_simpson,
    WeightedEdge,
    edge,
    weight,
    VertexDict,
    WeightDict,
    index,
    graph

using Base:
    nonmissingtype

using
    DataStructures,
    BioSequences,
    FASTX,
    TranscodingStreams,
    CodecZlib,
    LightGraphs

import FASTX.FASTA: identifier, description, sequence, seqlen
export identifier, description, sequence, seqlen

include("io.jl")
include("family.jl")
include("domain.jl")
include("matches.jl")
include("models.jl")
include("linclust.jl")
include("mda.jl")
include("sequence.jl")
include("minhash.jl")
include("sketch/graph.jl")

end # module
