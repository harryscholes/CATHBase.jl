struct DomainSection
    start::Int
    stop::Int
    function DomainSection(start, stop)
        start ≤ stop ? new(start, stop) : throw(ArgumentError("stop is not > start"))
    end
end

function DomainSection(s::AbstractString)
    if (m = match(r"\d+\-\d+", s)) isa RegexMatch
        return DomainSection(parse.(Int, split(m.match, '-'))...)
    else
        throw(ArgumentError("string cannot be parsed to a `DomainSection`"))
    end
end

start(ds::DomainSection) = ds.start
stop(ds::DomainSection) = ds.stop
Base.length(ds::DomainSection) = stop(ds) - start(ds) + 1
Base.string(ds::DomainSection) = string(start(ds))*"-"*string(stop(ds))
mass(ds::DomainSection) = sum(start(ds):stop(ds))
Base.hash(ds::DomainSection, h::UInt) = hash(:DomainSection, hash(start(ds), hash(stop(ds), h)))
Base.:(==)(x::DomainSection, y::DomainSection) = hash(x) == hash(y)

struct DomainBoundary
    boundary::Vector{DomainSection}
end

function DomainBoundary(s::AbstractString)
    if (m = match(r"^[\d\-\,]+$", s)) isa RegexMatch
        return DomainBoundary(DomainSection.(split(m.match, ',')))
    else
        throw(ArgumentError("string cannot be parsed to a `DomainBoundary`: "*s))
    end
end

DomainBoundary(start::Integer, stop::Integer) = DomainBoundary([DomainSection(start, stop)])

boundary(x::DomainBoundary) = x.boundary
iscontinuous(x::DomainBoundary) = length(boundary(x)) == 1
start(x::DomainBoundary) = iscontinuous(x) ? start(boundary(x)[1]) : error("Discontinuous domain")
stop(x::DomainBoundary) = iscontinuous(x) ? stop(boundary(x)[1]) : error("Discontinuous domain")
Base.length(x::DomainBoundary) = sum(length.(boundary(x)))
Base.getindex(seq::BioSequence, db::DomainBoundary) = *(map(ds->seq[start(ds):stop(ds)], boundary(db))...)
Base.string(b::DomainBoundary) = join(string.(boundary(b)), ',')
mass(db::DomainBoundary) = sum(mass.(boundary(db)))
centerofmass(db::DomainBoundary) = Int(floor(mass(db)/length(db)))
Base.hash(db::DomainBoundary, h::UInt) = hash(:DomainBoundary, hash(boundary(db), h))
Base.:(==)(x::DomainBoundary, y::DomainBoundary) = hash(x) == hash(y)

struct Domain
    superfamily::SuperFamily
    boundary::DomainBoundary
end

Domain(sf::AbstractString, b::AbstractString) = Domain(SuperFamily(sf), DomainBoundary(b))

superfamily(d::Domain) = d.superfamily
boundary(d::Domain) = d.boundary
Base.hash(d::Domain, h::UInt) = hash(:Domain, hash(superfamily(d), hash(boundary(d), h)))
Base.:(==)(x::Domain, y::Domain) = hash(x) == hash(y)

import Base: length, getindex

for f in (:iscontinuous, :start, :stop, :length, :getindex, :mass, :centerofmass)
    @eval ($f)(d::Domain) = ($f)(boundary(d))
end
