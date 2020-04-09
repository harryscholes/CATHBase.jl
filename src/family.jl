abstract type AbstractFamily end

identifier(x::AbstractFamily) = x.identifier

const SUPER_FAMILY_REGEX = r"^\d+.\d+.\d+.\d+$"

struct SuperFamily <: AbstractFamily
    identifier::String

    SuperFamily(id::AbstractString) = occursin(SUPER_FAMILY_REGEX, id) ? new(id) :
        throw(ArgumentError("Invalid SuperFamily id"))
end

Base.hash(x::SuperFamily, h::UInt) = hash(:SuperFamily, hash(identifier(x), h))
Base.:(==)(x::SuperFamily, y::SuperFamily) = hash(x) == hash(y)

macro sf_str(id)
    SuperFamily(id)
end

Base.string(xs::AbstractVector{SuperFamily}) = join(identifier.(xs), "-")

struct FunctionalFamily <: AbstractFamily
    identifier::String
end

Base.hash(x::FunctionalFamily, h::UInt) = hash(:FunctionalFamily, hash(identifier(x), h))
Base.:(==)(x::FunctionalFamily, y::FunctionalFamily) = hash(x) == hash(y)

macro ff_str(id)
    FunctionalFamily(id)
end
