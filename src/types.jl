


struct ObsData{T, K, V} <: AbstractIdData
    obs::T
    var::Symbol
    id::Dict{K, V}
end


struct Descriptives{T <: ObsData} <: AbstractIDResult{T}
    data::T
    result::AbstractDict
end

################################################################################
