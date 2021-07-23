


struct ObsData{T, K, V} <: AbstractIdData
    obs::T
    var::Symbol
    id::Dict{K, V}
end


struct Descriptives{T <: ObsData} <: AbstractIDResult{T}
    data::T
    result::Dict
end

#=
function indsdict!(d::Dict{T}, cdata::Tuple) where T
    @inbounds for (i, element) in enumerate(zip(cdata...))
        ind = ht_keyindex(d, element)
        if ind > 0
            push!(d.vals[ind], i)
        else
            v = Vector{Int}(undef, 1)
            v[1] = i
            d[element] = v
        end
    end
    d
end
=#
