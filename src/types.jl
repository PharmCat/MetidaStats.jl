


struct ObsData{T, K, V} <: AbstractIdData
    obs::T
    var::Symbol
    id::Dict{K, V}
end


struct Descriptives{T <: ObsData} <: AbstractIDResult{T}
    data::T
    result::Dict{:Symbol, :Float64}
end

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

function descriptive_(data, vars, sort)
    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)
    varlabel = :Variable
    while varlabel in sort
        varlabel = string(varlabel) * "_"
    end

    sdata = Vector{ObsData}(undef, length(d)*length(vars))

    i = one(Int)
    @inbounds for (k, v) in d
        for c in vars
            dict = Dict{Symbol, Any}(sort .=> k)
            dict[varlabel] = c
            sdata[i] = ObsData(view(Tables.getcolumn(data, c), v), varlabel,  dict)
            i += one(Int)
        end
    end
    return DataSet(identity.(sdata))
end
