
const STATLIST = [:n, :posn, :mean, :sd, :se, :median, :min, :max, :geom]

function sortbyvec!(a, vec)
    sort!(a, by = x -> findfirst(y -> x == y, vec))
end

ispositive(::Missing) = false
ispositive(x::AbstractFloat) = isnan(x) ? false : x > zero(x)
ispositive(x) = x > zero(x)
################################################################################
struct SkipNonPositive{T}
    x::T
end
skipnonpositive(itr) = SkipNonPositive(itr)

Base.IteratorEltype(::Type{SkipNonPositive{T}}) where {T} = Base.IteratorEltype(T)
Base.eltype(::Type{SkipNonPositive{T}}) where {T} = nonmissingtype(eltype(T))
function Base.iterate(itr::SkipNonPositive, state...)
    y = iterate(itr.x, state...)
    y === nothing && return nothing
    item, state = y
    while !ispositive(item)
        y = iterate(itr.x, state)
        y === nothing && return nothing
        item, state = y
    end
    item, state
end
Base.IndexStyle(::Type{<:SkipNonPositive{T}}) where {T} = IndexStyle(T)
Base.eachindex(itr::SkipNonPositive) =
    Iterators.filter(i -> ispositive(@inbounds(itr.x[i])), eachindex(itr.x))
Base.keys(itr::SkipNonPositive) =
    Iterators.filter(i -> ispositive(@inbounds(itr.x[i])), keys(itr.x))
Base.@propagate_inbounds function getindex(itr::SkipNonPositive, I...)
    v = itr.x[I...]
    !ispositive(v) && throw(ErrorException("the value at index $I is non positive"))
    v
end
function Base.length(itr::SkipNonPositive)
    n = 0
    for i in itr n+=1 end
    n
end
length2(x) = length(x)
function length2(itr::Base.SkipMissing)
    n = 0
    for i in itr n+=1 end
    n
end
################################################################################

function dataimport(data; vars, sort)
    dataimport_(data, vars, sort)
end

function dataimport_(data, vars, sort)
    cols   = Tables.columns(data)
    cdata  = Tuple(Tables.getcolumn(cols, y) for y in sort)
    d      = Dict{Tuple{eltype.(cdata)...}, Vector{Int}}()
    indsdict!(d, cdata)
    varlabel = :Variable
    while varlabel in sort
        varlabel = Symbol(varlabel, '_')
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

function statfunc(vec, logvec, s)
    #f in STATLIST ? eval(:($f(v))) : error("Unknown function")
end

function descriptives(data; kwargs...)
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:corrected in k)
        kwargs[:corrected] = true
    end
    if !(:skipmissing in k)
        kwargs[:skipmissing] = false
    end
    if !(:skipnonpositive in k)
        kwargs[:skipnonpositive] = false
    end

    if !(:stats in k)
        kwargs[:stats] = [:n, :mean, :sd, :se, :median, :min, :max]
    end

     kwargs[:stats] ⊆ STATLIST || error("Some statistics not known!")

    ds = Vector{Descriptives}(undef, length(data))
    i  = 1
    for d in data
        # skipmissing
        if kwargs[:skipmissing]
            logvec = vec = skipmissing(d.obs)
        else
            vec = d.obs
        end
        n_  = logn_ = length2(vec)
        # skipnonpositive
        if kwargs[:skipnonpositive]
            logvec = skipnonpositive(d.obs)
            logn_  = length2(logvec)
        elseif !kwargs[:skipmissing]
            logvec = d.obs
        end

        result = Dict{Symbol, Float64}()

        for s in kwargs[:stats]
            if s == :n
                result[s] = n_
            elseif s == :posn
                result[s] = logn_
            elseif s == :mean
                result[s] = sum(vec) / n_
            elseif s == :sd
                Base.ht_keyindex(result, :mean) > 0 || (result[:mean] = sum(vec) / n_)
                result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean])
            elseif s == :var
                Base.ht_keyindex(result, :mean) > 0 || (result[:mean] = sum(vec) / n_)
                result[:var] = var(vec; corrected = kwargs[:corrected], mean = result[:mean])
            elseif s == :se
                Base.ht_keyindex(result, :mean) > 0 || (result[:mean] = sum(vec) / n_)
                Base.ht_keyindex(result, :sd) > 0 || (result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]))
                result[:se] = result[:sd] / sqrt(n_)
            elseif s == :median
                result[:median] = median(vec)
            elseif s == :min
                result[:min] = minimum(vec)
            elseif s == :max
                result[:max] = maximum(vec)
            elseif s == :geom
                result[:geom] = sum(log, logvec) / logn_
            end
        end
        filter!(x -> x.first in kwargs[:stats], result)
        ds[i] =  Descriptives(d, result)
        i += 1
    end
    DataSet(identity.(ds))
end




function Base.show(io::IO, obj::DataSet{DS}) where DS <: Descriptives
    idset  = Set(keys(first(obj).data.id))
    resset = Set(keys(first(obj).result))
    if length(obj) > 1
        for i = 2:length(obj)
            union!(idset,  Set(keys(obj[i].data.id)))
            union!(resset, Set(keys(obj[i].result)))
        end
    end
    mt1 = metida_table((getid(obj, :, c) for c in idset)...; names = idset)
    ressetl = sortbyvec!(collect(resset), STATLIST)
    mt2 = metida_table((obj[:, c] for c in ressetl)...; names = ressetl)
    show(io, MetidaTable(merge(mt1.table, mt2.table)))
end
