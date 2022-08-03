
const STATLIST = [:n, :posn, :mean, :var, :bvar, :geom, :logmean, :logvar, :sd, :se, :cv, :geocv, :lci, :uci, :lmeanci, :umeanci, :median, :min, :max, :range, :q1, :q3, :iqr, :kurt, :skew, :harmmean, :ses, :sek, :sum]

#=
function sortbyvec!(a, vec)
    sort!(a, by = x -> findfirst(y -> x == y, vec))
end
=#
#=
ispositive(::Missing) = false
ispositive(x::AbstractFloat) = isnan(x) ? false : x > zero(x)
ispositive(x) = x > zero(x)
################################################################################
struct SkipNonPositive{T}
    x::T
end
skipnonpositive(itr) = SkipNonPositive(itr)

Base.IteratorEltype(::Type{SkipNonPositive{T}}) where {T} = Base.IteratorEltype(T)
#Base.eltype(::Type{SkipNonPositive{T}}) where {T} = nonmissingtype(eltype(T))
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
################################################################################
struct SkipNaNorMissing{T}
    x::T
end
skipnanormissing(itr) = SkipNaNorMissing(itr)

Base.IteratorEltype(::Type{SkipNaNorMissing{T}}) where {T} = Base.IteratorEltype(T)
#Base.eltype(::Type{SkipNaNorMissing{T}}) where {T} = nonmissingtype(eltype(T))
function Base.iterate(itr::SkipNaNorMissing, state...)
    y = iterate(itr.x, state...)
    y === nothing && return nothing
    item, state = y
    while isnanormissing(item)
        y = iterate(itr.x, state)
        y === nothing && return nothing
        item, state = y
    end
    item, state
end
Base.IndexStyle(::Type{<:SkipNaNorMissing{T}}) where {T} = IndexStyle(T)
Base.eachindex(itr::SkipNaNorMissing) =
    Iterators.filter(i -> isnanormissing(@inbounds(itr.x[i])), eachindex(itr.x))
Base.keys(itr::SkipNaNorMissing) =
    Iterators.filter(i -> isnanormissing(@inbounds(itr.x[i])), keys(itr.x))
Base.@propagate_inbounds function getindex(itr::SkipNaNorMissing, I...)
    v = itr.x[I...]
    !isnanormissing(v) && throw(ErrorException("The value at index $I is NaN or missing!"))
    v
end
function Base.length(itr::SkipNaNorMissing)
    n = 0
    for i in itr n+=1 end
    n
end
=#
################################################################################
#=
length2(x) = length(x)
function length2(itr::Base.SkipMissing)
    n = 0
    for i in itr n+=1 end
    n
end
=#
################################################################################
"""
    dataimport(data; vars, sort = nothing)

Import data.

* data - tabular data;
* vars - variables;
* sort - sort by categories.
"""
function dataimport(data; vars, sort = nothing)
    if isa(vars, Symbol) vars = [vars] end

    if !isnothing(sort) && isa(sort, Symbol) sort = [sort] end

    dataimport_(data, vars, sort)
end

function dataimport_(data, vars::AbstractVector, sort::AbstractVector)
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
function dataimport_(data, vars, sort::Nothing)
    sdata = Vector{ObsData}(undef, length(vars))
    for i in 1:length(vars)
        sdata[i] = ObsData(Tables.getcolumn(data, vars[i]), vars[i],  Dict(:Variable=>vars[i]))
    end
    DataSet(identity.(sdata))
end
"""
    descriptives(data, vars, sort = nothing; kwargs...)

* kwargs:
- `skipmissing` - drop NaN and Missing values, default = true;
- `skipnonpositive` - drop non-positive values (and NaN, Missing) for "log-statistics" - :geom, :geomean, :logmean, :logvar, :geocv;
- `stats` - default set `stats = [:n, :mean, :sd, :se, :median, :min, :max]`
Possible values for `stats` is: :n, :posn, :mean, :var, :geom, :logmean, :logvar, :sd, :se, :cv, :geocv, :median, :min, :max, :range, :q1, :q3, :iqr, :kurt, :skew, :harmmean, :ses, :sek, :sum

"""
function descriptives(data, vars = nothing, sort = nothing; kwargs...)
    if isa(vars, String) vars = [Symbol(vars)] end
    if isa(vars, Symbol) vars = [vars] end
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end
    if !isnothing(sort)
        vars = setdiff(vars, sort)
    end
    descriptives(dataimport_(data, vars, sort); kwargs...)
end
"""
    descriptives(data; vars = nothing, sort = nothing, kwargs...)

Try to include all collumns with numbers.
"""
function descriptives(data; vars = nothing, sort = nothing, kwargs...)
    if isnothing(vars)
        vars = Vector{Symbol}(undef, 0)
        coln   = Tables.columnnames(data)
        for i in coln
            if eltype(Tables.getcolumn(data, i)) <: Number push!(vars, i) end
        end
        if length(vars) == 0 error("No column found for descriptive statistics!") end
    end
    descriptives(data, vars, sort; kwargs...)
end
"""
    descriptives(data::DataSet{T}; kwargs...) where T <: ObsData

Descriptive statistics for `dataimport` structure, see  [`dataimport`](@ref).

"""
function descriptives(data::DataSet{T}; kwargs...) where T <: ObsData
    kwargs = Dict{Symbol, Any}(kwargs)
    k = keys(kwargs)

    if !(:corrected in k)
        kwargs[:corrected] = true
    end
    if !(:skipmissing in k)
        kwargs[:skipmissing] = true
    end
    if !(:skipnonpositive in k)
        kwargs[:skipnonpositive] = false
    end
    if !(:level in k)
        kwargs[:level] = 0.95
    end


    if !(:stats in k)
        kwargs[:stats] = [:n, :mean, :sd, :se, :median, :min, :max]
    else
        if isa(kwargs[:stats], Symbol)  kwargs[:stats] = [kwargs[:stats]] end
        if isa(kwargs[:stats], String)  kwargs[:stats] = [Symbol(kwargs[:stats])] end
    end

    kwargs[:stats] ⊆ STATLIST || error("Some statistics not known!")

    if any(x -> x in [:geom, :geomean, :logmean, :logvar, :geocv], kwargs[:stats]) makelogvec = true else makelogvec = true end

    if any(x -> x in [:lci, :uci, :lmeanci, :umeanci], kwargs[:stats])
        cicalk = true
    else
        cicalk = false
        q = NaN
    end

    ds = Vector{Descriptives}(undef, length(data))
    i  = 1
    for d in data
        # skipmissing
        if kwargs[:skipmissing]
            vec = skipnanormissing(d.obs)
        else
            vec = d.obs
        end
        n_ = length(vec)
        if cicalk
            if n_ > 1 q = quantile(TDist(n_ - 1), 1 - (1-kwargs[:level])/2) end
        end
        # skipnonpositive
        logstats = makelogvec #calk logstats
        if makelogvec
            if kwargs[:skipnonpositive]
                logvec = log.(skipnonpositive(d.obs))
            else
                if kwargs[:skipmissing]
                    itr = skipnanormissing(d.obs)
                    if any(!ispositive, itr)
                        logvec = Float64[]
                        logstats = false
                    else
                        logvec = log.(itr)
                    end
                else
                    if any(!ispositive, d.obs)
                        logvec = Float64[]
                        logstats = false
                    else
                        logvec = log.(d.obs)
                    end
                end
            end
        end
        logn_  = length(skipnonpositive(d.obs))

        result = Dict{Symbol, Float64}()

        for s in kwargs[:stats]

            if s == :n
                result[s] = n_
            elseif s == :posn
                result[s] = logn_
            elseif !(n_ > 0)
                result[s] = NaN
                continue
            elseif s == :mean
                result[s] = sum(vec) / n_
            elseif s == :sd
                haskey(result, :mean) || begin result[:mean] =  sum(vec) / n_ end
                result[s] = std(vec; corrected = kwargs[:corrected], mean = result[:mean])
            elseif s == :var
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                result[s] = var(vec; corrected = kwargs[:corrected], mean = result[:mean])
            elseif s == :bvar
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                result[s] = var(vec; corrected = false, mean = result[:mean])
            elseif s == :se
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                result[s] = result[:sd] / sqrt(n_)
            elseif s == :cv
                haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                result[s] = abs(result[:sd] / result[:mean] * 100)
            elseif s == :uci
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                result[s] = result[:mean] + q*result[:sd]
            elseif s == :lci
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                result[s] = result[:mean] - q*result[:sd]
            elseif s == :umeanci
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                haskey(result, :se) || begin result[:se] = result[:sd] / sqrt(n_) end
                result[s] = result[:mean] + q*result[:se]
            elseif s == :lmeanci
                haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                haskey(result, :sd) || begin result[:sd] = std(vec; corrected = kwargs[:corrected], mean = result[:mean]) end
                haskey(result, :se) || begin result[:se] = result[:sd] / sqrt(n_) end
                result[s] = result[:mean] - q*result[:se]
            elseif s == :median
                result[s] = median(vec)
            elseif s == :min
                result[s] = minimum(vec)
            elseif s == :max
                result[s] = maximum(vec)
            elseif s == :q1
                result[s] = quantile(vec, 0.25)
            elseif s == :q3
                result[s] = quantile(vec, 0.75)
            elseif s == :iqr
                result[s] = abs(quantile(vec, 0.75) - quantile(vec, 0.25))
            elseif s == :range
                result[s] = abs(maximum(vec) - minimum(vec))
            elseif s == :kurt
                haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                result[s] = kurtosis_(vec,  result[:mean])
            elseif s == :skew
                haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                result[s] = skewness_(vec,  result[:mean])
            elseif s == :harmmean
                result[s] = harmmean(vec)
            elseif s == :ses
                result[s] = sesvec(vec)
            elseif s == :sek
                result[s] = sekvec(vec)
            elseif s == :sum
                result[s] = sum(vec)
            elseif !logstats
                result[s] = NaN
                continue
            elseif s == :logmean
                result[s] = sum(logvec) / logn_
            elseif s == :geom || s == :geomean
                haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                result[s] = exp(result[:logmean])
            elseif s == :logvar
                haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                result[s] = var(logvec; corrected = kwargs[:corrected], mean = result[:logmean])
            elseif s == :geocv
                haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                haskey(result, :logvar) || begin result[:logvar] = var(logvec; corrected = kwargs[:corrected], mean = result[:logmean]) end
                result[s] = sqrt(exp(result[:logvar])-1)*100
            end
        end
        filter!(x -> x.first in kwargs[:stats], result)
        ds[i] =  Descriptives(d, result)
        i += 1
    end
    DataSet(identity.(ds))
end


function sesvec(vec)
    ses(length(vec))
end
function ses(n)
    return sqrt(6 * n *(n - 1) / ((n - 2) * (n + 1) * (n + 3)))
end

function sekvec(vec)
    n   = length(vec)
    sek(n; ses = ses(n))
end
function sek(n; ses = ses(n))
    return 2 * ses * sqrt((n * n - 1)/((n - 3) * (n + 5)))
end


function kurtosis_(v, m)
    n = length(v)
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment
    for i in v
        z = i - m
        z2 = z * z
        cm2 += z2
        cm4 += z2 * z2
    end
    cm4 /= n
    cm2 /= n
    return (cm4 / (cm2 * cm2)) - 3.0
end

function skewness_(v, m)
    n = length(v)
    cm2 = 0.0   # empirical 2nd centered moment (variance)
    cm3 = 0.0   # empirical 3rd centered moment
    for i in v
        z = i - m
        z2 = z * z
        cm2 += z2
        cm3 += z2 * z
    end
    cm3 /= n
    cm2 /= n
    return cm3 / sqrt(cm2 * cm2 * cm2)  # this is much faster than cm2^1.5
end

################################################################################
#
################################################################################

function MetidaBase.metida_table(obj::DataSet{DS}; sort = STATLIST, stats = nothing, id = nothing) where DS <: Descriptives
    idset  = Set(keys(first(obj).data.id))
    resset = Set(keys(first(obj).result))
    if length(obj) > 1
        for i = 2:length(obj)
            union!(idset,  Set(keys(obj[i].data.id)))
            union!(resset, Set(keys(obj[i].result)))
        end
    end
    if !isnothing(stats)
        stats ⊆ STATLIST || error("Some statistics not known!")
        if isa(stats, Symbol) stats = [stats] end
        ressetl = sortbyvec!(collect(intersect(resset, stats)), sort)
    else
        ressetl = sortbyvec!(collect(resset), sort)
    end
    if !isnothing(id)
        if isa(id, Symbol) id = [id] end
        id ⊆ idset || error("Some id not in dataset!")
        idset = intersect(idset, id)
    end
    mt1 = metida_table((getid(obj, :, c) for c in idset)...; names = idset)
    mt2 = metida_table((obj[:, c] for c in ressetl)...; names = ressetl)
    MetidaTable(merge(mt1.table, mt2.table))
end

################################################################################
# SHOW
################################################################################

function Base.show(io::IO, obj::DataSet{DS}) where DS <: Descriptives
    show(io, metida_table(obj))
end

function Base.show(io::IO, obj::DataSet{OD}) where OD <: ObsData
    println(io, "DataSet: observations")
    for i in obj
        println(io, "  Var: $(i.var); ID: $(i.id); Length: $(length(i.obs))")
    end
end


################################################################################
# To MetidaBase
#=
function cvfromsd(σ)
    return sqrt(exp(σ^2)-1)
end
=#
