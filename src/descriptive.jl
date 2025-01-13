
const STATLIST = [:n, 
:posn, 
:mean, 
:var, 
:bvar, 
:geom, 
:logmean, 
:logvar, 
:sd, 
:se, 
:cv, 
:geocv, 
:lci, 
:uci, 
:lmeanci, 
:umeanci, 
:median, 
:min, 
:max, 
:range, 
:q1, 
:q3, 
:iqr, 
:kurt,
:skew, 
:harmmean, 
:ses, 
:sek, 
:sum]

struct DStat
    sym::Symbol
    name
    corrected
    level
end

DStat(s::Symbol; corrected = nothing, level = nothing)              = DStat(s, string(s), corrected, level)
DStat(p::Pair{Symbol, <:Any}; corrected = nothing, level = nothing) = DStat(p[1], p[2], corrected, level)
DStat(s::String; corrected = nothing, level = nothing)              = DStat(Symbol(s), s, corrected, level)
DStat(p::Pair{String, <:Any}; corrected = nothing, level = nothing) = DStat(Symbol(p[1]), p[2], corrected, level)
DStat(ds::DStat) = d
make_stats(v::AbstractVector) = map(DStat, v)
make_stats(x) = [DStat(x)]
statssymbs(v::AbstractVector) = map(x->x.sym, v)
correctedval(a::Nothing, b) = b
correctedval(a::Bool, b) = a
qval(level::Nothing, q, n_) = q
qval(level::AbstractFloat, q, n_) = begin
    if n_ > 1 && a >= 0 && a <= 1 qv = quantile(TDist(n_ - 1), 1 - (1 - level) / 2) else qv = NaN end
    qv
end
levelval(level::Nothing, l) = l
levelval(level::AbstractFloat, l) = level 
"""
    dataimport(data; vars, sort = nothing)

Import data.

* data - tabular data;
* vars - variables;
* sort - sort by categories.
"""
function dataimport(data; vars, sort = nothing)
    if isa(vars, Symbol) vars = [vars] end
    if eltype(vars) <: Integer vars = Tables.columnnames(data)[vars] end
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
function dataimport_(data, vars, ::Nothing)
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
- `stats` - default set `stats = [:n, :mean, :sd, :se, :median, :min, :max]`;
- `corrected` - use corrected var (true);
- `level` - level for confidence intervals (0.95);

Possible values for `stats` is: 

* :n - number of observbations;
* :posn - positive (non-negative) number of observations;
* :mean - arithmetic mean;
* :var - variance;
* :bvar - variance with no correction;
* :geom - geometric mean; 
* :logmean - arithmetic mean for log-transformed data;
* :logvar - variance for log-transformed data ``σ^2_{log}``;
* :sd - standard deviation (or σ);
* :se - standard error; 
* :cv - coefficient of variation; 
* :geocv - coefficient of variation for log-transformed data (``CV = sqrt{exp(σ^2_{log})-1}``);
* :lci - lower confidence interval;
* :uci - upper confidence interval; 
* :lmeanci - lower confidence interval for mean; 
* :umeanci - lower confidence interval for mean; 
* :median - median,;
* :min - minimum; 
* :max - maximum; 
* :range - range; 
* :q1 - lower quartile;
* :q3 - upper quartile;
* :iqr - inter quartile range; 
* :kurt - kurtosis;
* :skew - skewness; 
* :harmmean - harmonic mean; 
* :ses standard error of skewness; 
* :sek - standard error of kurtosis; 
* :sum - sum.

"""
function descriptives(data, vars, sort = nothing; kwargs...)
    if isa(vars, String) vars = [Symbol(vars)] end
    if isa(vars, Symbol) vars = [vars] end
    if isa(sort, String) sort = [Symbol(sort)] end
    if isa(sort, Symbol) sort = [sort] end
    if eltype(vars) <: Integer vars = Tables.columnnames(data)[vars] end
    if !isnothing(sort)
        vars = setdiff(vars, sort)
        if length(sort) == 0 sort = nothing end
    end
    descriptives(dataimport_(data, vars, sort); kwargs...)
end
"""
    descriptives(data; vars = nothing, sort = nothing, kwargs...)

If `vars` is nothing - try to include all collumns with numbers.
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
    if eltype(vars) <: Integer vars = Tables.columnnames(data)[vars] end
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
    if :qts in k
        # Check quantiles
    end


    if !(:stats in k)
        stats          = [:n, :mean, :sd, :se, :median, :min, :max]
        kwargs[:stats] = make_stats(stats)
    else
        kwargs[:stats] = make_stats(kwargs[:stats])
        stats          = statssymbs(kwargs[:stats])
    end

    stats ⊆ STATLIST || error("Some statistics not known!")

    if any(x -> x in [:geom, :geomean, :logmean, :logvar, :geocv], stats) 
        makelogvec = true 
    else 
        makelogvec = false 
    end

    if any(x -> x in [:lci, :uci, :lmeanci, :umeanci], stats)
        cicalk = true
    else
        cicalk = false
    end
    #
    ds = Vector{Descriptives}(undef, length(data))
    i  = 1
    for d in data
        ds[i] =  Descriptives(d, kwargs[:stats], descriptives_(d.obs, kwargs, makelogvec, cicalk))
        i += 1
    end
    DataSet(identity.(ds))
    
end

function descriptives_(obsvec, kwargs, logstats, cicalk)
            # skipmissing
            if kwargs[:skipmissing]
                vec = skipnanormissing(obsvec)
            else
                vec = obsvec
            end
            n_ = length(vec)
            if cicalk
                if n_ > 1 && kwargs[:level] >= 0 && kwargs[:level] <= 1 q = quantile(TDist(n_ - 1), 1 - (1 - kwargs[:level]) / 2) else q = NaN end# add tdist / normal option # add multiple CI ?
            end
            # skipnonpositive
            # logstats = makelogvec #calk logstats
            if logstats
                if kwargs[:skipnonpositive]
                    logvec = log.(skipnonpositive(obsvec))
                else
                    if kwargs[:skipmissing]
                        itr = skipnanormissing(obsvec)
                        if any(!ispositive, itr)
                            logvec = Float64[]
                            logstats = false
                        else
                            logvec = log.(itr)
                        end
                    else
                        if any(!ispositive, obsvec)
                            logvec = Float64[]
                            logstats = false
                        else
                            logvec = log.(obsvec)
                        end
                    end
                end
            end
            logn_  = length(skipnonpositive(obsvec))
    
            result = OrderedDict{Symbol, Float64}()
    
            for s in kwargs[:stats]
    
                if s.sym == :n
                    result[s.sym] = n_
                elseif s.sym == :posn
                    result[s.sym] = logn_
                elseif !(n_ > 0)
                    result[s.sym] = NaN
                    continue
                elseif s.sym == :mean
                    result[s.sym] = sum(vec) / n_
                elseif s.sym == :sd
                    haskey(result, :mean) || begin result[:mean] =  sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    result[s.sym] = std(vec; corrected = corrected, mean = result[:mean])
                elseif s.sym == :var
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    result[s.sym] = var(vec; corrected = corrected, mean = result[:mean])
                elseif s.sym == :bvar
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    result[s.sym] = var(vec; corrected = false, mean = result[:mean])
                elseif s.sym == :se
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    result[s.sym] = result[:sd] / sqrt(n_)
                elseif s.sym == :cv
                    haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    result[s.sym] = abs(result[:sd] / result[:mean] * 100)
                elseif s.sym == :uci
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    qv = qval(s.level, q, n_) 
                    result[Symbol(string(s.sym)*"_$(levelval(s.level, kwargs[:level]))")] = result[:mean] + qv * result[:sd]
                elseif s.sym == :lci
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    qv = qval(s.level, q, n_) 
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    result[Symbol(string(s.sym)*"_$(levelval(s.level, kwargs[:level]))")] = result[:mean] - qv * result[:sd]
                elseif s.sym == :umeanci
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    haskey(result, :se) || begin result[:se] = result[:sd] / sqrt(n_) end
                    qv = qval(s.level, q, n_) 
                    result[Symbol(string(s.sym)*"_$(levelval(s.level, kwargs[:level]))")] = result[:mean] + qv * result[:se]
                elseif s.sym == :lmeanci
                    haskey(result, :mean) || begin result[:mean] = sum(vec) / n_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :sd) || begin result[:sd] = std(vec; corrected = corrected, mean = result[:mean]) end
                    haskey(result, :se) || begin result[:se] = result[:sd] / sqrt(n_) end
                    qv = qval(s.level, q, n_) 
                    result[Symbol(string(s.sym)*"_$(levelval(s.level, kwargs[:level]))")] = result[:mean] - qv * result[:se]
                elseif s.sym == :median
                    result[s.sym] = median(vec)
                elseif s.sym == :min
                    result[s.sym] = minimum(vec)
                elseif s.sym == :max
                    result[s.sym] = maximum(vec)
                elseif s.sym == :q1
                    result[s.sym] = quantile(vec, 0.25)
                elseif s.sym == :q3
                    result[s.sym] = quantile(vec, 0.75)
                elseif s.sym == :iqr
                    result[s.sym] = abs(quantile(vec, 0.75) - quantile(vec, 0.25))
                elseif s.sym == :range
                    result[s.sym] = abs(maximum(vec) - minimum(vec))
                elseif s.sym == :kurt
                    haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                    result[s.sym] = kurtosis_(vec,  result[:mean])
                elseif s.sym == :skew
                    haskey(result, :mean)  || begin result[:mean] = sum(vec) / n_ end
                    result[s.sym] = skewness_(vec,  result[:mean])
                elseif s.sym == :harmmean
                    result[s.sym] = harmmean(vec)
                elseif s.sym == :ses
                    result[s.sym] = sesvec(vec)
                elseif s.sym == :sek
                    result[s.sym] = sekvec(vec)
                elseif s.sym == :sum
                    result[s.sym] = sum(vec)
                elseif !logstats
                    result[s.sym] = NaN
                    continue
                elseif s.sym == :logmean
                    result[s.sym] = sum(logvec) / logn_
                elseif s.sym == :geom || s.sym == :geomean
                    haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                    result[s.sym] = exp(result[:logmean])
                elseif s.sym == :logvar
                    haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    result[s.sym] = var(logvec; corrected = corrected, mean = result[:logmean])
                elseif s.sym == :geocv
                    haskey(result, :logmean) || begin result[:logmean] = sum(logvec) / logn_ end
                    corrected = correctedval(s.corrected, kwargs[:corrected])
                    haskey(result, :logvar) || begin result[:logvar] = var(logvec; corrected = corrected, mean = result[:logmean]) end
                    result[s.sym] = sqrt(exp(result[:logvar]) - 1)*100
                end
            end
            filter!(x -> x.first in statssymbs(kwargs[:stats]) || occursin("ci", string(x.first)), result)

            #if any(:lci in keys(result)) Symbol(string(s)*@sprintf("%g", kwargs[:level]*100))

            result
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

function MetidaBase.metida_table_(obj::DataSet{DS}; sort = nothing, stats = nothing, id = nothing) where DS <: Descriptives
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
        if isnothing(sort)
            ressetl = sortbyvec!(collect(intersect(resset, stats)), collect(keys(first(obj).result)))
        else
            ressetl = sortbyvec!(collect(intersect(resset, stats)), sort)
        end
    else
        if isnothing(sort)
            ressetl = sortbyvec!(collect(resset), collect(keys(first(obj).result)))
        else
            ressetl = sortbyvec!(collect(resset), sort)
        end
    end
    if !isnothing(id)
        if isa(id, Symbol) id = [id] end
        id ⊆ idset || error("Some id not in dataset!")
        idset = intersect(idset, id)
    end
    mt1 = metida_table((getid(obj, :, c) for c in idset)...; names = idset)
    mt2 = metida_table((obj[:, c] for c in ressetl)...; names = ressetl)
    merge(mt1.table, mt2.table)
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
