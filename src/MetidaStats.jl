module MetidaStats

    using MetidaBase, Distributions, Printf, DataStructures

    import MetidaBase: AbstractIdData,
    AbstractIDResult,
    DataSet,
    Tables,
    indsdict!,
    metida_table,
    getid,
    MetidaTable,
    sortbyvec!,
    skipnanormissing, isnanormissing,
    sortbyvec!,
    skipnonpositive, ispositive,
    cvfromsd


    using StatsBase, Statistics


    import Base: getindex, names, length

    export descriptives, dataimport, cvfromsd, metida_table

    include("types.jl")
    include("descriptive.jl")

end
