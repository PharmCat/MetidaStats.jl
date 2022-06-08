module MetidaStats

    using MetidaBase

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

    export descriptives, dataimport, cvfromsd

    include("types.jl")
    include("descriptive.jl")

end
