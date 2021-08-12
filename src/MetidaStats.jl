module MetidaStats

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
    skipnonpositive, ispositive


    using MetidaBase, Statistics


    import Base: ht_keyindex, getindex, names

    export descriptives, dataimport, cvfromsd

    include("types.jl")
    include("descriptive.jl")

end
