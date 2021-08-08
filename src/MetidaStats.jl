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
    isnanormissing,
    sortbyvec!


    using MetidaBase, Statistics


    import Base: ht_keyindex, getindex

    export descriptives, dataimport, cvfromsd

    include("types.jl")
    include("descriptive.jl")

end
