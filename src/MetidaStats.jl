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
    isnanormissing


    import Base: ht_keyindex, getindex

    using MetidaBase, Statistics

    export descriptives, dataimport

    include("types.jl")
    include("descriptive.jl")

end
