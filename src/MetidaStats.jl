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



    using MetidaBase, Statistics

    import MetidaBase: sortbyvec!

    import Base: ht_keyindex, getindex

    export descriptives, dataimport

    include("types.jl")
    include("descriptive.jl")

end
