module MetidaStats

    import MetidaBase

    import MetidaBase: AbstractIdData, AbstractIDResult, DataSet, Tables, indsdict!, metida_table, getid, MetidaTable, isnanormissing
    import Base: ht_keyindex

    using Statistics

    export descriptives, dataimport

    include("types.jl")
    include("descriptive.jl")

end
