module MetidaStats

    import MetidaBase: AbstractIdData, AbstractIDResult, DataSet, Tables, indsdict!, metida_table, getid, MetidaTable

    import Base: ht_keyindex

    export descriptives, dataimport

    using Statistics

    include("types.jl")
    include("descriptive.jl")

end
