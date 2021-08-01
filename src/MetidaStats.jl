module MetidaStats

    using Statistics


    import MetidaBase: AbstractIdData, AbstractIDResult, DataSet, Tables, indsdict!, metida_table, getid, MetidaTable

    import Base: ht_keyindex

    export descriptives, dataimport

    include("types.jl")
    include("descriptive.jl")

end
