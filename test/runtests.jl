using MetidaStats
using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();




@testset "  Descripitive statistics                                  " begin
    ds  = CSV.File(path*"/csv/ds.csv") |> DataFrame
    di = MetidaStats.dataimport_(ds, [:var1, :var2], [:col, :row])
    des= MetidaStats.descriptives(di; skipmissing = true)
end
