using MetidaStats
using MetidaBase
using Test
using DataFrames, CSV

path     = dirname(@__FILE__)
io       = IOBuffer();


@testset "  Descripitive statistics                                  " begin
    ds  = CSV.File(path*"/csv/ds.csv") |> DataFrame
    di  = MetidaStats.dataimport(ds, vars = [:var1, :var2], sort = [:col, :row])
    des = MetidaStats.descriptives(di; skipmissing = true, skipnonpositive = true, stats = MetidaStats.STATLIST)

    di  = MetidaStats.dataimport(ds, vars = [:var1, :var2])


    mt  = MetidaBase.metida_table(des; stats = [:mean, :geom])

    mt  = MetidaBase.metida_table(des; stats = [:mean, :geom], id = [:Variable,:row])

    des2 = MetidaStats.descriptives(ds, [:var1, :var2], [:col, :row]; skipmissing = true, skipnonpositive = true, stats = MetidaStats.STATLIST)

    @test des[:, :mean] == des2[:, :mean]

    show(io, des)
end
