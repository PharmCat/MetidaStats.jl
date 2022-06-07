using MetidaStats
using MetidaBase
using StatsBase
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

    sort!(des2, [:col, :row, :Variable])

    @test des2[:, :n]  == [24
    22
    44
    44
    20
    20
    83
    83]

    @test des2[:, :mean]  ≈ [51.84339580144830000
    -3.2427526832130800
    58.06263961606340000
    1.8800364385656400
    51.84107641126520000
    0.4353627340913350
    47.25778648204360000
    -3.2515999935069800]

    @test des2[:, :sd]  ≈ [30.67890817616400000
    26.7334282254075000
    26.95184403646440000
    28.9591749375126000
    25.29979111399270000
    27.5437061096263000
    27.16599075359280000
    28.8185876098308000]

    @test des2[:, :se]  ≈ [6.26230590810839000
    5.6995860482981700
    4.06314336714229000
    4.3657598866360500
    5.65720527474328000
    6.1589599213400700
    2.98185487195489000
    3.1632509429411400]

    @test des2[1, :geom]  ≈ 39.7551197516893
    @test des2[3, :geom]  ≈ 48.4385401710072
    @test des2[5, :harmmean] ≈ 38.7909224529887
    @test des2[7, :harmmean] ≈ 15.394283582287837

    @test des2[1, :logvar]  ≈ 0.5413586421629522
    @test des2[1, :cv]  ≈ 46.41856487180453
    @test des2[1, :geocv]  ≈ 84.75493413206702

    di  = MetidaStats.dataimport(ds, vars = [:var1, :var2], sort = [:col, :row])
    sort!(di, [:col, :row, :Variable])


    des2[1, :skew] ≈ skewness(di[1].obs)


    des2[1, :kurt] ≈ kurtosis(di[1].obs)

end
