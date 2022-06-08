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

    mt  = MetidaStats.metida_table(des; stats = [:mean, :geom])

    mt  = MetidaStats.metida_table(des; stats = [:mean, :geom], id = [:Variable,:row])

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


    @test des2[3, :geom]  ≈ 48.4385401710072
    @test des2[3, :geom]  ≈ geomean(MetidaStats.skipnonpositive(des2[3].data.obs))
    @test des2[1, :harmmean] ≈ 27.416118945035485
    @test des2[5, :harmmean] ≈ 38.7909224529887
    @test des2[7, :harmmean] ≈ 15.394283582287837

    @test des2[1, :min]  ≈ 8.374838503411963
    @test des2[1, :max]  ≈ 94.74467112693401
    @test des2[1, :range] ≈ 86.36983262352204

    @test des2[1, :geom]  ≈ 39.7551197516893
    @test des2[1, :logmean]  ≈ 3.682738631591126
    @test des2[1, :logvar]  ≈ 0.7054485985860965

    @test des2[1, :cv]  ≈ 59.176116266880264
    @test des2[1, :geocv]  ≈ 101.23017254525884

    @test des2[2, :geom]  ≈ geomean(MetidaStats.skipnonpositive(des2[2].data.obs))
    @test des2[2, :logmean]  ≈ 2.631096578600917
    @test des2[2, :logvar]  ≈ 0.7574778424131346
    @test des2[2, :cv]  ≈ 824.4053998875627
    @test des2[2, :geocv]  ≈ 106.437302966393

    @test des2[3, :cv]  ≈ 46.41856487180453
    @test des2[3, :cv]  ≈ des2[3, :sd]/des2[3, :mean]*100

    di  = MetidaStats.dataimport(ds, vars = [:var1, :var2], sort = [:col, :row])
    sort!(di, [:col, :row, :Variable])


    des2[1, :skew] ≈ skewness(di[1].obs)

    des2[1, :kurt] ≈ kurtosis(di[1].obs)

    #dfdes = sort!(DataFrame(MetidaStats.metida_table(des2)), [:col, :row])
    des3 = MetidaStats.descriptives(ds, [:var1, :var2], [:col, :row]; skipmissing = true, skipnonpositive = false, stats = MetidaStats.STATLIST)
    sort!(des3, [:col, :row, :Variable])
    @test des3[2, :mean] ≈ des2[2, :mean]
    @test des3[2, :geom] === NaN
    des4 = MetidaStats.descriptives(ds, [:var1], [:col, :row]; skipmissing = false, skipnonpositive = false, stats = MetidaStats.STATLIST)
    sort!(des4, [:col, :row, :Variable])
    @test des4[1, :mean] ≈ des2[1, :mean]
end
