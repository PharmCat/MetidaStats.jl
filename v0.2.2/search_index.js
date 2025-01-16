var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#dataimport","page":"API","title":"dataimport","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaStats.dataimport","category":"page"},{"location":"api/#MetidaStats.dataimport","page":"API","title":"MetidaStats.dataimport","text":"dataimport(data; vars, sort = nothing)\n\nImport data.\n\ndata - tabular data;\nvars - variables;\nsort - sort by categories.\n\n\n\n\n\n","category":"function"},{"location":"api/#descriptives","page":"API","title":"descriptives","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MetidaStats.descriptives","category":"page"},{"location":"api/#MetidaStats.descriptives","page":"API","title":"MetidaStats.descriptives","text":"descriptives(data, vars, sort = nothing; kwargs...)\n\nkwargs:\nskipmissing - drop NaN and Missing values, default = true;\nskipnonpositive - drop non-positive values (and NaN, Missing) for \"log-statistics\" - :geom, :geomean, :logmean, :logvar, :geocv;\nstats - default set stats = [:n, :mean, :sd, :se, :median, :min, :max];\ncorrected - use corrected var (true);\nlevel - level for confidence intervals (0.95);\n\nPossible values for stats is: \n\n:n - number of observbations;\n:posn - positive (non-negative) number of observations;\n:mean - arithmetic mean;\n:var - variance;\n:bvar - variance with no correction;\n:geom - geometric mean; \n:logmean - arithmetic mean for log-transformed data;\n:logvar - variance for log-transformed data σ^2_log;\n:sd - standard deviation (or σ);\n:se - standard error; \n:cv - coefficient of variation; \n:geocv - coefficient of variation for log-transformed data (CV = sqrtexp(σ^2_log)-1);\n:lci - lower confidence interval;\n:uci - upper confidence interval; \n:lmeanci - lower confidence interval for mean; \n:umeanci - lower confidence interval for mean; \n:median - median,;\n:min - minimum; \n:max - maximum; \n:range - range; \n:q1 - lower quartile;\n:q3 - upper quartile;\n:iqr - inter quartile range; \n:kurt - kurtosis;\n:skew - skewness; \n:harmmean - harmonic mean; \n:ses standard error of skewness; \n:sek - standard error of kurtosis; \n:sum - sum.\n\n\n\n\n\ndescriptives(data; vars = nothing, sort = nothing, kwargs...)\n\nIf vars is nothing - try to include all collumns with numbers.\n\n\n\n\n\ndescriptives(data::DataSet{T}; kwargs...) where T <: ObsData\n\nDescriptive statistics for dataimport structure, see  dataimport.\n\n\n\n\n\n","category":"function"},{"location":"#MetidaStats","page":"Home","title":"MetidaStats","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = MetidaStats","category":"page"},{"location":"","page":"Home","title":"Home","text":"Metida descriptive statistics - provide tables with categirized descriptive statistics from tabular data.","category":"page"},{"location":"","page":"Home","title":"Home","text":"*This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n        \"index.md\",\n        \"api.md\",\n      ]\nDepth = 3","category":"page"},{"location":"#Example","page":"Home","title":"Example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using MetidaStats, CSV, DataFrames;\n\nds = CSV.File(joinpath(dirname(pathof(MetidaStats)), \"..\", \"test\", \"csv\",  \"ds.csv\")) |> DataFrame\n\nnothing; # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"For DataFrame ds:","category":"page"},{"location":"","page":"Home","title":"Home","text":"ds[1:5, :]","category":"page"},{"location":"#Import:","page":"Home","title":"Import:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"di  = MetidaStats.dataimport(ds, vars = [:var1, :var2], sort = [:col, :row])","category":"page"},{"location":"#Statistics:","page":"Home","title":"Statistics:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"des = MetidaStats.descriptives(di; skipmissing = true, skipnonpositive = true, stats = MetidaStats.STATLIST)","category":"page"},{"location":"#Make-DataFrame","page":"Home","title":"Make DataFrame","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"df = DataFrame(des)","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Textbooks:","category":"page"},{"location":"","page":"Home","title":"Home","text":"https://towardsdatascience.com/5-free-books-to-learn-statistics-for-data-science-768d27b8215","category":"page"},{"location":"","page":"Home","title":"Home","text":"Statistics for Julia:","category":"page"},{"location":"","page":"Home","title":"Home","text":"https://statisticswithjulia.org/","category":"page"}]
}
