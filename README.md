# MetidaStats.jl
Metida descriptive statistics.



| Status | Cover | Build | Docs |
|--------|-------|-------|------|
|[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)|[![codecov](https://codecov.io/gh/PharmCat/MetidaStats.jl/branch/main/graph/badge.svg?token=A9eyT9g0WZ)](https://codecov.io/gh/PharmCat/MetidaStats.jl)|![Tier 1](https://github.com/PharmCat/MetidaStats.jl/workflows/Tier%201/badge.svg) | [![Latest docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://pharmcat.github.io/MetidaStats.jl/dev/) [![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://pharmcat.github.io/MetidaStats.jl/stable/)|
## Install

```
import Pkg; Pkg.add(url = "https://github.com/PharmCat/MetidaStats.jl.git")
```

## Import DataFrame

```
data = CSV.File("somedata.csv") |> DataFrame

# variables to analyze
vars = [:Cmax, :AUClast]

# sorting variables 
sort = [:form, :period]

ds = dataimport(data; vars = vars, sort = sort)
```

## Get descriptive statistics

```
descriptives(ds, stats = [:n, :mean, :var])
```

## Or without dataimport step

```
descriptives(data; vars = vars, sort = sort, stats = [:n, :mean, :var])
```

Keywords:

- `skipmissing` - drop NaN and Missing values, default = true;
- `skipnonpositive` - drop non-positive values (and NaN, Missing) for "log-statistics" - :geom, :geomean, :logmean, :logvar, :geocv;
- `stats` - default set `stats = [:n, :mean, :sd, :se, :median, :min, :max]`;
- `corrected` - use corrected var (true);
- `level` - level for confidence intervals (0.95);

Possible values for `stats` is: 

* :n - number of observbations;
* :posn - positive (non-negative) number of observations;
* :mean - arithmetic mean;
* :var - variance;
* :bvar - variance with no correction;
* :geom - geometric mean; 
* :logmean - arithmetic mean for log-transformed data;
* :logvar - variance for log-transformed data;
* :sd - standard deviation (or σ);
* :se - standard error; 
* :cv - coefficient of variation; 
* :geocv - coefficient of variation for log-transformed data;
* :lci - lower confidence interval;
* :uci - upper confidence interval; 
* :lmeanci - lower confidence interval for mean; 
* :umeanci - lower confidence interval for mean; 
* :median - median;
* :min - minimum; 
* :max - maximum; 
* :range - range; 
* :q1 - lower quartile;
* :q3 - upper quartile;
* :iqr - inter quartile range; 
* :kurt - kurtosis;
* :skew - skewness; 
* :harmmean - harmonic mean; 
* :ses standard error of skewness; 
* :sek - standard error of kurtosis; 
* :sum - sum.