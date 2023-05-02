# MetidaStats

```@meta
CurrentModule = MetidaStats
```

Metida descriptive statistics.

*This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.


## Contents

```@contents
Pages = [
        "index.md",
        "api.md",
      ]
Depth = 3
```

## Example

```@example dsexample
using MetidaStats, CSV, DataFrames;

ds = CSV.File(joinpath(dirname(pathof(MetidaStats)), "..", "test", "csv",  "ds.csv")) |> DataFrame

nothing; # hide
```

For DataFrame `ds`:

```@example dsexample
ds[1:5, :]
```

### Import:

```
di  = MetidaStats.dataimport(ds, vars = [:var1, :var2], sort = [:col, :row])
```

### Statistics:

```
des = MetidaStats.descriptives(di; skipmissing = true, skipnonpositive = true, stats = MetidaStats.STATLIST)
```

### Make DataFrame

```
df = DataFrame(des)
```


## Reference


Textbooks:

https://towardsdatascience.com/5-free-books-to-learn-statistics-for-data-science-768d27b8215

Statistics for Julia:

https://statisticswithjulia.org/

