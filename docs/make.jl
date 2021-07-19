using Documenter, MetidaStats, Weave, PrettyTables, CSV, DataFrames
#using DocumenterLaTeX


makedocs(
        modules = [MetidaStats],
        sitename = "MetidaStats.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            ],
        )


deploydocs(repo = "github.com/PharmCat/MetidaStats.jl.git", devbranch = "main", forcepush = true
)
