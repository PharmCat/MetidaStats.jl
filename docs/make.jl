using Documenter, MetidaStats, CSV, DataFrames
#using DocumenterLaTeX


makedocs(
        modules = [MetidaStats],
        sitename = "MetidaStats.jl",
        authors = "Vladimir Arnautov",
        pages = [
            "Home" => "index.md",
            "API" => "api.md",
            ],
        )


deploydocs(repo = "github.com/PharmCat/MetidaStats.jl.git", devbranch = "main", forcepush = true
)
