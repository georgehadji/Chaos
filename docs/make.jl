# Build the user documentation: julia --project=docs docs/make.jl (output in docs/build, git-ignored).
using Documenter, Chaoslab

const LAYERS = [Chaoslab.Foundation, Chaoslab.Systems, Chaoslab.Tangent, Chaoslab.Diagnostics,
                Chaoslab.Lyapunov, Chaoslab.Indicators, Chaoslab.Experiments, Chaoslab.Analysis]

makedocs(; sitename = "Chaoslab", modules = [Chaoslab; LAYERS], checkdocs = :exports,
         format = Documenter.HTML(; prettyurls = false), repo = Remotes.GitHub("georgehadji", "Chaos"),
         pages = ["Home" => "index.md", "API" => "api.md"])
