using Documenter

makedocs(
    sitename = "CrystalCumulants.jl",
    authors = "Ethan Meitz & Alois Castellano",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ejmeitz.github.io/CrystalCumulants.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "API Reference" => "api.md",
        "Theory" => "theory.md",
        "Example" => "example.md",
    ],
)

deploydocs(
    repo = "github.com/ejmeitz/CrystalCumulants.jl.git",
    devbranch = "main",
)
