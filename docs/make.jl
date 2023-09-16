using Documenter

import PALEOboxes

using DocumenterCitations

bib = CitationBibliography(
        joinpath(@__DIR__, "src/paleo_references.bib"),
        style=:authoryear,
    )

makedocs(;
        sitename="PALEOboxes Documentation", 
        pages = [
            "Home" => "index.md",
            "Design" => [
                "DesignOverview.md",
                "CreateInitializeLoop.md"
            ],
            "Reference" => [
                "DomainsVariablesFields.md",
                "Solver API.md",
                "Reaction API.md",
                "ReactionCatalog.md",
            ],
            "References.md",
            "indexpage.md",
        ],
        format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
        ),
        plugins = [bib],
        # repo = "https://github.com/PALEOtoolkit/PALEOboxes.jl/blob/master/{path}#{line}"
    )

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOboxes.jl.git",
)
