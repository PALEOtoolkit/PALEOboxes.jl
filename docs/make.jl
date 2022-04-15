using Documenter

import PALEOboxes

using DocumenterCitations

bib = CitationBibliography("src/paleo_references.bib")

makedocs(bib, sitename="PALEOboxes Documentation", 
        pages = [
            "Home" => "index.md",
            "Design" => [
                "DesignOverview.md",
                "CreateInitializeLoop.md"
            ],
            "Reference" => [
                "PALEOboxes.md",
            ],
            "References.md",
            "indexpage.md",
        ],
        format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
        ),
        # repo = "https://github.com/PALEOtoolkit/PALEOboxes.jl/blob/master/{path}#{line}"
    )

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOboxes.jl.git",
)