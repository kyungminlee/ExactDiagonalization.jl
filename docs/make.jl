using Documenter
using ExactDiagonalization

makedocs(
    modules=[ExactDiagonalization],
    doctest=true,
    sitename="ExactDiagonalization.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
    pages = [
      "Home" => "index.md",
      "Types" => Any[
        "Hilbert space" => "hilbertspace.md",
        "Operator" => "operator.md",
        "Representation" => "representation.md",
        "Symmetry" => "symmetry.md"
      ],
      "API" => "api.md",
    ]
  )

deploydocs(
    #deps=Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/kyungminlee/ExactDiagonalization.jl.git",
  )
