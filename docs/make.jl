using Documenter
using ExactDiagonalization

makedocs(
    modules=[ExactDiagonalization],
    doctest=true,
    sitename="ExactDiagonalization.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
  )

deploydocs(
    deps=Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo = "github.com/kyungminlee/ExactDiagonalization.jl.git",
  )