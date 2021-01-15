# ExactDiagonalization.jl

**Documentation**: [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] 
**Build Status**: [![Build Test Submit][githubaction-img]][githubaction-url]
**Code Coverage**: [![Code Coverage][codecov-img]][codecov-url]

`ExactDiagonalization.jl` is a library for constructing quantum many-body Hamiltonians. It aims to provide
- convenient and efficient representation of a generic lattice Hamiltonian and wave function
- reduction of the Hilbert space dimension using symmetry

## Installation

To install, type the following in Julia's package Pkg REPL-mode:
```julia-repl
(v1.3) pkg> registry add https://github.com/kyungminlee/KyungminLeeRegistry.jl.git
(v1.3) pkg> add ExactDiagonalization
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://kyungminlee.org/ExactDiagonalization.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: http://kyungminlee.org/ExactDiagonalization.jl/dev

[githubaction-img]: https://github.com/kyungminlee/MathExpr.jl/workflows/Build/badge.svg
[githubaction-url]: https://github.com/kyungminlee/ExactDiagonalization.jl/actions?query=workflow%3ABuild

[codecov-img]: https://codecov.io/gh/kyungminlee/ExactDiagonalization.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/ExactDiagonalization.jl
