# ExactDiagonalization.jl

| **Documentation** | **Build Status** | **Code Coverage** |
|:-----------------:|:----------------:|:-----------------:|
| [![**STABLE**][docs-stable-img]][docs-stable-url] [![**DEV**][docs-dev-img]][docs-dev-url] | [![Build Status][travis-img]][travis-url] [![Build Status][appveyor-img]][appveyor-url] | [![Code Coverage][codecov-img]][codecov-url] [![Code Coverage][coveralls-img]][coveralls-url] |

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

[travis-img]: https://travis-ci.org/kyungminlee/ExactDiagonalization.jl.svg?branch=master
[travis-url]: https://travis-ci.org/kyungminlee/ExactDiagonalization.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/r5270ufhu14imba0?svg=true
[appveyor-url]: https://ci.appveyor.com/project/kyungminlee/exactdiagonalization-jl

[codecov-img]: https://codecov.io/gh/kyungminlee/ExactDiagonalization.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/ExactDiagonalization.jl

[coveralls-img]: https://coveralls.io/repos/github/kyungminlee/ExactDiagonalization.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/kyungminlee/ExactDiagonalization.jl?branch=master
