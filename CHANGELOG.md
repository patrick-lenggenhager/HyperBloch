# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.1] - 2024-03-05

## Added
- Add a change log (this file).
- Add link to getting-started guide and to Wolfram Community post to README.md.
- Add warning message if NCAlgebra paclet version 6+ is not installed.
- Add option for `AbelianBlochHamiltonianExpression` to return the Bloch Hamiltonian as a
`SparseArray`.

## Changed
- Improve performance of `AbelianBlochHamiltonianExpression`.
- Update references in readme.
- Load example data from files in example data directory instead of embedding the data
in the documentation notebooks.

## Fixed
- Fixed `Det::luc` error message being raised by an unnecessary check in `Det` for
showing triangle tessellations.


## [0.9.0] - 2023-11-29

Initial release.