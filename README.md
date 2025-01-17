# Backboner

[![Latest Release](https://img.shields.io/github/release/MurrellGroup/Backboner.jl.svg)](https://github.com/MurrellGroup/Backboner.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://MurrellGroup.github.io/Backboner.jl/dev/)
[![Build Status](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/Backboner.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/Backboner.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/Backboner.jl)

Backboner is a Julia package that offers a suite of tools for storing protein backbone atom positions, estimating oxygen atom positions, assigning secondary structure, and more.

## Installation

Backboner is a registered Julia package, and can be installed with the Julia package manager:

```julia
using Pkg
Pkg.add("Backboner")
```

## Types and functions

The `Protein` type wraps a vector of `Chain`s, which in turn wraps the `Backbone{4}` type (4, because it stores the positions of 4 atoms per residue: N, CA, C, O). The `Backbone{N}` type has the `N` type parameter in order to remain flexible. It allows one pass only the N, CA, and C atoms of a backbone, such that the O atom positions can added in using the `add_oxygens` function.

The secondary structure of an entire chain is described by a `Vector{Char}`, where '-' stands for coil/loop, 'H' for helix, and 'E' for strand. For assignment of secondary structure, this package uses the [AssigningSecondaryStructure.jl](https://github.com/MurrellGroup/AssigningSecondaryStructure.jl) package, which implements a simplified version of the DSSP algorithm.

Protein backbones can be loaded from a PDB file using the `pdb_to_protein` function, which returns a `Protein` instance. Inversely, a `Protein` instance can be written to a PDB file using the `protein_to_pdb` function.

## Example

```julia
julia> using Backboner

julia> protein = pdb_to_protein("test/data/1ZAK.pdb")
2-element Protein{Float32}:
 Chain A with 220 residues
 Chain B with 220 residues

julia> chain = protein["A"]
Chain A with 220 residues

julia> chain.backbone
3×4×220 Backbone{4, Float32}:
[:, :, 1] =
 22.346  22.901  23.227  22.689
 17.547  18.031  16.793  15.72
 23.294  21.993  21.163  21.448

;;; … 

[:, :, 220] =
 21.808  22.263  21.085  19.939
 13.861  13.862  14.233  13.851
  2.734   1.355   0.446   0.791
```