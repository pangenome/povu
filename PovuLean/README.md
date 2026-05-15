# PovuLean Scaffold

This directory contains the initial Lean4 workspace for the povu proof effort.
The formalization architecture is documented in
[`docs/lean4-proof/architecture.md`](../docs/lean4-proof/architecture.md).

## Build

The Lean toolchain is pinned in [`../lean-toolchain`](../lean-toolchain).
From the repository root, build the proof workspace with:

```sh
lake build
```

`lake build` is the canonical gate for scaffold and future proof-module
changes. Lake build output is written under `.lake/` and is intentionally not
tracked.

## Module Layout

All trusted proof modules live under namespace `PovuLean`.

- [`../PovuLean.lean`](../PovuLean.lean) is the root import aggregator.
- [`Core.lean`](Core.lean) aggregates foundational graph, path, finite, and DFS
  modules owned by `lean4-core-graph-model`.
- [`GFA.lean`](GFA.lean) aggregates semantic GFA modules owned by
  `lean4-gfa-spec`.
- [`Algorithms.lean`](Algorithms.lean) aggregates flubble, hairpin, and
  flubble-tree algorithm modules owned by their downstream proof tasks.
- [`VCF.lean`](VCF.lean) aggregates semantic VCF modules owned by
  `lean4-vcf-semantics`.

Downstream tasks should add substantive modules under their owned module
families, then import completed trusted modules from the appropriate
aggregator. Draft modules with placeholder proof obligations must stay outside
the trusted import path and must be documented by the owning task.

This scaffold deliberately contains no placeholder proof terms.
