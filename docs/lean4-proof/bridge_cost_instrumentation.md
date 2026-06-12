# Bridge Cost Instrumentation

Date: 2026-06-12
Task: `bridge-cost-instrumentation`

This note documents the opt-in C++ stage-cost trace points added for later
proof-bridge work. The counters are implementation evidence only. They do not
prove that current povu satisfies the Lean linear-time contracts; they provide
stable names and coarse size proxies that future checkers can compare against
Lean `StageCost` records and `ExtractionPipelineCosts` fields.

## Enabling

Instrumentation is disabled by default and does not write to `stdout` or
`stderr` during ordinary `povu`, C++ test, or Lean/Rust conformance runs.

Set `POVU_STAGE_COST_TRACE=1` to enable it:

```bash
POVU_STAGE_COST_TRACE=1 povu gfa2vcf -i input.gfa -P ref-prefix
```

When enabled, the process writes a compact report to `stderr` at exit:

```text
povu-stage-cost-trace schema=povu.stage-cost.v1
povu-stage-cost contract=componentDecomposition trace=component-decomposition calls=1 input_items=42 output_items=1 elapsed_ns=1234
```

The report format is intentionally line-oriented so conformance harnesses can
keep using their existing `stdout` artifacts. The default harness output is
unchanged unless the environment variable is explicitly set.

## Counter Names

| Trace point | Stage-contract name | Lean cost slot | Current implementation hook | Size proxy |
| --- | --- | --- | --- | --- |
| Component decomposition | `componentDecomposition` | `SemanticDecomposePipelineCosts.componentDecomposition` and `componentDecompositionLinear` | `bd::VG::componetize` | input graph vertices + edges; output component count |
| Tip/dummy augmentation | `tipDummyAugmentation` | `SemanticDecomposePipelineCosts.tipDummyAugmentation` and `tipDummyAugmentationLinear` | `pst::Tree::from_bd` | component vertices + edges + tips; output tree vertices + tree edges + backedges |
| Traversal-frame construction | `traversalFrameConstruction` | `SemanticDecomposePipelineCosts.traversalFrameConstruction` and `traversalFrameConstructionLinear` | `compute_eq_class_stack` | spanning-tree vertices + tree edges; output candidate stack entries |
| Cycle-class assignment | `cycleClassAssignment` | `SemanticDecomposePipelineCosts.cycleClassAssignment` and `cycleClassAssignmentLinear` | `simple_cycle_equiv` over `handle_vertex` | spanning-tree vertices + tree edges + initial backedges; output labeled tree/backedge count |
| Flubble boundary extraction | `boundaryEmission` | `ExtractionPipelineCosts.boundaryEmission` in `ConditionalExtractionCostContract` | `compute_eq_class_metadata` | candidate stack entries; output close-after-gap boundary candidates |
| Flubble hierarchy construction | `hierarchyConstruction` | `ExtractionPipelineCosts.hierarchyConstruction` in `ConditionalExtractionCostContract` | `add_flubbles` | candidate stack entries; output PVST vertex count |

The trace names are human-readable labels. The `contract=` names are the stable
bridge identifiers and should be used by downstream proof/checker code.

## Scope And Limits

The counters aggregate per-process calls, input size proxies, output size
proxies, and wall-clock nanoseconds. Multi-component and multi-threaded runs
therefore produce totals for all observed calls in the process.

These counters are deliberately not a formal cost model:

- container operation costs for `std::vector`, `std::map`, `std::unordered_map`,
  `std::unordered_set`, `std::set`, `std::list`, and `std::stack` are still
  external assumptions for a runtime proof;
- elapsed time is diagnostic only and should not be used as a Lean `Nat` cost
  witness;
- the size proxies name the implementation surfaces that must later be related
  to Lean graph, frame, candidate-stack, boundary, and hierarchy sizes.

The intended proof-bridge workflow is to use the `contract=` names to align C++
execution traces with the Lean records in `PovuLean/Complexity/Decomposition.lean`
and `PovuLean/Complexity/Flubble.lean`, then separately prove or assume the
operation-level bounds needed for each stage.
