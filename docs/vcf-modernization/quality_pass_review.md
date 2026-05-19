# VCF Modernization Quality Pass Review

Task: `vcf-modern-quality-pass`
Date: 2026-05-18

## Evaluation Scope

This review evaluates the practical VCF modernization task batch before
downstream workers start. The grade applies to the batch readiness after this
quality pass, not to the downstream implementation or research results that
have not run yet.

The task rubric was not underspecified. It named five concrete acceptance
criteria: downstream scope/dependencies/validation/file ownership, Rust-first
targeting, vg/vcfbub/vcfwave study ordering, Lean/Rust guardrail blocking, and
fan-out convergence into `vcf-modern-synthesis-check`.

## Changes Made

- Added explicit `## Dependency expectations` and `## File ownership /
  boundaries` sections to all eight downstream tasks in scope.
- Kept Rust as the implementation target and documented C++ as context only.
- Preserved the ordering where `vcf-modern-vg-pipeline-study`,
  `vcf-modern-rust-output-inventory`, and `vcf-modern-lean-rust-guardrail`
  precede `vcf-modern-output-spec`.
- Added a direct dependency from `vcf-modern-lean-rust-guardrail` to
  `vcf-modern-synthesis-check`, making the Lean/Rust guardrail an explicit
  final-claim blocker instead of only a transitive prerequisite.
- Confirmed all implementation/design branches feed into
  `vcf-modern-synthesis-check`.

## Dimension Scores

| Dimension | Score | Rationale |
| --- | ---: | --- |
| Downstream scope clarity | 0.94 | Each downstream task now has a concrete deliverable, task-specific validation, and explicit ownership boundaries. Residual risk is that Rust emitter and conformance corpus may still need coordination if fixture files become shared. |
| Dependency graph quality | 0.93 | The graph enforces research/inventory/guardrail before the output spec, then fans out to emitter, corpus, and untangling, then joins in synthesis. The added direct guardrail-to-synthesis edge improves final-claim blocking. |
| Rust-first targeting | 0.96 | Rust is named as the implementation target throughout. C++ is limited to documentation, divergence notes, or explicitly logged unavoidable fixture context. |
| Upstream pipeline study gating | 0.92 | `vcf-modern-vg-pipeline-study` blocks `vcf-modern-output-spec`, and implementation depends on the spec. The study task requires upstream versions/commits or documentation dates, which is the right control for fast-moving external behavior. |
| Lean/Rust correctness guardrail | 0.94 | The guardrail is represented before the spec, before conformance work, and now directly before synthesis. The guardrail task also requires filing deeper intermediate-export work if the existing harness is insufficient. |
| Fan-out and synthesis convergence | 0.95 | The batch permits parallel research/design/implementation branches but all practical claims converge on `vcf-modern-synthesis-check`. The synthesis task owns the readiness report and must file unresolved mismatches before completion. |
| Validation specificity | 0.90 | Validation is concrete for research, design, implementation, and synthesis tasks. Remaining risk is normal for pre-work quality passes: exact commands for the later Rust/Lean conformance checks must be confirmed by the guardrail and corpus workers. |

Overall calibrated grade: **0.93**

Confidence: **0.82**

## Validation Checklist

- [x] Every downstream task has clear scope, dependencies, validation criteria,
  and file ownership expectations.
- [x] Rust is the implementation target; C++ changes or known prior C++ fixes
  are documented only unless explicitly justified.
- [x] vg deconstruct / vcfbub / vcfwave study precedes output-spec and
  implementation work.
- [x] Lean/Rust structure conformance is represented as a blocking guardrail
  before final claims about practical VCF output.
- [x] The graph permits fan-out and all branches feed into
  `vcf-modern-synthesis-check`.

## Residual Risks

- `vcf-modern-synthesis-check` had its direct dependencies changed, so WG
  cleared its current agent assignment and should reassign it when ready.
- The conformance corpus and Rust emitter tasks must avoid editing the same
  files concurrently. Their ownership sections now separate practical
  conformance fixtures from Rust implementation, but workers should file a
  blocking task if that boundary becomes false.
- No code was changed in this quality pass, so cargo and Lean validation were
  not applicable here; they remain required on downstream code/proof tasks.
