# Lean4 Language Targeting Strategy

Date: 2026-05-15
Task: `lean4-language-targeting-plan`

This document defines how future povu language targets should be built from the
Lean-validated core without turning the proof into background inspiration for
independent rewrites.  It is a strategy artifact only: it does not scaffold a
new target, generate target-language code, or change the Lean proof modules,
conformance harness, or fixture corpus.

## Readiness Gate

Language-target work starts only after the proof and conformance milestones are
complete.  The current gate is satisfied by the completed upstream chain:

- `lean4-proof-synthesis-check` integrated the trusted proof path and exposed
  `PovuLean.Pipeline.semanticGfaToVcf_correct`.
- `docs/lean4-proof/proof_obligation_check.md` records that `lake build`
  succeeded for the integrated Lean project and that the trusted `PovuLean`
  path had no `sorry` or `admit` tokens at that checkpoint.
- `lean4-conformance-harness` added the Rust/Lean boundary harness documented
  in `docs/lean4-proof/conformance.md`.
- `lean4-e2e-validation-corpus` expanded the fixture suite documented in
  `docs/lean4-proof/e2e_validation.md` and is complete in WG before this task.

These gates do not mean a future port is automatically verified.  They mean a
target has a proof-owned semantic oracle and a concrete corpus to be tested
against.  Any new language target must consume that oracle or produce an
equivalent proof-backed artifact before release.

## Lineage

This strategy inherits the trust boundary from the earlier Lean proof plan:

- `docs/lean4-proof/architecture.md` chose "Lean as executable oracle" as the
  near-term strategy and warned that future targets must either call
  extracted/ported verified algorithms or pass the same conformance suite.
- `docs/lean4-proof/proof_obligation_check.md` narrowed the trusted theorem
  boundary to semantic GFA documents, certified traversal/class/scan witnesses,
  supported hierarchy inputs, reference call sets, and semantic VCF emission.
- `docs/lean4-proof/conformance.md` and
  `docs/lean4-proof/e2e_validation.md` define the external implementation
  checks that bridge actual `povu gfa2vcf` behavior to Lean semantic VCF
  records.

The design lineage for future targets is therefore:

1. Lean definitions and theorems own the semantic contract.
2. The conformance harness owns evidence that an implementation matches that
   contract at the system boundary.
3. Language targets inherit proof guarantees only through extraction,
   proof-preserving reuse, or conformance against the Lean oracle.

## Non-Goals

- Do not treat Lean as a prose spec that can be reinterpreted independently in
  each language.
- Do not add a second target's fixtures, adapters, or implementation before a
  first target has passed its release gate or a human explicitly changes the
  target order.
- Do not broaden the trusted computing base silently.  Parsers, FFI layers,
  serializers, generated-code compilers, runtimes, and harness adapters are
  trusted or externally validated only when documented for that target.
- Do not claim that byte-level GFA parsing, current C++ control flow, thread
  scheduling, VCF text formatting, or a future port are proved by
  `semanticGfaToVcf_correct`.  Those remain outside the theorem unless a later
  proof explicitly brings them in.

## Strategy Options

| Strategy | Trust boundary | Proof guarantee preserved | Maintenance cost | Validation requirements | Fit |
| --- | --- | --- | --- | --- | --- |
| Lean extraction or code generation | Trust Lean kernel and theorem statements, plus the extraction/codegen path, target compiler, runtime, FFI, and any handwritten wrappers. | Strongest path if the generated artifact is directly derived from proved Lean definitions and the wrapper boundary is small.  The proof applies to the semantic core, not automatically to IO or formatting wrappers. | Medium to high upfront cost.  Build reproducibility, generated-code review, runtime integration, and performance tuning become permanent maintenance items. | `lake build`; no trusted-path placeholders; reproducible generated artifact; documented generated-code TCB; target integration tests; full Lean conformance corpus; wrapper-specific tests for parsing, formatting, and error boundaries. | Best long-term option for a small pure core or a service-style oracle.  Not the first production target until extraction performance and wrapper boundaries are known. |
| Manual port from the Lean reference algorithm | Trust the human port, target compiler/runtime, and target libraries.  Lean proves the source semantics, not the copied code. | Preserves the proof as an executable specification and review guide, but the port itself is unverified unless separately proved or exhaustively connected by conformance. | High ongoing cost.  Every Lean semantic change must be reflected in the port and its adapter tests.  Divergence risk grows with each optimized target-specific rewrite. | A target contract mapping target inputs/outputs to Lean semantic objects; full fixture corpus; negative-boundary tests; property or differential tests where generators exist; mismatch classification; CI that runs target tests and Lean conformance together. | Acceptable for performance-sensitive languages only after the conformance gate is mature.  Must not be marketed as independently proved. |
| Lean oracle with per-language conformance suite | Trust Lean as the oracle and the target as a black-box implementation under test.  Harness adapters and normalization rules are part of the external TCB. | Preserves the proof guarantee at the comparison boundary: a released target is evidence-backed to match the proved semantics on the maintained corpus, but it is not a proof of the target implementation. | Medium.  The main cost is keeping corpus fixtures, adapters, diagnostics, and target runners in sync with the Lean semantic boundary. | Full existing fixture corpus; target-specific runner; semantic normalization identical to the documented harness rules; controlled failures for unsupported inputs; diagnostics that expose differing semantic fields; CI blocking release on mismatches. | Best first strategy.  It keeps the proof central while allowing practical language work before extraction is production-ready. |
| Verified core behind FFI or service boundary | Trust the verified or conformance-checked core, the ABI/service protocol, serialization, target binding layer, and deployment runtime. | Preserves core guarantees if target code delegates semantic decisions to the checked core instead of reimplementing them.  Binding correctness still requires tests. | Medium.  Avoids duplicate algorithms but adds ABI/versioning, lifetime, memory, and packaging maintenance. | Core conformance suite; binding-level tests for all fixture classes; ABI compatibility tests; failure propagation tests; semantic output comparison after crossing the boundary. | Good for language bindings that should expose povu behavior without becoming algorithm ports. |

The recommended near-term pattern is the third option, with the fourth option
for bindings: use Lean as the release oracle and require each target to pass a
target-specific conformance suite before release.  Extraction or generation can
be introduced later for stable pure functions, but it should tighten the trust
boundary rather than replace the conformance gate.

## Target Release Gate

A target language can be released only when all of the following are true.

1. Proof baseline is current.
   - `lake build` succeeds.
   - Trusted `PovuLean/**` modules remain free of `sorry` and `admit`.
   - The target names the exact Lean theorem or reference function boundary it
     claims to implement, starting with
     `PovuLean.Pipeline.semanticGfaToVcf_correct` for GFA-to-VCF semantics.

2. Target contract is explicit.
   - The target states whether it is an extraction, generated artifact,
     binding/FFI layer, manual port, or black-box implementation under
     conformance.
   - The accepted input subset, unsupported-input behavior, and normalization
     rules match `docs/lean4-proof/conformance.md` unless the target explicitly
     defines a narrower release surface.
   - The target maps its outputs to Lean semantic objects: semantic VCF records,
     genotype calls, allele strings, traversal strings, ordering keys, INFO
     fields, and controlled failure classes.

3. Conformance suite passes.
   - The existing end-to-end fixtures in `tests/lean4_conformance/fixtures/`
     pass for the target surface or are documented as out of scope before any
     release claim.
   - The target covers positive VCF cases, header-only/no-variant behavior,
     deterministic ordering where required, hairpin/SUBR behavior, nested
     flubble hierarchy fields, and unsupported or malformed input boundaries.
   - Mismatches block release until classified as implementation bug, Lean spec
     bug, unsupported input, formatting-only difference, or intentionally
     external behavior.  Formatting-only differences still require documented
     normalization.

4. Target-native tests pass.
   - The target's unit and integration tests pass under its normal package
     manager.
   - Binding targets also test memory ownership, error propagation, and ABI
     stability at the boundary.
   - Manual ports include regression tests for every Lean fixture and at least
     one target-native test for each public API that exposes proof-owned
     behavior.

5. Release language is precise.
   - "Verified" is reserved for code paths that are either extracted from proved
     Lean definitions with a documented extraction TCB, or separately proved.
   - "Conforms to the Lean-validated semantics" is the right claim for manual
     ports and black-box implementations that pass the oracle suite.
   - "Binding to the checked povu core" is the right claim for wrappers that
     delegate to a conformance-checked implementation.

## First Target Recommendation: Rust

Rust should be the first language target track.  This does not mean starting a
native Rust rewrite of the algorithms.  The first Rust milestone should be a
Lean-backed conformance contract for the existing Rust binding surface, then a
small adapter/runner that proves the Rust-facing API reaches the same semantic
boundary as the current CLI or C++ core.

Rust is the right first target because:

- `povu-rs/` already exists, so the target track can validate a real binding
  surface instead of bootstrapping a new language ecosystem.
- The current Lean conformance harness is already Rust, so target-specific
  fixture orchestration can reuse Cargo, Rust diagnostics, and the existing
  normalization model.
- Rust's type system is useful at the FFI and semantic-adapter boundary, where
  ownership, controlled failures, and record normalization need to be explicit.
- A Rust target can start as a binding/conformance target and later choose
  between generated Lean-derived code or a manual native port if performance or
  deployment needs justify the larger trust boundary.

Concrete prerequisites for the first Rust target task:

- The proof and conformance milestones listed in the readiness gate are done.
- The Rust task must name a single release surface before implementation:
  binding-level API, CLI wrapper, or native algorithm module.  The recommended
  first surface is binding-level API because it minimizes algorithm divergence.
- The Rust task must define how Rust outputs are converted into the same
  semantic comparison objects used by `tests/lean4_conformance/src/main.rs`.
- The Rust task must reuse the existing fixture corpus before adding any new
  Rust-only fixtures.
- The Rust task must not edit `PovuLean/**`, the existing Lean conformance
  fixture corpus, or non-Rust language targets.

## Follow-Up Policy

Because the Lean proof/conformance prerequisites are complete, WG may add
Rust-only follow-up work after this plan.  Such tasks should depend on
`lean4-language-targeting-plan` and should include validation that runs both the
Rust target tests and the Lean-backed conformance gate relevant to the chosen
surface.

No Python, TypeScript, WebAssembly, C API, or other target task should be added
from this plan until one of these happens:

- the Rust target passes its release gate;
- a human explicitly changes the target priority;
- extraction/code generation becomes stable enough to create a shared generated
  core that multiple target bindings can consume without duplicating algorithm
  logic.

## Rust Target Validation Template

Any Rust follow-up task created from this plan should include at least this
validation shape:

```markdown
## Validation
- [ ] The Rust target contract cites `docs/lean4-proof/language_targets.md`
      and preserves the `PovuLean.Pipeline.semanticGfaToVcf_correct` boundary.
- [ ] The selected Rust surface maps its output to the same semantic VCF
      comparison objects used by `tests/lean4_conformance/src/main.rs`, or the
      task explicitly documents why the surface is narrower.
- [ ] `lake build` passes.
- [ ] `cargo test --manifest-path tests/lean4_conformance/Cargo.toml` passes.
- [ ] The relevant Rust target tests pass, for example
      `cargo test --manifest-path povu-rs/Cargo.toml` when the work touches
      `povu-rs/`.
- [ ] The full target fixture command is documented and passes, or every
      missing fixture is documented as an unreleased blocker rather than a
      skipped release criterion.
```

Implementation tasks should add a failing target-side test before changing the
target implementation or adapter.  Documentation-only Rust contract tasks may
replace that item with a completeness review against this strategy document and
the existing conformance docs.

## Decision Summary

Future language targets should be target implementations under a Lean-owned
contract, not independent implementations with similar wording.  The initial
release path is:

1. Keep Lean as the semantic oracle.
2. Require per-target conformance at the GFA-to-VCF boundary.
3. Start with Rust because the repository already has Rust bindings and a Rust
   conformance harness.
4. Defer new language branches until Rust has a passing target gate or a human
   reprioritizes the graph.
5. Consider extraction or generated code later as a way to reduce divergence,
   while keeping conformance as the release guard.
