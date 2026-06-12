# Flubble Linear-Time Batch Quality Pass

Date: 2026-06-12
Task: `quality-pass-flubble`

This note records the staging decisions applied to the flubble linear-time proof
batch before downstream implementation work starts. The WG task graph remains
the source of truth for executable task descriptions and dependencies.

## Scope Boundaries

- The existing semantic flubble detector is treated as a verified reference
  specification, not as evidence of linear runtime.
- Count bounds and hierarchy output-size bounds must land before any
  full-linear-time theorem depends on emitted flubble or hierarchy sizes.
- Cost-model work introduces vocabulary and named assumptions only. It must not
  assert current C++ or Rust runtime behavior.
- A verified indexed detector may support a one-pass Lean cost argument, but
  only under explicit class-index/map operation contracts.
- C++ and Rust conformance tests are implementation evidence. They are not a
  formal runtime proof unless a later machine-checkable bridge connects the
  implementation to the Lean cost contracts.

## Dependency Staging

The downstream graph now enforces this order for the core proof chain:

1. `lean-add-flubble`: define the cost/asymptotic vocabulary and named stage
   assumptions.
2. `lean-prove-flubble`: prove count bounds for the current semantic detector.
3. `lean-prove-flubble-2`: prove hierarchy output-size bounds.
4. `lean-implement-indexed`: implement/prove an indexed detector after count
   and hierarchy bounds are available.
5. `lean-prove-flubble-3`: prove flubble extraction and hierarchy construction
   linear only under the available detector and cost assumptions.
6. `lean-compose-conditional`: compose a conditional semantic Lean
   decomposition theorem from explicit stage contracts.
7. `bridge-povu-implementation`: audit C++/Rust implementation conformance
   separately from formal Lean proof claims.
8. `synthesis-flubble-linear`: summarize what is formally proved,
   conditionally proved, conformance-tested, and still open.

## Validation Requirements Added

Every downstream Lean code/proof task now requires:

- `lake build`;
- a placeholder scan with `rg -n '\b(sorry|admit|axiom|opaque|unsafe)\b'`
  over the relevant trusted Lean modules;
- explicit documentation of any remaining assumptions instead of hidden axioms;
- wording that prevents current povu implementation runtime from being claimed
  without a concrete cost model and implementation bridge.

The bridge and synthesis tasks additionally require separating:

- formally proved Lean theorems;
- conditional Lean theorems;
- C++/Rust conformance test evidence;
- open implementation or proof obligations.

