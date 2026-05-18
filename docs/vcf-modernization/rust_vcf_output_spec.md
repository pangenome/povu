# Rust VCF Output Spec

Task: `vcf-modern-output-spec`

Date: 2026-05-18

## Scope and Inputs Consumed

This document specifies the first practical Rust VCF output contract for povu.
It is a contract for what the Rust emitter must produce, what stays in the
semantic call model before formatting, and what is intentionally delegated to
downstream normalization or untangling tools.

The spec consumes all three predecessor artifacts:

- `docs/vcf-modernization/vg_vcf_pipeline_study.md`: the upstream
  vg/vcfbub/vcfwave pipeline keeps graph-aware raw VCF separate from popping,
  left-normalization, WFA decomposition, clustering, and richer nested modes.
  This spec adopts that boundary. povu emits raw graph records with traversal
  provenance and hierarchy; vcfbub-like popping, `bcftools norm`, and
  vcfwave-like decomposition are downstream profiles or external steps.
- `docs/vcf-modernization/rust_vcf_output_inventory.md`: current `povu-rs`
  advertises `gfa_to_vcf` and `GraphAnalysis::write_vcf`, but both fail before
  producing VCF. Current C++ output conventions are compatibility context, not
  a byte-for-byte Rust target. This spec therefore defines a Rust semantic data
  model and writer behavior instead of treating current FFI stubs as output
  behavior.
- `docs/vcf-modernization/lean_rust_guardrail.md`: Lean proves semantic VCF
  emission from accepted semantic GFA, certified flubble/hairpin witnesses,
  supported hierarchy input, and well-formed `VCF.VariantCall`s. It does not
  prove current Rust/C++ intermediate structures or byte-level VCF formatting.
  This spec keeps the core output contract at that Lean `VariantCall` boundary
  and treats Rust extraction, headers, decimal formatting, and external
  normalization as conformance obligations.

The output contract is intentionally not a claim that current Rust binding VCF
generation is implemented, and not a claim that current implementation PVST or
flubble structures are Lean-conformant. Those claims require the separate
`vcf-modern-structure-export-conformance` work.

## Design Decisions

Default Rust output is a raw graph VCF. It preserves povu's proved structure
semantics as directly as VCF allows:

- one semantic `VariantCall` maps to one VCF record;
- supported sources are Lean-backed flubble-tree nodes and hairpin boundaries;
- graph traversal provenance is emitted in `AT`;
- flubble hierarchy is emitted as `ES` and `LV`;
- hairpin reverse substitutions are emitted as `SUBR` without flubble hierarchy
  fields;
- record order follows Lean's deterministic order key, not C++ chunk/map
  generation order;
- REF/ALT alleles are VCF-anchored but not left-normalized across repeats;
- complex allele splitting, large-parent popping, left alignment, and
  repetitive untangling are separate downstream steps.

The Rust emitter should support two output surfaces over the same semantic
document:

- `to_string()` or equivalent string output for tests and embedded callers;
- `write_path(path)` or equivalent file output for `povu::gfa_to_vcf`.

Both surfaces must serialize the same `VcfDocument`. The first supported Rust
implementation should emit a single combined VCF. C++-style split output under
an output directory is a compatibility feature, not part of the core semantic
contract.

## Supported VCF Subset

The Rust emitter writes VCF 4.2-style records. The semantic subset is the Lean
subset in `PovuLean/VCF/Spec.lean`, with explicit byte-formatting choices added
here.

Required columns:

| Column | Rust contract |
| --- | --- |
| `CHROM` | Selected reference coordinate name for the record. For current conformance fixtures this is the full reference path name such as `HG1#1#chr1` or `ref`. The emitter must not invent a coordinate system if the input lacks VCF-capable reference/path metadata. |
| `POS` | One-based positive coordinate of the emitted REF allele's first base. For `INS` and `DEL`, this is the anchor base coordinate. For `SUB` and `SUBR`, this is the first replaced reference base. |
| `ID` | Stable graph source id. Default format is the current oriented boundary/source string, for example `>0>3`, `>1>4`, or `>1>5`. The id is derived from canonical source endpoints after any explicit node-id translation. Duplicate source ids within one output are invalid. |
| `REF` | Nonempty reference allele string supplied by the semantic `VariantCall`. |
| `ALT` | One or more nonempty alternate allele strings. All alternates in a record must share the same `REF` and `VARTYPE`. Symbolic alleles such as `*` or `<INV>` are not in the first core subset. |
| `QUAL` | Always `60`. |
| `FILTER` | Always `PASS`. |
| `INFO` | Ordered, deterministic fields described below. Unknown, lossy, or downstream-only transformation state must not be hidden in existing fields. |
| `FORMAT` | Always `GT`. No `GQ`, `DP`, `AD`, or other per-sample fields are emitted in the first Rust contract. |

Core INFO fields, in serialization order:

| Field | Number/Type | Semantics |
| --- | --- | --- |
| `AC` | `Number=A`, `Type=Integer` | Alternate allele counts from the semantic call. |
| `AF` | `Number=A`, `Type=Float` | Alternate count divided by total allele count. The semantic model stores exact integer counts; text rendering must be deterministic and more precise than current C++ one-decimal rendering. Use `0.5` for 1/2 and a stable decimal such as `0.3333333333333333` for 1/3, not `0.3`. |
| `AN` | `Number=1`, `Type=Integer` | Total called alleles represented by the record: reference allele count plus all alternate counts. This corrects the current C++ header's `Type=String` declaration. |
| `NS` | `Number=1`, `Type=Integer` | Number of sample columns with at least one non-missing genotype allele. |
| `AT` | `Number=R`, `Type=String` | Traversal provenance for REF followed by each ALT. This is graph provenance, not a normalized coordinate path. |
| `VARTYPE` | `Number=1`, `Type=String` | One of `DEL`, `INS`, `SUB`, `SUBR`. |
| `TANGLED` | `Number=1`, `Type=String` | `T` or `F`. `T` marks a graph-faithful but representation-complex record that should not be mistaken for a primitive normalized variant. |
| `ES` | `Number=1`, `Type=String` | Present only for non-`SUBR` flubble records. It is the enclosing/source site id used by current povu/Lean fixtures. |
| `LV` | `Number=1`, `Type=Integer` | Present only for non-`SUBR` flubble records. `0` is top-level. Nested children use increasing levels. |

Header policy:

- emit exactly one `##fileformat=VCFv4.2`;
- emit `##source=povu-rs` or `##source=povu` with an additional
  implementation/version line when available;
- emit one `##FORMAT=<ID=GT,...>` line, not the duplicate `GT` lines currently
  produced by C++;
- emit INFO metadata matching the field types above;
- emit `##contig` lines for selected reference coordinate names when lengths
  are known; if lengths are unknown, emit ids without length or omit contigs
  consistently by configuration;
- do not emit dynamic `##fileDate` by default in deterministic test mode.
  A compatibility mode may include it, but semantic conformance must ignore it.

## Semantic Data Model Required in Rust

The Rust implementation needs explicit types equivalent to the Lean semantic
input, not ad hoc VCF strings:

- `VariantSource`: `Flubble { boundary/source id, parent?, level, order keys }`
  or `Hairpin { boundary/source id, order keys }`;
- `AlleleConstruction`: `Deletion { anchor, deleted }`,
  `Insertion { anchor, inserted }`, `Substitution { reference, alternate }`,
  and `ReverseSubstitution { reference, alternate }`;
- `AlternateAllele`: construction, traversal string, and count;
- `GenotypeColumn`: sample name plus ordered phased alleles
  (`missing`, `ref`, or one-based alternate index);
- `VariantCall`: source, chrom, contig order, positive POS, id, REF,
  reference traversal, alternates, variant type, tangled flag, optional ES/LV,
  reference allele count, and genotype columns;
- `Record`: the formatted semantic record produced from `VariantCall`.

The writer must validate the same well-formedness obligations as Lean before
serializing:

- `CHROM`, `ID`, `REF`, reference traversal, and all ALT/traversal strings are
  nonempty;
- position is positive;
- each ALT refers to an existing one-based genotype allele number;
- all alternates share the record REF and VARTYPE;
- source supports the variant type: flubbles emit only `DEL`, `INS`, or `SUB`;
  hairpins emit only `SUBR`;
- `ES`/`LV` are present for flubble records and absent for `SUBR`;
- total allele count is positive;
- sample names are nonempty and genotype columns have at least one phase.

If these checks fail, Rust must return a structured error. It must not emit an
empty VCF that looks successful and must not panic across the public API.

## Allele Anchoring and Normalization

The core Rust emitter uses graph-derived VCF anchoring only:

- `DEL(anchor, deleted)`: `REF = anchor + deleted`, `ALT = anchor`, `POS` is the
  anchor coordinate.
- `INS(anchor, inserted)`: `REF = anchor`, `ALT = anchor + inserted`, `POS` is
  the anchor coordinate.
- `SUB(reference, alternate)`: `REF = reference`, `ALT = alternate`, `POS` is
  the first replaced reference base.
- `SUBR(reference, alternate)`: same REF/ALT rule as `SUB`, but the source must
  be a hairpin and `VARTYPE=SUBR`.

The emitter must not left-shift indels across repeats, trim shared prefixes or
suffixes beyond the semantic construction supplied by the call, or split one
complex graph record into primitive records. Those operations rewrite POS,
REF, ALT, record count, and sometimes genotypes; they belong to downstream
normalization.

Recommended downstream profiles:

| Profile | Responsibility |
| --- | --- |
| `raw-graph` | Default. Emit Lean-aligned records with `AT`, `ES`/`LV`, anchored REF/ALT, and complex alleles intact. |
| `top-level-only` | A filter over semantic records that keeps `LV=0` flubble records and all `SUBR` records. It is not equivalent to vcfbub popping. |
| `popped` | vcfbub-like two-pass filtering by level, max allele length, and max REF length, with child rescue for popped parents. This should operate after raw semantic records are created. |
| `normalized` | External or later-library `bcftools norm`/vcfwave-style rewriting. It must record provenance that POS/REF/ALT and genotypes may have changed. |

Core Rust output must keep `AT` stable enough that a downstream normalizer can
trace rewritten alleles back to graph traversals.

## Nested and Complex Flubbles

Nested flubbles are represented directly in raw graph mode. The Rust emitter
must not flatten a nested site by inventing a primitive coordinate unless a
separate normalization pass can prove the rewrite.

Rules:

- emit a record for each well-formed semantic flubble call selected by the
  caller's `ReferenceCallSet`;
- use `LV=0` for top-level flubble records and increasing `LV` for nested
  records;
- use `ES` as the current povu/Lean enclosing/source site id;
- serialize child records at their semantic reference coordinate if the child
  has one on the selected reference path;
- if a nested child exists only inside an inserted haplotype and has no ordinary
  coordinate on the selected reference, the first Rust contract must either
  reject that child for VCF output with a structured unsupported error or keep
  it out of the VCF while preserving the structure for downstream graph-aware
  export. It must not silently place it at a misleading top-level coordinate.

For genotype semantics at nested sites, the internal Rust model must distinguish
at least these cases even if VCF text collapses some of them:

- called reference allele;
- called alternate allele;
- no traversal through this nested site;
- no-call or unavailable sample data;
- future explicit spanning deletion/star allele.

Default VCF serialization uses `.` for no traversal through a nested child,
matching current Lean/C++ conformance for `nested-deletion`. A future
`star-spanning` compatibility mode may emit `*` alleles like `vg -R`, but that
is not in the first core subset because Lean does not yet model `*` as a
variant allele class.

Complex flubbles:

- If all alternates in a site share the same REF and variant type, they may be
  emitted as one multi-ALT record.
- If a flubble contains heterogeneous allele classes or ambiguous internal
  decomposition, Rust should emit the smallest graph-faithful semantic record
  it can justify. Usually that means a `SUB` record over the reference span,
  `TANGLED=T`, full REF/ALT traversal strings, and no claim that the record is
  a primitive SNP/indel.
- Rust must not decompose one complex flubble into multiple SNP/indel records in
  the core emitter. vcfwave-like WFA decomposition is downstream.
- Subflubble families such as tiny, parallel, concealed, midi, or smothered are
  not separate VCF variant types in this spec. They may affect call selection
  only after they are mapped to a Lean-supported flubble source and covered by
  conformance.

## Hairpin and SUBR Records

Hairpin/SUBR output is graph-native inversion-like output, not a symbolic SV
encoding.

Rules:

- source must be a verified hairpin boundary;
- `VARTYPE=SUBR`;
- `REF` is the forward reference sequence over the supported region;
- `ALT` is the reverse-orientation alternate sequence supplied by the semantic
  call;
- `AT` contains the forward reference traversal followed by the reverse
  alternate traversal, for example `>1>2>3>4>5,<5<4<3<2<1`;
- `ES` and `LV` must be absent;
- no `<INV>`, `SVTYPE=INV`, `END`, or symbolic structural-variant fields are
  emitted in the first core subset.

This intentionally preserves povu's current `SUBR` convention and the Lean
semantic boundary. A downstream vcfwave-like or SV-oriented export may add an
inversion annotation, but it must not replace the graph-native `SUBR` semantics
without explicit provenance.

## Deterministic Ordering and Stable IDs

Records must be sorted by Lean's semantic order key:

1. selected contig/reference order;
2. `POS`;
3. source primary key;
4. source secondary key.

For flubbles, the source keys are the canonical open and close boundary edge
ids. For hairpins, they are the canonical outer and inner hairpin boundary edge
ids. The Rust extractor must preserve the same canonical orientation used to
form the semantic `VariantCall`.

This is an intentional improvement over current C++ generation order, which is
influenced by reference map order, RoV chunking, and SNE final-chunk behavior.
The current `two-ordered-substitutions` fixture is the minimal ordering guard,
but the Rust emitter should apply the order key to all records.

Default IDs are the current povu oriented source strings. They are stable across
runs when input graph ids and reference selection are stable. If a future mode
translates graph ids, the translation must happen before `ID`, `AT`, `ES`, and
any compatibility parent tag are generated so all provenance fields agree.

## Sample, Genotype, and Phase Handling

The first Rust contract supports `FORMAT=GT` only.

Sample columns:

- PanSN-style path metadata should group haplotypes into sample columns when
  available, with haplotype phases ordered deterministically by haplotype id and
  path order.
- If path names are not PanSN-compatible, Rust may use the full path name as a
  haploid sample column, matching existing documentation. This fallback must be
  explicit in metadata or logs.
- If no VCF-capable reference/path metadata are available, VCF emission must
  fail with a structured error instead of fabricating samples.

Genotypes:

- allele `0` is REF;
- alternate alleles are one-based indexes into ALT;
- `.` is missing for one phase;
- phases are serialized with `|`;
- haploid columns have one allele, for example `0` or `1`;
- if every phase in a column is missing, serialize a single `.`, matching
  current C++ compatibility;
- if only some phases are missing, keep them in place, for example `0|.`;
- `NS` counts sample columns with at least one non-missing phase;
- `AN` and `AC` come from semantic call counts and must agree with genotype
  semantics.

Unphased `/` genotypes, genotype likelihoods, depths, quality fields, and
sample-level filters are omitted from this contract.

## Compatibility Goals and Intentional Divergences

### vg deconstruct

Aligned behavior:

- VCF is reference-coordinate output selected before variant export.
- REF/ALT indels are anchored.
- `AT`-style traversal provenance is emitted by default.
- Nested levels are visible through `LV`.
- Missing genotype at nested sites is the default for haplotypes that do not
  traverse the child site.

Intentional divergences:

- Rust does not implement experimental `vg -n` fields such as `PA`, `PL`, `PR`,
  or `RL` in the first contract, because the predecessor study found release
  status differences.
- Rust does not emit off-reference nested child records by default. A nested
  child inside an insertion must be explicitly supported by a later raw graph
  compatibility mode or omitted/rejected with provenance.
- Rust does not implement `vg -R` star allele semantics in the first core
  subset.
- Rust does not implement `vg -u` untangled traversal coordinates or `-L`
  lossy traversal clustering in the core emitter.

### vcfbub

Aligned behavior:

- Raw graph records carry enough hierarchy and length information for a later
  parent-popping pass.
- A `popped` profile should use the vcfbub two-pass idea: drop large parents but
  rescue eligible children of popped parents.

Intentional divergence:

- Core Rust output uses current povu/Lean `ES` plus `LV`, not vg/vcfbub `PS` as
  a required field. A `vg-nesting-tags` compatibility mode may add `PS` as a
  derived parent-id alias once parent ids are available and fixture-confirmed.
  Until then, external `vcfbub` should not be assumed to consume default Rust
  VCF unchanged.

### vcfwave

Aligned behavior:

- Rust preserves `AT` provenance so a downstream decomposition pass can keep
  graph-origin information.
- Rust treats inversions/hairpins as a special case rather than blindly
  decomposing them into many primitive substitutions.

Intentional divergence:

- Core Rust output does not run WFA/BiWFA or split complex alleles. vcfwave-like
  decomposition is downstream and must update counts/genotypes/provenance
  explicitly.

### bcftools norm

Aligned behavior:

- Rust emits valid anchored alleles that are suitable for later reference-based
  normalization when a matching FASTA exists.

Intentional divergence:

- Rust does not left-align or normalize against a FASTA in core output.
  Reference-FASTA-dependent normalization is outside the flubble/hairpin proof
  boundary.

### Current C++ povu output

Preserved:

- `QUAL=60`, `FILTER=PASS`, `FORMAT=GT`;
- INFO order `AC`, `AF`, `AN`, `NS`, `AT`, `VARTYPE`, `TANGLED`, then `ES` and
  `LV` for non-`SUBR`;
- `SUBR` omits `ES`/`LV`;
- all-missing genotype columns serialize as `.`;
- current oriented source strings remain default IDs.

Changed intentionally:

- Rust should not duplicate the `GT` FORMAT header line.
- Rust should declare `AN` as `Integer`.
- Rust should render `AF` deterministically with useful precision, not one
  decimal place.
- Rust should use Lean semantic ordering, not chunk/ref-map generation order.
- Rust should expose structured API errors for unsupported input instead of
  process exits or successful empty VCFs.
- Rust output is semantically compatible with Lean first, not byte-compatible
  with C++ headers.

## Proposed Conformance Fixtures

These fixtures should be used by `vcf-modern-conformance-corpus` and the Rust
emitter task. They are semantic expectations, not byte-for-byte header tests.
Every record below also has `QUAL=60`, `FILTER=PASS`, and `FORMAT=GT`.

| Fixture | Expected semantic records |
| --- | --- |
| `minimal-substitution` | One flubble `SUB`: `CHROM=HG1#1#chr1`, `POS=2`, `ID=>0>3`, `REF=C`, `ALT=G`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>1,>2;VARTYPE=SUB;TANGLED=F;ES=>0>3;LV=0`, genotypes `HG1=0`, `HG2=1`. |
| `insertion-flubble` | One anchored `INS`: `CHROM=HG1#1#chr1`, `POS=1`, `ID=>0>1`, `REF=A`, `ALT=AG`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>0,>0>2;VARTYPE=INS;TANGLED=F;ES=>0>1;LV=0`, genotypes `HG1=0`, `HG2=1`. |
| `deletion-flubble` | One anchored `DEL`: `CHROM=HG1#1#chr1`, `POS=1`, `ID=>0>1`, `REF=AG`, `ALT=A`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>0>2,>0;VARTYPE=DEL;TANGLED=F;ES=>0>1;LV=0`, genotypes `HG1=0`, `HG2=1`. |
| `nested-deletion` | One nested inner `DEL`: `CHROM=HG1#1#chr1`, `POS=2`, `ID=>1>4`, `REF=CT`, `ALT=C`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>1>3,>1;VARTYPE=DEL;TANGLED=F;ES=>1>4;LV=1`, genotypes `HG1=0`, `HG2=1`, `HG3=.`. This protects default missing-for-nontraversal semantics. |
| `hairpin-inversion-subr` | One hairpin `SUBR`: `CHROM=ref`, `POS=2`, `ID=>1>5`, `REF=ACGTA`, `ALT=TACGT`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>1>2>3>4>5,<5<4<3<2<1;VARTYPE=SUBR;TANGLED=F`, genotypes `ref=0`, `alt=1`. `ES` and `LV` must be absent. |
| `linear-no-variant` | Header and sample columns for `HG1`, `HG2`, with no records. This protects successful empty callsets when reference/sample metadata exist and no flubbles are selected. |
| `two-ordered-substitutions` | Two records in raw order: first the `minimal-substitution` record at `POS=2`, then `CHROM=HG1#1#chr1`, `POS=4`, `ID=>3>6`, `REF=A`, `ALT=C`, `INFO AC=1;AF=0.5;AN=2;NS=2;AT=>4,>5;VARTYPE=SUB;TANGLED=F;ES=>3>6;LV=0`, genotypes `HG1=0`, `HG2=1`. This protects Lean order-key sorting. |

Additional fixtures recommended after the core emitter exists:

- multi-ALT same-type flubble with deterministic ALT and genotype indexing;
- non-terminating `AF`, for example one alternate allele among three called
  alleles, to prevent one-decimal precision loss;
- nested child inside an insertion, expected to be a structured unsupported
  VCF child unless an off-reference raw mode is explicitly enabled;
- all-missing diploid column serializing as a single `.`, plus partial missing
  phased genotype such as `0|.`;
- complex/tangled flubble emitted as graph-faithful `SUB` with `TANGLED=T`,
  then optionally decomposed by a downstream vcfwave profile.

## Migration Notes

For Rust users:

- `povu::gfa_to_vcf` and `GraphAnalysis::write_vcf` currently return
  not-implemented errors. After this spec is implemented, callers should expect
  real output only when the input graph has VCF-capable labels, references, and
  genotype metadata.
- Manual in-memory Rust graph construction cannot produce VCF until path and
  reference metadata are supported. Topology-only graphs can still be analyzed,
  but VCF output must fail cleanly.
- The Rust API should prefer returning a semantic document or structured error
  before writing bytes, so tests can inspect `VariantCall`/record semantics.

For users expecting old C++ output:

- Row semantics for the covered fixtures should match Lean/current C++
  expectations, including anchored alleles, `AT`, `ES`/`LV`, and `SUBR` without
  `ES`/`LV`.
- Header bytes will not be identical: duplicate `GT` metadata is removed,
  `AN` metadata is corrected, and `fileDate` may be absent in deterministic
  mode.
- Record order is defined by the semantic order key. If this differs from a
  historical chunking accident, Rust order wins.
- `AF` may be more precise than C++ one-decimal output.
- `vg -n`, `vg -R`, vcfbub parent popping, vcfwave decomposition, and
  `bcftools norm` are not silently baked into default Rust output. They are
  separate compatibility or normalization profiles.

For downstream WG tasks:

- `vcf-modern-rust-emitter` should implement the semantic types and validation
  checks before text formatting.
- `vcf-modern-conformance-corpus` should turn the fixture table above into
  Rust API tests, not only C++ CLI harness tests.
- `vcf-modern-untangling-design` should treat `raw-graph` records as input and
  specify popping/decomposition/normalization as explicit rewrites with
  provenance.
- `vcf-modern-synthesis-check` must continue to avoid structure-level claims
  until the canonical structure export comparison is complete.
