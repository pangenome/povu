# povu GFA-to-VCF pipeline inventory

Task: `lean4-povu-pipeline-inventory`

Scope: this is an inventory of the current C++ implementation only. No
implementation code was changed for this task.

## Source files and module roles inspected

### CLI and orchestration

- `app/main.cpp`: constructs `core::config`, dispatches subcommands after CLI
  parsing.
- `app/cli/cli.cpp`, `app/cli/cli.hpp`: defines the `gfa2vcf`, `decompose`,
  `call`, `info`, `prune`, and `vcf` subcommands and maps CLI flags onto
  `core::config`.
- `app/subcommand/gfa2vcf.cpp`, `app/subcommand/gfa2vcf.hpp`: implements the
  one-shot GFA-to-VCF flow by creating a temporary PVST forest directory,
  running `decompose`, then running `call`, then deleting the temporary
  directory.
- `app/subcommand/decompose.cpp`, `app/subcommand/decompose.hpp`: reads GFA,
  splits the graph into connected components, builds spanning trees, computes
  flubbles/subflubbles, and writes one `.pvst` file per nontrivial component.
- `app/subcommand/call.cpp`, `app/subcommand/call.hpp`: reads the graph, PVST
  forest, and reference selection, initializes VCF output, generates VCF record
  batches, and writes records.
- `include/povu/common/app.hpp`: owns `core::config`, including flags that
  affect graph loading, decomposition, reference selection, VCF output,
  streaming, and region filtering.

### IO and external parsing

- `src/mto/from_gfa.cpp`, `include/mto/from_gfa.hpp`: wraps `liteseq` GFA
  parsing and builds the internal bidirected graph.
- `src/mto/to_pvst.cpp`, `include/mto/to_pvst.hpp`: serializes PVST trees to
  `.pvst`.
- `src/mto/from_pvst.cpp`, `include/mto/from_pvst.hpp`: parses `.pvst` files
  back into `pvst::Tree`.
- `src/mto/to_vcf.cpp`, `include/mto/to_vcf.hpp`: emits VCF headers, contig
  lines, and records to stdout or per-sample files.
- `src/mto/common.cpp`, `include/mto/common.hpp`: file listing, text-file
  reading, and output-directory creation helpers.
- `cmake/deps.cmake`: pins `liteseq` at commit
  `3678c8b7d1c5e1bfaa9e3f07133e6e3f00567015`; the current GFA byte-level
  parser is therefore external to this repository.

### Graph, spanning tree, and PVST structures

- `src/povu/graph/bidirected.cpp`, `include/povu/graph/bidirected.hpp`:
  `bd::VG`, `bd::Vertex`, `bd::Edge`, reference metadata, vertex-to-reference
  step matrix, connected-component extraction, and graph summaries/output.
- `src/povu/refs/refs.cpp`, `include/povu/refs/refs.hpp`: `pr::Refs`,
  `pr::Ref`, PanSN/raw reference matching, sample/ploidy metadata, and genotype
  column metadata.
- `src/povu/graph/types.cpp`, `include/povu/graph/types.hpp`: graph-side enums
  for vertex ends, orientations, oriented IDs, walks, and reference steps.
- `src/povu/graph/spanning_tree.cpp`,
  `include/povu/graph/spanning_tree.hpp`: converts a bidirected connected
  component into a spanning tree with tree edges and backedges.
- `src/povu/graph/bracket_list.cpp`,
  `include/povu/graph/bracket_list.hpp`: bracket stack/list metadata used by
  simple-cycle equivalence and flubble classification.
- `src/povu/graph/tree_utils.cpp`, `include/povu/graph/tree_utils.hpp`:
  branch descriptors, Euler-tour metadata, LCA support, depth/pre/post maps,
  and backedge interval indexing used by subflubble algorithms.
- `include/povu/graph/pvst.hpp`: PVST vertex families, clans, route
  parameters, bounds, and `pvst::Tree`.
- `include/povu/common/constants.hpp`, `include/povu/common/core.hpp`: shared
  IDs, invalid sentinels, PVST symbols/version, VCF separators, and integer
  aliases.

### Flubble, hairpin, and subflubble computation

- `src/povu/algorithms/flubbles.cpp`,
  `include/povu/algorithms/flubbles.hpp`: simple-cycle equivalence,
  hairpin-boundary detection, equivalence-class stack construction, endpoint
  normalization, and base flubble PVST construction.
- `src/povu/algorithms/tiny.cpp`, `include/povu/algorithms/tiny.hpp`:
  identifies leaf flubbles that satisfy tiny-flubble trunk/branch criteria and
  retags them as `tiny`.
- `src/povu/algorithms/parallel.cpp`,
  `include/povu/algorithms/parallel.hpp`: identifies leaf flubbles that satisfy
  parallel/overlap criteria and retags them as `parallel`.
- `src/povu/algorithms/concealed.cpp`,
  `include/povu/algorithms/concealed.hpp`: computes concealed subflubbles,
  attaches them under flubbles, and records boundary/location metadata.
- `src/povu/algorithms/smothered.cpp`,
  `include/povu/algorithms/smothered.hpp`: computes smothered subflubbles under
  concealed vertices.
- `src/povu/algorithms/midi.cpp`, `include/povu/algorithms/midi.hpp`: computes
  midi bubbles from pairs of concealed children.

### Variant and VCF semantics

- `src/ita/variation/color.cpp`, `include/ita/variation/color.hpp`: chooses
  PVST vertices to call by checking whether all selected references traverse
  both boundary vertices.
- `src/ita/variation/rov.cpp`, `include/ita/variation/rov.hpp`: builds regions
  of variation (`ir::RoV`) from selected PVST vertices and optional genomic
  region filters.
- `src/ita/graph/graph.cpp`, `include/ita/graph/graph.hpp`: finds/sorts walks
  within an RoV using BFS first and a haplotype-lap sort fallback.
- `src/ita/variation/overlay.cpp`, `include/ita/variation/overlay.hpp`: builds
  depth matrices, compares reference/alternate haplotype rows, computes
  context-bounded minimal RoVs, records no-coverage/matches-reference data, and
  collects inversion seed pins.
- `src/ita/variation/sne.cpp`, `include/ita/variation/sne.hpp`: seed-and-extend
  support for inversion-like records from pin pairs.
- `src/ita/graph/slice_tree.cpp`, `include/ita/graph/slice_tree.hpp`: interval
  tree for inversion extension output.
- `src/ita/genomics/allele.hpp`: haplotype slices, allele traversal strings,
  REF/ALT DNA string extraction, variant context, and trek/minimal-RoV data
  structures.
- `src/ita/genomics/vcf.cpp`, `include/ita/genomics/vcf.hpp`: converts treks
  and inversion slice trees into `iv::VcfRec` objects and genotype fields.
- `src/ita/genomics/genomics.cpp`, `include/ita/genomics/genomics.hpp`: top
  level VCF-record producer for `call`; chunks RoVs and pushes `VcfRecIdx`
  batches to the bounded queue consumed by `mto::to_vcf`.

### Tests and existing docs

- `tests/main_tests.cc`: includes all current test sources into one GTest
  executable.
- `tests/integration_tests/pvst_tests.cc`: in-memory graph fixture for flubble
  and PVST hierarchy checks.
- `tests/unit_tests/spanning_tree_tests.cc`: in-memory graph fixture for
  spanning-tree construction.
- `tests/integration_tests/args_tests.cc`: placeholder GTest assertions only.
- `tests/data/LPA.gfa`: GFA fixture present in the repository, but not wired
  into the current GTest executable.
- `tests/CMakeLists.txt`: builds `test_povu` against `povulib`, `liteseq`,
  `fmt`, and `taywee::args`.
- `docs/README.md`: existing high-level user docs for `decompose`, `call`,
  PVST, and GFA input assumptions.
- `docs/vcf.md`: existing user-facing VCF semantics notes.
- `README.md`: repository build, test, and quick-start commands.

## CLI entry points and behavior-affecting options

### Global options

All subcommands share these flags from `app/cli/cli.cpp`:

- `--version`: prints `VERSION` and exits.
- `-v`, `--verbosity <int>`: sets `core::config::verbosity_`. It affects debug
  logging and whether `gfa2vcf` prints temporary-directory/decompose/call
  progress. With `--hairpins`, boundaries are printed to stderr independently
  of VCF records.
- `-t`, `--threads <int>`: sets `core::config::thread_count_`. `decompose`
  uses it to divide connected components among threads. `call` uses it for the
  `povu::thread::thread_pool` allocated in VCF generation.

No range validation is visible in the CLI layer for verbosity or thread count.

### `povu gfa2vcf`

`gfa2vcf` is the current end-to-end entry point:

```bash
povu gfa2vcf -i input.gfa -P SAMPLE_PREFIX
povu gfa2vcf -i input.gfa -r refs.txt
povu gfa2vcf -i input.gfa -P SAMPLE_PREFIX -o out_dir
```

Options affecting GFA-to-VCF behavior:

- `-i`, `--input-gfa <gfa>`: required GFA file path. The implementation passes
  a file path to `liteseq`; stdin is not a supported path in the current code.
- `-h`, `--hairpins`: enables hairpin-boundary reporting during decomposition.
  It prints `Boundary: <b1> <b2>` lines to stderr, but does not change the
  current PVST file format or VCF records.
- `-s`, `--subflubbles`: enables tiny, parallel, concealed, midi, and smothered
  PVST processing after base flubbles are found.
- `-c`, `--chunk-size <size>`: sets RoV batch size in VCF generation. Default
  is `100`.
- `-q`, `--queue-length <size>`: sets the bounded producer/consumer queue
  capacity. Default is `4`.
- Output destination:
  - `-o`, `--output-dir <dir>`: emits split VCF files under the directory,
    using selected sample labels as file bases.
  - If `-o` is absent, `gfa2vcf_handler` sets `stdout_vcf = true`, so a single
    combined VCF stream is written to stdout.
  - `--stdout` exists in the shared output group, but the handler only needs it
    when an explicit output-dir alternative is present; absence of `-o` already
    selects stdout.
- Reference source, exactly one source in the CLI group:
  - `-r`, `--prefix-list <path>`: file containing reference name prefixes, one
    per line.
  - `-P`, `--path-prefix <NAME>`: may be repeated; all paths matching the
    prefix/sample rule are selected.
  - positional `refs`: alternate CLI parameter source for reference prefixes.

`gfa2vcf` sets `inc_vtx_labels = true` and `inc_refs = true` before reading
the graph. It does not expose `call`'s `-f/--forest-dir` or
`-g/--restrict`; the temporary PVST forest and full-graph call are implicit.

### `povu decompose`

`decompose` is the first half of the two-step workflow:

```bash
povu decompose -i input.gfa -o pvst_dir
```

Options affecting decomposition:

- `-i`, `--input-gfa <gfa>`: required.
- `-o`, `--output-dir <dir>`: directory for `.pvst` files. Defaults to `.`.
- `-h`, `--hairpins`: print hairpin boundaries to stderr.
- `-s`, `--subflubbles`: compute and serialize subflubble annotations.
- Global `-t/--threads` affects component-level parallel decomposition.

Unlike `call`/`gfa2vcf`, the `decompose` CLI handler does not set
`inc_vtx_labels` or `inc_refs`. The decomposition algorithms do not need GFA
sequence labels or reference paths, so the graph can be read without them.

### `povu call`

`call` is the second half of the two-step workflow:

```bash
povu call -i input.gfa -f pvst_dir -P SAMPLE_PREFIX
povu call -i input.gfa -f pvst_dir -r refs.txt --stdout
povu call -i input.gfa -f pvst_dir -P SAMPLE_PREFIX -g ref:start-end
```

Options affecting VCF generation:

- `-i`, `--input-gfa <gfa>`: required.
- `-f`, `--forest-dir <dir>`: directory containing `.pvst` files. Defaults to
  `.`.
- `-g`, `--restrict <ref:start-end>`: optional region filter. The parser treats
  `start` as inclusive and `end` as exclusive, requires numeric coordinates,
  and looks up `ref` by exact reference tag via `g.get_ref_id`.
- `-c`, `--chunk-size <size>` and `-q`, `--queue-length <size>`: same streaming
  controls as `gfa2vcf`.
- Output and reference-source options are the same as `gfa2vcf`.

`call` also sets `inc_vtx_labels = true` and `inc_refs = true` before reading
the graph, because VCF construction needs sequence labels, reference walks,
contig metadata, and genotype columns.

## End-to-end dataflow: GFA bytes to VCF records

1. `app/main.cpp` calls `cli::cli`, which mutates `core::config`.
2. For `gfa2vcf`, `povu::subcommands::gfa2vcf::do_gfa2vcf` creates a temporary
   directory using `mkdtemp("/tmp/povu_gfa2vcf_XXXXXX")`.
3. `gfa2vcf` copies the app config, switches the copy to `decompose`, sets the
   temporary output directory, and calls `decompose::do_decompose`.
4. `decompose::do_decompose` calls `mto::from_gfa::to_bd`.
5. `mto::from_gfa::to_bd` builds a `liteseq::gfa_config_cpp` with the input
   file path, `inc_vtx_labels`, and `inc_refs`; `liteseq::gfa_new` reads and
   parses the GFA bytes.
6. `to_bd` creates `bd::VG` from the `liteseq::gfa_props` pointer, then:
   - iterates `gfa->vtx_arr_size`, skipping null vertices;
   - adds each vertex by numeric `v->id`;
   - stores `v->seq` as the vertex label only when `inc_vtx_labels` is true;
   - iterates link records and maps `liteseq` left/right sides to
     `pgt::v_end_e::l`/`pgt::v_end_e::r`;
   - adds references, vertex-to-reference step indices, and genotype metadata
     when `inc_refs` is true;
   - marks tips from missing left/right incident edges.
7. `decompose` calls `bd::VG::componetize` and deletes the original graph.
   Components with fewer than three vertices are skipped.
8. Each retained component is converted to `pst::Tree::from_bd`:
   - each bidirected graph vertex can become two tree vertices, one for each
     side;
   - if graph tips exist, a dummy root is added and tips get simplifying
     root behavior;
   - internal left/right sides are joined by black tree edges;
   - graph adjacency contributes gray tree edges or backedges depending on DFS
     discovery state;
   - self-loops and tips have explicit backedge handling.
9. `povu::flubbles::find_flubbles` mutates the spanning tree:
   - `simple_cycle_equiv` computes `hi` values, bracket lists, capping and
     simplifying backedges, and equivalence classes;
   - hairpin boundaries are collected during this pass and optionally printed;
   - `compute_eq_class_stack` and `compute_eq_class_metadata` derive endpoint
     pairs for flubbles;
   - a `pvst::Tree` is created with a dummy root and `pvst::Flubble` vertices.
10. If `--subflubbles` is set, `decompose_component` computes `ptu::tree_meta`
    and runs, in order, `find_tiny`, `find_parallel`, `find_concealed`,
    `find_midi`, and `find_smothered`.
11. `mto::to_pvst::write_pvst` writes one `<component_id>.pvst` file per
    processed component.
12. `gfa2vcf` copies the original app config again, switches the copy to
    `call`, sets `forest_dir` to the temporary PVST directory, and calls
    `call::do_call`.
13. `call::do_call` starts three reads:
    - graph read with `mto::from_gfa::to_bd(app_config)`;
    - PVST forest read with `read_pvsts`, which reads all `.pvst` files in
      `forest_dir` and computes heights;
    - reference-prefix file read if `-r/--prefix-list` was used.
14. After graph and refs are available, `call` resolves each requested prefix
    with `g->get_refs_in_sample(prefix)`, builds `sample_to_ref_ids`, and
    builds the union `vcf_ref_ids`.
15. `mto::to_vcf::VcfOutput` is constructed for stdout or split files, and
    `init_vcfs` writes common VCF headers, contig lines for selected prefixes,
    and a genotype column header.
16. `call` starts a producer thread running `ig::gen_vcf_rec_map`; the main
    thread consumes `iv::VcfRecIdx` batches from a bounded queue and writes
    them via `mto::to_vcf::write_vcfs`.
17. `ig::gen_vcf_rec_map` optionally parses `-g/--restrict`, creates all RoVs
    with `ir::gen_rov`, processes them in chunks, and for each chunk:
    - calls `po::overlay_generic` for each RoV to build context-bounded treks;
    - may run `ise::sne` to convert inversion seed pins into context-free
      inversion slice trees;
    - calls `iv::gen_vcf_records` to produce VCF records grouped by reference
      haplotype.
18. `mto::to_vcf::write_vcfs` writes each record with:
    - `CHROM` from the reference tag;
    - `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, and
      genotype columns from `iv::VcfRec`;
    - stdout or split stream selection from `VcfOutput`.
19. `gfa2vcf` deletes the temporary PVST directory with `fs::remove_all`.

## Core data structures and module boundaries

- `core::config`: single mutable configuration object shared across CLI and
  subcommands. It carries task type, GFA path, decomposition flags, graph read
  flags, streaming parameters, output mode, reference-selection source,
  reference prefixes, optional region filter, and VCF analysis options.
- `bd::VG`: in-memory bidirected graph. It owns:
  - `Vertex` objects keyed by numeric GFA segment ID and indexed internally;
  - `Edge` objects keyed by internal indices and incident on left/right vertex
    ends;
  - `tips_` as side/id pairs;
  - `vertex_to_step_matrix_`, indexed by graph vertex index then reference ID,
    containing sorted step indices;
  - `pr::Refs`, wrapping `liteseq` reference pointers and genotype metadata;
  - the `liteseq::gfa_props *`, freed in the destructor.
- `pr::Refs`: wraps path/reference records. PanSN refs populate sample ploidy
  metadata by hap ID; raw refs fall back to undefined ploidy. Genotype column
  names and hap-to-column/phase mappings are generated once after references
  are loaded.
- `pst::Tree`: spanning tree over vertex sides rather than graph vertices.
  Tree vertices store graph vertex IDs, side type, DFS/pre/post numbers,
  parent/child edges, incoming/outgoing backedges, and `hi`. Tree edges carry
  black/gray color and equivalence class. Backedges carry source/target,
  type, and equivalence class.
- `pvst::Tree`: tree of variation structures. Vertices are polymorphic
  `VertexBase` values with families: `dummy`, `flubble`, `tiny`, `parallel`,
  `concealed`, `smothered`, `midi`, and `sne_exp`. Route-capable vertices
  return route parameters with oriented start/end IDs and route direction.
- `ptu::tree_meta`: derived metadata used only for subflubble detection:
  Euler tour arrays, first occurrences, depth, pre/post maps, backedge lists,
  offset tables, LoA/HiD values, and LCA support.
- `ir::RoV`: variant-region wrapper around a PVST vertex plus a sorted list of
  graph vertex IDs inside the region.
- `ia::hap_slice`: a slice of a `liteseq::ref_walk`, carrying ref ID, start
  step, and length. It is the source of REF/ALT traversal strings, DNA strings,
  and POS calculation.
- `ia::minimal_rov`: one context-bounded reference allele with matching-ref
  haplotypes and grouped alternate allele slices.
- `ia::trek`: all minimal-RoV calls for one RoV, keyed by reference haplotype
  and context; also tracks no-coverage and matches-reference haps.
- `iv::VcfRec`: a single VCF record before serialization. It stores reference
  haplotype, position, ID, enclosing flubble, reference slice, alternate slice
  groups, genotype columns, NS, PVST height, variant type, and tangled flag.
- `iv::VcfRecIdx`: batch container mapping reference haplotype IDs to vectors
  of `VcfRec`.
- `mto::to_vcf::VcfOutput`: output abstraction for combined stdout or split
  per-sample files.

## GFA parsing and normalization rules

- The repository docs state current GFA input support as GFA 1.0 and require
  unique numeric segment names. The C++ graph types use `pt::id_t`/`pt::idx_t`
  unsigned 32-bit aliases, and `bd::VG` maps numeric vertex IDs to internal
  indices.
- Actual GFA byte parsing is delegated to `liteseq::gfa_new`; this repository
  only passes `file path`, `include vertex labels`, and `include references`.
- Sequence labels are loaded into `bd::Vertex::label_` only when
  `inc_vtx_labels` is true. `call` and `gfa2vcf` set this flag; `decompose`
  does not.
- Reference/path metadata is loaded only when `inc_refs` is true. `call` and
  `gfa2vcf` set this flag; `decompose` does not.
- Link sides from `liteseq` are mapped directly:
  - `liteseq::vtx_side_e::LEFT` -> `pgt::v_end_e::l`;
  - any other side -> `pgt::v_end_e::r`.
- For every reference step, `to_bd` records `v_id -> ref_idx -> step_idx` in
  `vertex_to_step_matrix_`, inserting step indices in sorted order and
  avoiding duplicates.
- Reference selection uses `pr::Refs::get_refs_in_sample`:
  - PanSN references match when `liteseq::get_sample_name(r)` equals the
    requested sample/prefix string;
  - raw references match when the requested string is a prefix of the raw
    sample name;
  - all references also get a fallback prefix match against the full reference
    tag.
- Genotype metadata groups references sharing a sample. For PanSN refs, hap ID
  determines phase column order. For raw/non-PanSN refs, ploidy is undefined
  and a single genotype entry is used for each sample.
- Tips are inferred after vertices/edges/references are loaded:
  - isolated vertices are reported as a left-side tip only;
  - vertices missing only left edges are left tips;
  - vertices missing only right edges are right tips.
- `bd::VG::componetize` copies graph vertices, labels, edges, and tips into
  connected components. It does not copy reference metadata, which is acceptable
  for decomposition because PVST construction uses topology, not references.
- Components with fewer than three vertices are skipped during decomposition.
- There is no visible local validation of GFA overlap strings or non-numeric
  segment names in `mto::from_gfa.cpp`; those behaviors are inherited from
  `liteseq`.

## Flubble, hairpin, and PVST output semantics

- `pst::Tree::from_bd` represents each bidirected vertex side separately. A
  graph vertex visited from one side creates two adjacent tree vertices: the
  visited side and its complement. The side-pair edge is black; graph traversal
  edges that discover new vertices become gray tree edges; already-seen
  non-parent relationships become backedges.
- If graph tips exist, a dummy root is inserted. Tip-like vertices can receive
  simplifying backedges to the root during simple-cycle equivalence handling.
- `simple_cycle_equiv` computes `hi`, bracket lists, equivalence classes, and
  additional capping/simplifying backedges. Tree edge classes are assigned from
  the top bracket's recent class on retreat.
- Hairpin boundaries are side effects of `simple_cycle_equiv`. When a bracket
  list becomes empty on a non-dummy tree vertex, the current boundary starts at
  that graph vertex ID. If a simplifying backedge remains on top while inside a
  hairpin, the boundary end is extended. Boundaries are pushed when a leaf or
  root ends the hairpin.
- `--hairpins` prints collected boundaries as `Boundary: <b1> <b2>` to stderr.
  Hairpin boundaries are not serialized into `.pvst` and are not present in
  VCF `INFO`.
- Base flubbles are generated from pairs of equivalent-class stack entries.
  `compute_ai_zi` derives the spanning-tree interior bounds from the tree edges
  associated with the pair. `gen_fl` normalizes endpoints: if both endpoints
  are reverse oriented, it reverses the endpoint order and makes both forward.
  Other orientation combinations are preserved.
- The PVST always starts with a dummy root vertex. Base flubbles are added as
  nested children according to repeated equivalence-class stack encounters.
- With `--subflubbles`, leaf flubbles may be retagged:
  - `tiny` for selected short trunk/branch cases;
  - `parallel` for selected branch/trunk overlap cases.
- `concealed`, `midi`, and `smothered` vertices are added under existing
  flubbles/concealed vertices. They carry route parameters and bounds needed by
  later RoV and VCF generation.
- `.pvst` serialization writes five tab-separated columns:
  1. record symbol;
  2. file vertex index;
  3. vertex label from `as_str()`;
  4. comma-separated child indices or `.`;
  5. route, `L` for start-to-end, `R` for end-to-start, or `.`.
- PVST symbols and families:
  - `H`: header, version `0.0.3`;
  - `D`: dummy;
  - `F`: generic flubble;
  - `T`: tiny;
  - `O`: parallel/overlap;
  - `C`: concealed;
  - `S`: smothered;
  - `M`: midi.
- `from_pvst` only reconstructs route-visible data from serialized labels and
  routes. Some richer in-memory metadata, such as computed `ai`, `zi`, `m`,
  `n`, and concealed/smothered bounds, is not fully represented after PVST
  round-trip; downstream VCF generation uses route parameters, height, family,
  and tree structure.

## Variant/RoV computation semantics

- `call` chooses reference haplotypes from CLI prefixes and builds
  `vcf_ref_ids`, the set of reference IDs to call against.
- `ir::gen_rov` processes each PVST:
  - `ic::color_pvst` starts at the PVST root and skips subflubble-clan
    vertices as direct call targets;
  - a route-capable PVST vertex is selected when all requested reference IDs
    have at least one step at both boundary graph vertices;
  - if both a parent and child are selected, the parent is removed so more
    specific selected children win;
  - `find_walks` is run for each selected PVST vertex.
- `find_walks` first builds and sorts a BFS tree between oriented route
  endpoints. If BFS sorting fails, it falls back to `gen_sort`, which derives a
  sorted vertex list from haplotype laps between the route endpoints.
- RoVs smaller than three sorted vertices are skipped with a warning.
- `call --restrict ref:start-end` only includes PVST vertices whose route
  endpoint vertices have loci inside the requested half-open interval on the
  exact reference tag. It does not currently scan every internal vertex in the
  RoV for overlap.
- `overlay_generic` builds a depth matrix with one row per haplotype/reference
  and one column per sorted RoV vertex. It compares each selected reference row
  to every other haplotype row:
  - empty alternate rows become no-coverage entries;
  - rows with no variant contexts and no inversion pattern become
    matches-reference entries;
  - bounded differences become `minimal_rov` contexts with reference and
    alternate `hap_slice` values;
  - local inversions add pin pairs to the shared `pin_cushion` instead of
    immediately producing context-bounded records;
  - tangled matrices are untangled into multiple matrices before trek creation.
- `iv::gen_vcf_records` emits context-bounded records from treks and
  context-free inversion records from SNE slice trees. Variant types are:
  - `INS`: alternate slice length greater than reference slice length;
  - `DEL`: alternate slice length less than reference slice length;
  - `SUB`: alternate slice length equal to reference slice length;
  - `SUBR`: reverse-substitution/inversion-like records from SNE.

## VCF record construction rules

### Header and output routing

- `mto::to_vcf` emits VCF version `4.2`, `##fileDate=<today>`,
  `##source=povu`, `GT` format metadata, and INFO metadata for `AC`, `AT`,
  `AN`, `AF`, `NS`, `VARTYPE`, `TANGLED`, and `LV`.
- The current header writes `##FORMAT=<ID=GT,...>` twice.
- `init_vcfs` writes contig lines only for references matched by the selected
  sample/prefix strings.
- The final `#CHROM` line includes all genotype column names in the graph,
  not only the selected references.
- Stdout mode writes all headers and records to one stream. Split-file mode
  creates `<sample_prefix>.vcf` under `--output-dir` and routes all selected
  reference IDs for that sample/prefix to that file.

### Per-record fields

- `CHROM`: `g.get_ref_by_id(ref_id).tag()`.
- `POS`: computed from `ia::hap_slice::comp_pos`:
  - `DEL`/`INS`: `loci[ref_start_idx + 1] - 1`;
  - `SUB`/`SUBR`: `loci[ref_start_idx + 1]`.
- `ID`:
  - context-bounded records use the oriented context boundary string, for
    example `>u>v`;
  - SNE/SUBR records construct an oriented start/end ID from the reference
    slice ends and forwardize the ID when both ends are reverse.
- `REF`:
  - `INS`/`DEL` include the first base of the left anchor step, then internal
    slice bases after dropping the terminal context step;
  - `SUB` uses the slice DNA string after dropping both context endpoints;
  - `SUBR` uses the full slice.
- `ALT`: one comma-separated value per unique alternate walk group, using the
  same DNA extraction rules as `REF`.
- `QUAL`: fixed string `60`.
- `FILTER`: fixed string `PASS`.
- `FORMAT`: fixed string `GT`.
- Genotypes:
  - `0` denotes the reference allele;
  - alternate allele numbers start at `1`;
  - phased PanSN samples are joined with `|`;
  - columns where every phase is missing are serialized as `.`;
  - otherwise missing phases remain `.` inside the phased genotype.

### INFO fields

- `AC`: count of haplotype slices in each alternate allele group; reference
  alleles are not included.
- `AF`: `AC / AN`, formatted with one decimal digit.
- `AN`: count of reference-allele haplotypes plus all alternate haplotype
  slices for the record.
- `NS`: number of sample genotype columns with at least one non-missing allele.
- `AT`: reference allele traversal plus alternate traversals, using oriented
  graph vertex IDs. For insertions/deletions, terminal context is omitted from
  the traversal string according to `hap_slice::as_str`.
- `VARTYPE`: one of `DEL`, `INS`, `SUB`, `SUBR`.
- `TANGLED`: `T` or `F`.
- `ES`: enclosing PVST vertex string for non-`SUBR` records only.
- `LV`: `pvst_height - 1` for non-`SUBR` records only. A top-level child of
  the dummy root has `LV=0`.

## Edge cases, implicit behavior, and undefined areas

- GFA parsing details such as malformed lines, duplicate segment names,
  non-numeric segment names, and overlap-string interpretation are delegated to
  `liteseq`; they are not specified or tested in this repository.
- `decompose` reads GFA without sequence labels and references, while `call`
  and `gfa2vcf` read with both enabled. Proof interfaces should distinguish
  topology-only graph construction from VCF-capable graph construction.
- Isolated graph vertices are recorded as left tips only, not both left and
  right tips.
- `bd::VG::componetize` starts DFS at vertex index `0`; behavior for an empty
  graph is not guarded locally.
- Decomposition skips components with fewer than three vertices. No empty PVST
  placeholder is written for skipped components.
- `thread_count(app_config, components.size())` divides by configured thread
  count. There is no visible guard against `--threads 0`.
- `read_pvsts` exits if no `.pvst` files are found in the forest directory.
- Reference prefix selection relies partly on prefix matching against raw names
  and tags. Overlapping prefixes can map the same reference to multiple output
  streams during header setup, but `VcfOutput::ref_id_to_ofs_idx_` stores one
  stream per reference ID.
- If a requested reference prefix matches no graph references, debug builds
  assert, but release builds have no obvious fatal check before VCF
  initialization and RoV coloring.
- `ic::has_any_refs` currently has the same all-references behavior as
  `has_all_refs`, but it is not used by the active `color_pvst` implementation.
- `color_pvst` skips subflubble-clan vertices as direct call targets. It can
  still call flubble-like vertices that were retagged as `tiny` or `parallel`.
- `from_pvst` reconstructs only serialized route data for PVST vertices. Bounds
  and some derived fields become invalid sentinels after reading `.pvst`.
- `call --restrict` uses only route endpoint loci to decide region overlap, so
  a variant whose internal vertices overlap a region but whose endpoints do not
  will be skipped.
- RoV walk enumeration has pressure valves: `MAX_FLUBBLE_STEPS = 1000` and
  `MAX_UNBLOCK_CTR = 10000`. Exceeding these can silently prune work with
  warnings/continues rather than a proof-level result.
- VCF records are emitted in generation order by chunk and reference map
  iteration; there is no final sort by genomic position.
- `AF` is formatted to one decimal place, so allele frequencies are rounded
  aggressively.
- The VCF header duplicates the `GT` FORMAT metadata line.
- For `SUBR` records, `ES` and `LV` are intentionally omitted by the writer.
- Current SNE/inversion generation appears gated by `if (chunk_num ==
  CHUNK_SIZE)` in `gen_vcf_rec_map`. With the default chunk size of `100`, this
  means SNE is only invoked on chunk number 100 rather than after each chunk or
  at the end. This looks like a behavior mismatch for inversion/SUBR output and
  should be fixed in a follow-up task, not in this inventory.

## Current tests, fixtures, and coverage gaps

### Existing test commands

The repository build/test path is CMake, not root Cargo:

```bash
cmake -B build -DPOVU_ENABLE_TESTING=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build -- -j 3
ctest --test-dir build
```

Validation run during this inventory:

- `cmake -B build -DPOVU_ENABLE_TESTING=ON -DCMAKE_BUILD_TYPE=Debug`: passed.
- `cmake --build build -- -j 3`: passed and built `bin/povu`.
- `ctest --test-dir build`: passed, 5/5 tests.
- Documentation scan for non-ASCII text in this file: no matches found.

No command failures were encountered while validating this inventory.

The dependency task `lean4-proof-quality-pass` already noted that root
`cargo build` is not applicable because this repository root has no
`Cargo.toml`; it also noted a pre-existing `povu-rs` FFI header-path failure
and created follow-up task `fix-povu-rs`.

### Existing tests and fixtures

- `PVSTTest.VertexCount`, `PVSTTest.HasVertices`, and
  `PVSTTest.VertexHierarchy` exercise in-memory graph -> connected component
  -> spanning tree -> base flubble PVST generation. Expected PVST labels are
  `.`, `>1>7`, and `>4>6`.
- `SpanningTreeTest.ConstructSpanningTree` exercises in-memory bidirected graph
  -> spanning tree construction without assertions on tree shape.
- `HelloTest.BasicAssertions` is a placeholder assertion test.
- `tests/data/LPA.gfa` is a substantial GFA fixture with GFA header, segments,
  and links. It is not currently used by `test_povu`.

### Missing coverage relevant to Lean proof and conformance

- No current test reads `tests/data/LPA.gfa` through `mto::from_gfa::to_bd`.
- No test covers malformed GFA, duplicate/non-numeric segment names, link side
  normalization, isolated vertices, self-loops, multiple components, or skipped
  small components.
- No test asserts `--hairpins` stderr output or absence/presence of hairpin
  data in PVST/VCF artifacts.
- No tests cover `--subflubbles`, tiny/parallel retagging, concealed,
  smothered, or midi PVST serialization/round-trip behavior.
- No test exercises `gfa2vcf` end-to-end, either to stdout or split output.
- No test exercises `call` reading `.pvst` files from disk.
- No VCF test asserts header fields, duplicate `GT` metadata, contig selection,
  genotype columns, `INFO` fields, REF/ALT anchoring, `AT`, `ES`, `LV`, `SUBR`,
  or output ordering.
- No test covers PanSN genotype grouping, raw-reference fallback behavior,
  reference-prefix ambiguity, missing prefix behavior, or selected-vs-all
  genotype columns.
- No test covers `--restrict` parsing or endpoint-only overlap semantics.
- No test covers chunking/queue behavior, including the current SNE gating
  condition.
- No test covers release-vs-debug behavior for assertions on missing refs.

## Follow-up candidates identified

- Fix or specify SNE/SUBR generation gating in `ita::genomics::gen_vcf_rec_map`
  so inversion-like records are generated at the intended time.
- Add a conformance fixture that drives `povu gfa2vcf` from a small GFA with
  paths to VCF and checks headers, REF/ALT anchoring, genotypes, and INFO
  fields.
- Add parser-boundary tests for numeric segment IDs, duplicate IDs, link-side
  orientation, tips, isolated vertices, and self-loops.
- Add PVST round-trip tests for every serialized family symbol and route
  direction.
- Specify and test behavior for missing or ambiguous reference prefixes.
