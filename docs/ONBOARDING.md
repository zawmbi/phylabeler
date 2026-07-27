# PhyLabeler Onboarding Guide

## What This Project Does

PhyLabeler checks whether clades in a phylogenetic tree are monophyletic with
respect to the NCBI taxonomy database.

At a high level:

1. Load taxonomy data from `names.dmp` and `nodes.dmp`, from an un-extracted
   NCBI archive, or from a cached pickle.
2. Parse a Newick tree into an in-memory node structure.
3. Resolve every tip label in the tree to an NCBI taxonomic identifier.
4. For each internal non-root clade, compute the taxonomic MRCA shared by all
   tips in that clade.
5. Test whether that MRCA is exclusive to the clade or also appears in taxa
   outside the clade.
6. Mark the clade as monophyletic, non-monophyletic, or unresolved.
7. Write a labeled Newick tree and a text report.

## Fast Reading Order

Use this order if you want to understand the modern codebase quickly:

1. `main.py`
   Learn how the project decides between CLI and GUI and where the main
   workflows start.
2. `archive_io.py`
   Understand how compressed and archived inputs are detected and streamed.
   Everything else opens files through this module.
3. `taxonomy_db.py`
   Understand how NCBI taxonomy is loaded, parsed, cached, and queried.
4. `tree_parser.py`
   Understand the internal tree model and how bipartitions are derived.
5. `monophyly.py`
   Understand the analysis algorithm and the result objects.
6. `gui.py`
   See how the analysis engine is exposed through the Tkinter interface.

Only after that, look at the legacy files:

- `LabelPhy.py`
- `TreeCode.py`
- `Renamer.py`
- `NameNavigate.py`

These older modules are useful for project history, but they are not the main
implementation you should study first.

## Core Modules

### `main.py`

This is the project dispatcher.

- Parses command-line arguments.
- Chooses GUI or CLI mode.
- Loads taxonomy.
- Parses one tree or many trees.
- Runs the monophyly checker.
- Writes labeled tree and report outputs.

Important functions:

- `main()`
- `run_cli(args)`
- `_process_single(args, checker)`
- `_process_batch(args, checker)`
- `run_gui()`

## `archive_io.py`

This file is the single entry point for reading input files, so no other
module has to care whether something is compressed or packed in an archive.

### Main responsibilities

- Detect gzip/bz2/xz compression and tar/zip containers from magic bytes,
  never from the file extension
- `open_text()` for single files, decompressing transparently
- `stream_members()` to read archive members without unpacking them to disk
- `extract_matching()` when a caller genuinely needs files on disk, writing
  members flat under their basename so archive paths cannot escape

### Why magic bytes

Users rename downloads and strip suffixes. A `taxdump.tar.gz` saved as
`taxonomy.gz` is still a gzipped tar, and this module treats it as one.

## `taxonomy_db.py`

This file builds the project knowledge base from NCBI taxonomy.

### Main responsibilities

- Parse `names.dmp` and `nodes.dmp`, plain or gzipped
- Parse both straight out of an un-extracted `taxdump.tar.gz` / `taxdmp.zip`
  in a single streaming pass, via `load_archive()`
- Store lookup dictionaries
- Save and load a pickle cache
- Download fresh taxonomy from NCBI
- Provide lineage lookups used by the monophyly engine

### Key data structures

- `name_to_code`
  Maps taxon name to NCBI tax ID.
- `code_to_names`
  Maps tax ID to all known names.
- `code_to_scientific`
  Maps tax ID to the scientific name.
- `code_to_parent`
  Maps tax ID to parent tax ID.
- `code_to_rank`
  Maps tax ID to taxonomic rank.

### Tip-name resolution behavior

`name_to_code` is keyed by a normalised form (lowercase, underscores as
spaces, whitespace collapsed), so every lookup is O(1). `resolve_tip()` then
tries, in order:

1. the whole label
2. the label with any `|`, `/`, `@`, `#` suffix stripped
3. progressively fewer trailing tokens, down to the genus

Each match records *how* it matched, because a genus-level fallback is a much
weaker claim than a species-level hit and the report says so.

### Why scientific names take priority

names.dmp holds ~4.9M rows, only 2.9M of which are scientific names; the rest
are synonyms, authorities and strain designations. Letting the last row win
made `"Bacteria"` resolve to a stick-insect genus rather than the domain.
Scientific names now win, and genuine homonyms are collected in `homonyms`
so the UI can flag them instead of quietly picking one.

## `tree_parser.py`

This file defines the in-memory tree representation and the Newick parser.

### `Node`

Each node stores:

- `label`
- `branch_length`
- `parent`
- `children`
- `is_tip`
- `monophyly_status`
- `taxonomy_code`
- `taxonomy_name`

### Important tree methods

- `get_tip_labels()`
  Collect all descendant tip labels.
- `get_tip_nodes()`
  Collect all descendant tip node objects.
- `get_bipartitions()`
  Return a post-order list of internal non-root clades.
- `to_newick()`
  Rebuild the tree as a Newick string.
- `count_tips()`, `count_internal()`, `depth()`
  Provide structural summaries.

### Why bipartitions matter here

The monophyly checker analyzes each internal clade independently. In practice,
the clade is represented by the set of tip labels below a node.

## `monophyly.py`

This is the project's analysis engine.

### `MonophylyResult`

Each analyzed clade gets a result object containing:

- the node reference
- the tip labels in the clade
- the taxonomic MRCA code
- the MRCA scientific name
- the MRCA rank
- the status
- the candidate names for that MRCA
- the intruder taxa outside the clade

### Two views of the same data

`analyze()` returns both, and they answer different questions.

**Taxon-level (`check_taxa`)** is the report a biologist reads. It runs one
bottom-up pass (`_taxon_pass`) that scores *every* taxon at *every* node at
once, carrying two things per node: how many tips below it matched the
taxonomy, and how many of those belong to each taxon. From that fall out:

1. **clusters** — the maximal clades containing only members. A child is one
   when the child is pure and its parent is not. Exactly one cluster means
   monophyletic.
2. **home clade** — adopting node `v` as "the group" costs
   `foreign_inside + members_outside`, which equals
   `M + resolved[v] - 2 * members[v]`. `M` is fixed per taxon, so minimising
   `resolved[v] - 2 * members[v]` finds the cheapest clade in the whole tree.
3. **intruders / outliers / unknown_inside** — the tips inside the home clade
   that don't belong, the members that sit outside it, and the tips that
   matched nothing.

`conflict_cost` is `intruders + outliers`: the fewest tips that must move.

### Why the cost is minimised over every node

The obvious shortcut — `min(intruders under the MRCA, members outside the
largest cluster)` — evaluates the same objective at only two candidate nodes,
so it is an upper bound. On the bundled test tree it rates **rosids at 59**
when evicting one tip and gathering one other makes it a clade: the true cost
is **2**. `test_monophyly.py::test_conflict_costs` pins this.

### Fragmentation is not severity

`cluster_count` says how scattered the members are, not how wrong the tree is.
A single intruder nested deep in a ladder splits Fabaceae's 25 members into 8
pure clades while remaining a one-tip fix. Never sort on it.

### Unplaced tips are wildcards

Tips that matched no taxon are counted in neither total, so they never spoil a
cluster and never become intruders — otherwise one unmatched sequence would
appear to break every group it sits inside.

**Node-level (`_check_nodes`)** is what colours the drawing and supplies the
taxonomy labels written into the output Newick. For each internal node it finds
the taxonomic MRCA of the node's tips and checks that no tip outside the node
shares it.

### Why the node-level view is not the report

It produces one row per internal node, so a single conflict reappears at every
nested node beneath it — a 232-tip tree gave four consecutive
`Mesangiospermae` rows describing one problem. Groups, not nodes, are what the
taxonomy names.

### Important conceptual point

Monophyly here is implemented as an exclusivity test relative to the other tips
present in the input tree, not as a separate phylogenetic inference procedure.

## `gui.py`

This file wraps the core engine in a Tkinter desktop interface.

### Main UI sections

- Taxonomy loading controls
- Tree file list
- Tree visualization canvas
- Analysis buttons
- Results notebook

### GUI workflow

1. Start the app.
2. Load taxonomy from cache, local files, or NCBI download.
3. Load one tree or many trees.
4. Parse the selected tree.
5. Run the monophyly analysis.
6. View the results:
   - Summary
   - Non-monophyletic clades
   - Monophyletic clades
   - Unresolved tips
7. Export a labeled tree or report.

### Tree drawing model

`_layout_tree()` computes coordinates; `_draw_tree()` renders them.

- y-position comes from tip order
- x-position comes from node depth (cladogram, tips aligned in a column) or
  from cumulative branch length when "Branch lengths" is ticked (phylogram)
- branches are coloured by monophyly status
- clade names are drawn to the **right** of their node on an opaque chip;
  drawing them to the left clipped every shallow clade off the canvas
- support values are drawn small and grey, in their own visual channel
- tip colour encodes match quality: grey unresolved, amber matched loosely,
  red when the tip intrudes into the selected group

Selecting a row in the results tables highlights that clade and scrolls to it.

## CLI Workflow

### Single-tree mode

Typical path:

1. load taxonomy
2. parse one tree
3. check monophyly across all internal clades
4. label monophyletic internal nodes
5. write:
   - `*.labeled.nwk`
   - `*.report.txt`

### Batch mode

Typical path:

1. scan a directory for Newick-like files
2. parse each tree
3. run the same checker on each tree
4. write one labeled output per tree
5. write a `batch_summary.txt`

## Output Files

### Labeled Newick

The labeled tree preserves the original topology and branch lengths, but
monophyletic internal nodes receive cleaned taxonomy labels suitable for Newick.

### Report

The report summarizes:

- total internal nodes analyzed
- monophyletic count
- non-monophyletic count
- unresolved count
- a sorted list of monophyletic clades
- a sorted list of non-monophyletic clades with intruders
- unresolved tips, if any

## Recommended Study Path For PI Meetings

Use this sequence when explaining the project:

1. Problem statement
   This tool checks whether gene-tree clades agree with NCBI taxonomy.
2. Inputs
   Newick tree plus NCBI taxonomy files or cache.
3. Internal model
   Parsed tree plus taxonomy lookup tables.
4. Core method
   MRCA sharing inside the clade, then exclusivity against the complement.
5. Outputs
   Labeled tree and report.
6. Limitations
   Depends heavily on tip-name matching and the taxa present in the tree.

## Known Caveats

- Tip matching drops trailing name tokens but does no fuzzy/edit-distance
  matching; misspellings stay unresolved.
- `conflict_cost` counts tips, not tree edits. It is the minimum number of
  tips that must move, which is the right work-queue ordering but is not a
  branch-length or likelihood-aware measure of how wrong the tree is.
- A star tree reports its groups as non-monophyletic even though a polytomy
  does not actually contradict them — it is unresolved, not in conflict. The
  distinction is not currently made.
- Groups with a single sampled tip are omitted: they are monophyletic by
  definition and say nothing.
- The modern UI and CLI share the same core analysis engine.
- The legacy files reflect an older implementation style and should not be
  treated as the primary architecture.

## Tests

Both suites are self-contained, need no data files and no test runner:

```
python test_monophyly.py     # resolution, parsing, monophyly
python test_archive_io.py    # compressed and archived inputs
```

## Practical Questions To Ask While Reading

- How are tip labels normalized before lookup?
- What exactly counts as a bipartition in this implementation?
- How is the MRCA computed from lineage lists?
- Why does a clade fail monophyly even if its internal members share a common
  taxonomy?
- Which parts are UI-only and which parts are analysis-core?
- What assumptions are made about the input tree labels?

## Suggested Next Improvements

If you want to extend the project later, likely high-value directions are:

- stronger tip-name normalization and alias handling
- search/filter tools in the results UI
- richer tree rendering for large trees
- explicit export of unresolved-name diagnostics
- clearer separation of legacy and modern code paths
