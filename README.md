# phylabeler
Program To label phylogenetic trees

Necessary files for the program can be downloaded here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

## Running

```
python main.py                    # GUI
python main.py --tree FILE ...    # CLI (see --help)
python validate.py                # verify the analysis is correct
python test_monophyly.py          # engine tests
python test_archive_io.py         # input-handling tests
```

## Verifying it

`validate.py` re-derives the tool's answers by brute force, from a different
representation of the tree, and reports any disagreement. It needs no network
and no NCBI download:

```
python validate.py --trials 1000     # more random trees
python validate.py --seed 7          # a different random sample
python validate.py --tree mytree.nwk # also check your own tree
```

It checks the reported "tips to move" against exhaustive search over every
subset of tips, recovers deliberately planted errors, and confirms Newick
round-trips lose nothing. See **docs/VALIDATION.md** for what each check
proves, what it doesn't, and how to cross-check against MonoPhy or ETE3.

## Built-in help

Every panel carries an **ⓘ** — hover for a one-line explanation, click to open
the help window at the relevant section. **F1** or the **Help** menu opens the
same window directly.

Help text size is adjustable with the **A− / A+ / Reset** buttons, with
**Ctrl/Cmd +**, **−**, **0**, or with **Ctrl+scroll**. The choice applies to
every tab and persists while the app is open.

The **Results** panel has the same control. It scales the Summary and Tip
Matching text *and* the row and heading font of the three tables, so switching
tabs never lands on a control that does nothing. Results and help sizes are
independent.

The window has five tabs:

- **Getting Started** — the four-step workflow and how tip names are matched
- **The Tree View** — colours, zoom, rectangular vs circular, and why labels
  disappear when they do
- **Reading Results** — what "tips to move", the conflict types, alternative
  repairs and the merged rows mean
- **Applications** — what kinds of study this tool is for, what organisms it
  covers, and what it can't tell you
- **Project Guide** — the architecture notes from `docs/ONBOARDING.md`

## Reading the results

Results are reported **one row per taxonomic group**, not one per tree node.
The question answered is "is *Fabaceae* monophyletic in this tree, and if not,
what is sitting inside it?" — so a group broken by a single rogue sequence
appears once, rather than repeating at every nested node beneath it.

- **Conflicting Tips** — start here. Tips ranked by how many groups removing
  them would repair. On the bundled test tree one sequence accounts for 10 of
  the 22 conflicts. A tip here is where to look, not proof of a bad sequence:
  outdated taxonomy, hybridization and incomplete lineage sorting cause
  conflicts too.
- **Broken Groups** — sorted by **Tips to move**: the fewest tips that would
  have to move for the group to be a clade. That is the real difficulty, and
  it is not the intruder count — a group whose two sampled tips merely are not
  sister sweeps a large enclosing clade but is still a one-tip problem.
  **Conflict** says whether foreign tips are *inside* the group (`intruded`),
  the group's own members are *elsewhere* (`scattered`), or both.
  **Support** says whether a supported branch contradicts the group, or only
  weak branches do (see below); weak-only rows sort last.
- **Merge groups with identical tips** — when several nested ranks contain
  exactly the same sampled tips the tree cannot tell them apart, so they share
  one expandable row (`Apiales … Apioideae`). Untick to see every rank.
- **Show rank** / **Conflict** filter the table; ranks include `clade` and
  `no rank`, which is where most deep conflicts live. **Conflict** can also
  show supported conflicts only, or weak-branch-only ones.
- **Tip Matching** — which tips failed to match, which matched only after
  dropping name tokens (so they sit at a coarser rank than the label implies),
  which names are ambiguous across unrelated taxa, and which taxon IDs NCBI
  has merged into newer ones. Tips that match nothing are reported, never
  counted as intruders.

**One tip explains it** names the tip whose removal alone repairs the group. When several tips would
repair a group equally well the cell reads `Cucurbita_pepo or 2 others`, and
selecting the row lists them all — with only two sampled members of a group
the tree genuinely cannot say which one moved, and naming one as the answer
would be a guess dressed as a finding.

### Weak branches, species, and placeholder taxa

An unresolved node fails to confirm a group but does not contradict it.
Internal branches of length 1e-05 or less — the placeholder length RAxML and
IQ-TREE write where the data hold no signal — are therefore treated as
polytomies: a group that fits the tree once they collapse is reported as
"contradicted only by weak branches" rather than broken, and never used to
blame a tip. `--min-support 70` (on the tree's own scale; `0.95` for
posteriors) does the same for poorly supported branches, and
`--min-branch-length -1` turns the length test off. In the GUI both live under
**Analysis > Weak Branch Settings** and are remembered between runs.

Species and infraspecific conflicts get their own section in the report:
non-monophyletic species are common biology, so they weigh less than a broken
genus or family. NCBI placeholder bins ("incertae sedis", "unclassified",
"environmental samples") are never tested — they only collect organisms of
unknown placement.

Bootstrap/support values in the input are preserved. Labelled trees are written
as `Fabaceae/98` by default; `--label-format taxon` or `support` writes just one
or the other.

## The tree view

Branches are green where the tree agrees with the taxonomy and red where it
conflicts; tip colour reports match quality. Zoom runs 10%–400% in discrete
steps (`−` / `+`, or Ctrl/Cmd + scroll), **Fit** sizes the tree to the panel,
**Expand** gives it the whole column, and **F11** is full screen.

**View** switches between a rectangular cladogram and a circular fan. The fan
is the one to use on large trees: the bundled 232-tip tree needs about nine
screenfuls of scrolling as a rectangle and under three as a circle, and that
gap widens with tree size.

Labels never overlap — bigger clades win the spot, shifted names get a leader
line, and anything that still won't fit is dropped with the count shown beside
the display toggles rather than vanishing silently.

## Taxonomy input

The NCBI taxonomy can be supplied in whatever form you have it in. There is no
need to unpack the download first:

```
python main.py --download                             # fetch from NCBI and cache
python main.py --tree t.nwk --taxdump taxdump.tar.gz
python main.py --tree t.nwk --taxdump taxdmp.zip
python main.py --tree t.nwk --names names.dmp.gz --nodes nodes.dmp.gz
python main.py --tree t.nwk --cache                   # reuse the cached database
```

Accepted: `taxdump.tar.gz`, `taxdump.tgz`, `taxdump.tar`, `taxdmp.zip`, and
plain or gzipped `names.dmp` / `nodes.dmp`. The format is detected from the file
contents rather than its name, so a renamed archive still works. In the GUI, use
**Taxonomy DB > Load Archive**, or just pick the archive in **Load from Files**.

## Tree input

Trees may be plain or compressed (`tree.nwk.gz`, `tree.tre.gz`), and batch mode
accepts either a directory or an archive of trees.

Tips may be species names (`Homo_sapiens`) or bare NCBI taxon IDs (`9606`).
Retired IDs are followed to their current taxon via `merged.dmp` and listed
under "Outdated Taxon IDs" in the report; reports show ID tips with their
name, e.g. `312418 (Macrothumia kuhlmannii)`.

```
python main.py --batch trees/ --cache
python main.py --batch trees.tgz --cache
```
