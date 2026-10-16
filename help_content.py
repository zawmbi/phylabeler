"""
help_content.py - In-app documentation.

Kept as plain data so the GUI stays layout code and the wording can be edited
without touching widget plumbing. Two audiences:

- PANEL_HELP: one-line hover tooltips, for someone who is looking at a control
  right now and wants to know what it does.
- HELP_SECTIONS: the Help window, for someone who wants to know what the tool
  is for and how to read what it produced.

Every key in PANEL_HELP must correspond to a `self._panel(...)` call in
gui.py: `_panel` indexes this dict directly, so a missing key is a crash at
startup rather than a missing tooltip.
"""

# Panel key -> (tooltip shown on hover, Help window section to open on click).
PANEL_HELP = {
    "taxonomy": (
        "The NCBI taxonomy your tree is checked against.\n\n"
        "Load Files: pick names.dmp and nodes.dmp, or an un-extracted\n"
        "taxdump.tar.gz / taxdmp.zip.\n"
        "Download NCBI: fetch it once (~60 MB) and cache it.\n"
        "Load Cache: reuse the copy already downloaded.",
        "Getting Started",
    ),
    "trees": (
        "The Newick trees to analyse.\n\n"
        "Files may be plain or compressed (.nwk.gz, .tre.gz, .bz2, .xz). You\n"
        "can add several at once, a whole directory, or an archive of trees\n"
        "via File > Open Trees from Archive. Tips may be species names or\n"
        "bare NCBI taxon IDs.",
        "Getting Started",
    ),
    "analysis": (
        "Run the monophyly check on the selected tree, or on every loaded\n"
        "tree at once.\n\n"
        "The filters below control which rows the results tables show; they\n"
        "never change the analysis itself. Which branches count as too weak\n"
        "to contradict a group is set under Analysis > Weak Branch Settings.",
        "Reading Results",
    ),
    "search": (
        "Ring every tip AND clade name containing this text, and scroll the\n"
        "tree until at least one of them is on screen - in either view.\n\n"
        "Useful for checking where a suspect sequence actually sits.",
        "Reading Results",
    ),
    "tree_view": (
        "Branches are coloured by monophyly status: green agrees with the\n"
        "taxonomy, red conflicts, grey is unresolved.\n\n"
        "Zoom in discrete steps with - and +, or Ctrl/Cmd + scroll. Fit sizes\n"
        "the whole tree to the panel; Expand gives the tree the whole column.\n"
        "View switches between a rectangular cladogram and a circular fan,\n"
        "which suits large trees far better.",
        "The Tree View",
    ),
    "results": (
        "One row per taxonomic group, not per tree node.\n\n"
        "Broken Groups is sorted so supported, cheapest fixes come first;\n"
        "rows contradicted only by weak branches come last. Conflicting Tips\n"
        "names the tips behind the most conflicts - a bad sample, or outdated\n"
        "taxonomy, hybridization or lineage sorting.",
        "Reading Results",
    ),
}


_GETTING_STARTED = """\
GETTING STARTED

PhyLabeler asks one question: do the clades in your tree agree with the NCBI
taxonomy, and where they disagree, what is responsible?

FOUR STEPS

1. Load the taxonomy.
   First run only. "Download NCBI" fetches the taxonomy dump (~60 MB
   compressed) and caches it; afterwards the app loads that cache
   automatically on startup and you can ignore this panel.

   If you already have the files, "Load Files" accepts any of these forms:
     - names.dmp and nodes.dmp
     - gzipped copies of them (names.dmp.gz)
     - the un-extracted taxdump.tar.gz, taxdump.tgz or taxdmp.zip
   There is no need to unpack a download first. The format is read from the
   file's contents, so a renamed archive still works.

2. Add one or more trees.
   "Add Tree(s)" takes Newick files (.nwk, .tre, .tree, .newick), plain or
   compressed - tree.tre.gz, tree.nwk.bz2 and so on. For many trees use
   File > Open Multiple Trees (a directory) or File > Open Trees from Archive
   (a .tar.gz / .zip). Click a file in the list to display it - adding does
   not select.

3. Run the check.
   "Run Monophyly Check" analyses the selected tree. "Run Batch Analysis"
   runs every loaded tree and writes a summary line for each. The buttons
   grey out while work is in progress.

4. Read the results, starting with the Conflicting Tips tab.

Optional: Analysis > Weak Branch Settings decides which branches are too weak
to count against a group (see "Reading Results"). The default suits most
RAxML and IQ-TREE output; set a support threshold if your tree has bootstrap
or posterior values.

TIP NAMES

Tips are matched to NCBI by name. "Homo_sapiens", "homo sapiens" and
"Homo sapiens" all match. Labels carrying extra tokens still resolve:
"Homo_sapiens|COX1" and "Homo_sapiens_KX12345" both find Homo sapiens, and
"Homo_notaspecies" falls back to the genus.

Tips may also be bare NCBI taxon IDs ("9606"), as written by many
NCBI-backed pipelines. An ID that NCBI has since merged into a newer one is
followed to the current taxon and listed under "Outdated taxon IDs" in the
Tip Matching tab. A numeric label is never matched as a name - a few strain
and isolate names in NCBI are plain numbers, and pairing an ID with one of
those would place the tip in an unrelated organism. Reports print ID tips
with their name, e.g. "312418 (Macrothumia kuhlmannii)".

Anything matched loosely is listed in the Tip Matching tab, because a
genus-level match is a weaker claim than a species-level one and the analysis
should not quietly pretend otherwise. Misspellings are left unresolved on
purpose - there is no fuzzy matching, because silently matching the wrong
organism is worse than reporting a failure.

MAKING THINGS READABLE

Both the help window and the Results panel have their own A- / A+ / Reset
controls, and respond to Ctrl/Cmd with plus, minus and zero. The Results one
also scales the row and heading font of the tables. The two sizes are
independent and persist while the app is open.

The tree has its own zoom, described under "The Tree View".

COLOUR THEMES

View > Colour Theme offers seven palettes. Your choice is remembered between
runs, so whichever one you pick is the one you get next launch:

    Dark (default)              the original
    Light (cream)               warm and low-glare, good for long sessions
    Coffee (warm brown)         the same idea as the default, without the grey
    Black and white             no hue at all - for printing and photocopying
    Terminal green              high-contrast phosphor green
    Colour-blind safe           orange and blue instead of red and green
    Colour-blind safe (dark)    the same, on a dark background

The last three are not decoration. Conflict status normally rides on red
against green, which is exactly the pair that red-green colour blindness
collapses - around 8% of men - so the colour-blind themes swap in orange and
blue, which stay distinct under every common form. The black-and-white theme
has no hue to work with at all, so it marks conflicting branches with a
dashed line instead: that survives greyscale printing, which colour never
does.

EXPORTING

File > Export Labeled Tree writes Newick with taxonomy names added to the
clades, keeping the bootstrap values your input had (Anthocerotophyta/100).
File > Export Report writes the full text report.

Your input files are never modified. The tool only ever writes where you tell
it to.
"""


_TREE_VIEW = """\
THE TREE VIEW

WHAT THE COLOURS MEAN

The key sits under the tree itself and is drawn from whichever theme you are
using, so it always shows the marks actually on screen rather than the ones
the default theme would have used. Untick "Legend" in the Show row to reclaim
the space; the setting is remembered between runs.

Branches are coloured by the status of the clade below them:
    Monophyletic      agrees with the taxonomy
    Not monophyletic  conflicts with it
    Not assessed      no taxon to judge it against

Tip names carry their own colour, which reports match quality rather than
conflict:
    exact             matched a taxon exactly
    Inexact match     matched only after dropping name tokens, so it sits at
                      a coarser rank than the label suggests
    Unmatched tip     matched nothing and is excluded from the analysis
    Intruding tip     intrudes into the group you currently have selected

The hues themselves vary by theme, which is the point: the colour-blind
palettes carry conflict on orange against blue rather than red against green,
and the black-and-white theme has no hue to spend at all, so it marks a
conflicting branch with a dashed line. The legend shows the dash.

The small grey numbers are the support values from your file. Clade names sit
on opaque chips so they stay readable over the branches.

ZOOM AND FRAMING

Zoom moves in discrete steps from 10% to 400% using - and +, so every click
is a visible change. Ctrl/Cmd + scroll wheel does the same; the plain wheel
scrolls. Zooming holds your position in the tree rather than throwing you
somewhere else.

    Fit       size the whole tree to the panel
    Reset     back to 100%
    Expand    give the tree the entire right-hand column
    F11       full screen

Selecting a clade in the results tables highlights it in the tree and scrolls
to it, in both the rectangular and the circular view. In the circular view the
selection is the wedge of angles the clade occupies, carried out past the tips
so the names in the group fall inside it.

FIND

The Find box sits on the Results bar, next to the tables it searches
alongside. It matches clade names as well as tip names, so typing "hypn"
finds the Hypnanae chip drawn on the canvas rather than reporting nothing.

Every match is ringed on the tree, in both the rectangular and the circular
view, and the view scrolls until at least one of them is on screen - along
both axes, because a clade's name sits far to the left of the tips it spans
and a purely vertical scroll could leave every hit off the edge. Matching a
clade also selects the tips underneath it. The band behind a selection is
drawn only when the matches form one unbroken run: a search scattered across
the tree has no single region to shade, and pretending otherwise would claim
everything between the hits was selected too.

PANELS

The dividers between panels are the controls. Drag one to resize, and
double-click it to collapse the panel it points at or restore it to the size
it had. That applies to the controls column and to the split between the tree
and the results, so any panel can be closed or opened without a button for it.

The file list grows with the number of trees loaded rather than reserving the
whole column up front, so a single tree no longer leaves most of the left-hand
side empty.

Selecting a row in the results tables frames that clade: if it is too big for
the current zoom the view zooms out to fit it, but it never zooms in, so a
zoom you chose deliberately is respected. Use Fit to go the other way.

RECTANGULAR OR CIRCULAR

    Rectangular   a cladogram; tips aligned in a column, names easy to read
    Circular      a fan; tips spread around a circle

The circular view exists because a column does not scale. The same 232-tip
tree needs about nine screenfuls of scrolling as a rectangle and under three
as a fan, and that gap widens the bigger the tree gets. Names radiate outward
and flip on the left half so nothing reads upside down.

"Branch lengths" switches the horizontal axis (or the radius) from tree depth
to cumulative evolutionary distance - a cladogram becomes a phylogram. Leave
it off when you care about reading names, on when you care about distances.

WHEN LABELS DISAPPEAR

Labels are placed so they never sit on top of one another. Bigger clades are
placed first, so when two nested groups compete for the same spot the more
inclusive one keeps it. A name that has to shift gets a thin leader line back
to its own branch.

If a name still cannot fit it is dropped rather than drawn over its
neighbour, and the count is reported next to the display toggles, e.g.
"tip names need more zoom - 34 clade names hidden". Zooming in brings them
back. Tip names specifically disappear when adjacent tips are closer together
than the text is tall; clade names hold a legible size at any zoom.
"""


_READING_RESULTS = """\
READING RESULTS

Results are reported one row per taxonomic GROUP, not one per tree node. A
group broken by a single rogue sequence appears once, instead of repeating at
every nested node beneath it.

CONFLICTING TIPS  -  start here

Tips ranked by how many groups removing them would repair. This is usually the
shortest path to a conclusion: on the bundled test tree one sequence accounts
for 10 of the 22 conflicts, so "look at this accession" replaces "investigate
22 findings". Only supported conflicts count here; a group contradicted only
by weak branches never puts a tip on this list.

A tip at the top of this list is where to look, not proof of a bad sequence.
The usual causes are:
    a misidentified or contaminated sample   often sits on a near-zero
                                             terminal branch, identical to
                                             sequences from another family
    an outdated NCBI classification          e.g. a recently split genus still
                                             nested in its parent
    hybridization or allopolyploidy          e.g. polyploid wheats inside
                                             Aegilops
    incomplete lineage sorting               common among close species
    an error in the tree itself

    Groups repaired   how many conflicts this tip alone is responsible for
    Role              whether it sits inside groups it does not belong to
                      (intruder), is a member found in the wrong place
                      (member), or both
    Intrudes into     how many broken groups it merely sits inside

Ranking by how many broken groups a tip merely sits inside does not work: a
group whose members are scattered drags a huge enclosing clade along, so every
tip inside it would score alike. Asking whether removing the tip actually
repairs the group avoids that entirely.

Click a row to highlight that tip in the tree.

BROKEN GROUPS

Tips to move
    The headline number: the fewest tips that would have to move for the
    group to become a clade. This is the real difficulty, and it is NOT the
    intruder count. A group whose two sampled tips simply are not sister
    sweeps a large enclosing clade - 159 foreign tips in one real case - but
    is still a one-tip problem. Rows are sorted cheapest-first, so the top of
    the table is a work queue.

Conflict
    intruded   foreign tips sit inside an otherwise intact group
    scattered  the group's own members are in different parts of the tree
    both       some of each

Support
    supported            at least one well-supported branch contradicts the
                         group - a real conflict
    weak branches only   not a clade as drawn, but it fits the tree once weak
                         branches are collapsed into polytomies, so the data
                         do not actually reject it. These rows sort last.

    An unresolved node fails to confirm a group but does not contradict it.
    By default a branch of length 1e-05 or less is weak: RAxML and IQ-TREE
    write 1e-06 or 1e-08 where the data hold no signal, and on real trees
    14-24% of internal branches are like that. Add a support threshold under
    Analysis > Weak Branch Settings (or --min-support on the command line).
    The settings are remembered between runs; re-run the analysis after
    changing them.

One tip explains it
    The tip whose removal would repair the group. Sometimes more than one tip
    would do equally well, and the cell then reads "Cucurbita_pepo or 2
    others" - selecting the row lists every candidate in the status bar.
    That is not indecision: with only two sampled members of a group the tree
    genuinely cannot say which of them moved, and reporting one as though it
    were the answer would be a guess dressed as a finding.

    Foreign tips are preferred over the group's own, because a sequence
    sitting where it does not belong is the actionable case. The group's own
    tips are only listed when no intruder alone repairs it.

    Blank means the conflict needs two or more moves and no single tip is at
    fault.

Intruders / Outliers
    The raw evidence: foreign tips inside the group, and members of the group
    sitting outside it.

FILTERS AND MERGING

Show rank      narrow to genus, family, order and so on. Includes "clade" and
               "no rank", which is where most of the deep conflicts live.
Conflict       narrow to supported conflicts, weak-branch-only ones, or to
               intruded, scattered or both. Use "intruded" to hunt foreign
               sequences, "scattered" to find sampling problems.
Merge groups with identical tips
               When several nested ranks contain exactly the same sampled
               tips the tree cannot tell them apart, so they share one
               expandable row (Apiales ... Apioideae). Click the triangle to
               see what was merged. Untick to get every rank back.

Filters change what is displayed, never what was computed.

WHAT IS NOT TESTED

NCBI placeholder bins - "Salvia incertae sedis", "unclassified Poaceae",
"environmental samples" - only collect organisms whose placement is unknown.
Nobody claims they are monophyletic, so they are never tested and never used
to label a clade. Groups with a single sampled tip, and groups covering every
tip, say nothing and are skipped too.

SPECIES-LEVEL CONFLICTS

The exported report lists species and infraspecific conflicts in their own
section. Non-monophyletic species are common in real data - incomplete lineage
sorting, introgression, or a species sampled alongside its own subspecies - so
they weigh less than a broken genus or family.

TIP MATCHING

Which tips matched nothing, which matched only after dropping name tokens
(so they sit at a coarser rank than their label suggests), which names are
ambiguous because unrelated taxa share them, and which taxon IDs are outdated
(merged by NCBI into a newer ID, which was used instead). An ID NCBI has
deleted outright cannot be matched and is listed under "No match".

Tips that match nothing are reported but never counted as intruders. A
sequence we could not place is not evidence that a group is broken, and
blaming it would let one unmatched tip appear to break every group it happens
to sit inside.
"""


_APPLICATIONS = """\
WHAT THIS IS FOR

PhyLabeler finds disagreement between a tree and the NCBI taxonomy, and points
at what is responsible. That serves several quite different purposes.

1. SEQUENCE QC BEFORE A REAL ANALYSIS

The strongest use. Run it on your gene trees and read the Conflicting Tips
tab: you get a ranked list of sequences worth a second look before
concatenation, coalescent analysis or divergence dating.

The diagnostic pattern is a biologically impossible placement on a
WELL-SUPPORTED branch. A poorly supported misplacement is phylogenetic noise,
which is why conflicts across weak branches are reported separately.
A strongly supported one is a data problem - a mislabelled accession, a sample
swap, contamination, or a mis-assembled contig. Catching that here costs
seconds; catching it after an MCMC run does not.

2. PHYLOGENOMIC PIPELINE SCREENING AT SCALE

Batch mode over hundreds of gene trees, then ask whether the same taxon tops
the conflicting-tips list in many trees or in one.

  many trees  -> a bad sample or a bad assembly; fix or drop it
  one or two  -> possibly real biology - paralogy, horizontal transfer,
                 introgression, incomplete lineage sorting - and worth
                 looking at rather than discarding

3. SYSTEMATICS AND TAXONOMIC REVISION

The classic question: which named groups does the molecular data contradict?
With adequate sampling, a genus or family that is genuinely non-monophyletic
is a finding and a candidate for revision.

4. REFERENCE DATABASE AND BARCODING CURATION

The same logic applied to a reference set rather than your own samples: do the
sequences actually cluster by the taxonomy they claim?

WHAT ORGANISMS

The whole NCBI taxonomy - roughly 2.9 million taxa. Bacteria, archaea, fungi,
plants, animals, protists and viruses. If it has an NCBI entry, this works on
it.

Two cautions. Prokaryotic taxonomy is unstable and horizontal gene transfer is
everywhere, so bacterial and archaeal trees will show more conflicts and each
deserves more scepticism - a non-monophyletic bacterial genus is quite often
simply true. Viruses work mechanically, but viral taxonomy is structured very
differently and recombination makes "the tree" a shakier concept to begin
with.

Organisms with no NCBI entry - undescribed species, bare "sp." labels, strain
designations that were never registered - cannot be matched at all, and are
listed in the Tip Matching tab rather than guessed at.

CHECKING THE TOOL ITSELF

A validation script ships with it:

    python validate.py

It re-derives the answers by brute force from a different representation of
the tree and reports any disagreement. Among other things it checks the
reported "tips to move" against an exhaustive search over every subset of
tips, and plants a known error in a clean tree to confirm the tool names the
right sequence. No network and no NCBI download required. `--tree yours.nwk`
runs the invariant checks against your own data.

See docs/VALIDATION.md for what each check establishes, what it does not, and
how to cross-check the monophyly calls against MonoPhy in R or ETE3 in Python.

WHAT IT CANNOT TELL YOU

It finds disagreement. It cannot tell you which of three things is at fault:
a wrong sequence, a wrong tree, or a taxonomy that has not caught up with the
biology - and real biology (hybridization, lineage sorting) produces genuine
conflicts that are not errors at all. Deep groups showing as broken often reflect sampling and rooting
rather than a real result, and the tool flags them the same way it flags a
single contaminated tip. Interpretation stays with you.

WEAK BRANCHES AND POLYTOMIES

An unresolved node does not contradict a group - it simply fails to confirm
it. PhyLabeler therefore treats weak branches as polytomies: by default any
internal branch of length 1e-05 or less, the placeholder length RAxML and
IQ-TREE write where the data hold no signal. A group that fits the tree once
those collapse is reported as "weak branches only" rather than broken, sorted
after the supported conflicts, and never used to blame a tip.

Support is not used by default because trees mix scales (0-1 posteriors,
0-100 bootstraps). Set a threshold under Analysis > Weak Branch Settings, or
--min-support 70 on the command line, to treat poorly supported branches the
same way. --min-branch-length -1 switches the length test off.

Species and infraspecific groups are listed separately in the report.
Non-monophyletic species are common biology (incomplete lineage sorting,
introgression, subspecies sampled with their species), so weigh them less
than a broken genus or family.

And rooting matters: change the root and you change which groups look
monophyletic, so root deliberately before running this.

IT DOES NOT INFER TREES

PhyLabeler reads a tree you already have. Build it with RAxML, IQ-TREE,
MrBayes, ASTRAL or anything else that writes Newick.
"""


HELP_SECTIONS = [
    ("Getting Started", _GETTING_STARTED),
    ("The Tree View", _TREE_VIEW),
    ("Reading Results", _READING_RESULTS),
    ("Applications", _APPLICATIONS),
]
