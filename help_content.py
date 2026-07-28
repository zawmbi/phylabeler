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
        "Files may be plain or gzipped. You can add several at once, a whole\n"
        "directory, or an archive of trees via File > Open Trees from Archive.",
        "Getting Started",
    ),
    "analysis": (
        "Run the monophyly check on the selected tree, or on every loaded\n"
        "tree at once.\n\n"
        "The filters below control which rows the results tables show; they\n"
        "never change the analysis itself.",
        "Reading Results",
    ),
    "search": (
        "Highlight every tip whose name contains this text and scroll the\n"
        "tree to the first match.\n\n"
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
        "Broken Groups is sorted so the cheapest fixes come first. Start with\n"
        "the Conflicting Tips tab: it names the sequences worth rechecking.",
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
   "Add Tree(s)" takes Newick files, plain or gzipped. For many trees use
   File > Open Multiple Trees (a directory) or File > Open Trees from Archive
   (a .tar.gz / .zip). Click a file in the list to display it - adding does
   not select.

3. Run the check.
   "Run Monophyly Check" analyses the selected tree. "Run Batch Analysis"
   runs every loaded tree and writes a summary line for each. The buttons
   grey out while work is in progress.

4. Read the results, starting with the Conflicting Tips tab.

TIP NAMES

Tips are matched to NCBI by name. "Homo_sapiens", "homo sapiens" and
"Homo sapiens" all match. Labels carrying extra tokens still resolve:
"Homo_sapiens|COX1" and "Homo_sapiens_KX12345" both find Homo sapiens, and
"Homo_notaspecies" falls back to the genus.

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

View > Colour Theme offers six palettes, and your choice is remembered
between runs:

    Dark (default)      the original
    Light (cream)       warm and low-glare, good for long sessions
    Dark (warm brown)   the same idea as the default, without the grey
    Black and white     no hue at all - for printing and photocopying
    Terminal green      high-contrast phosphor green
    Colour-blind safe   orange and blue instead of red and green

The last two are not decoration. Conflict status normally rides on red
against green, which is exactly the pair that red-green colour blindness
collapses - around 8% of men - so the colour-blind theme swaps in orange and
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

Branches are coloured by the status of the clade below them:
    green   agrees with the taxonomy
    red     conflicts with it
    grey    could not be assessed

Tip names carry their own colour, which reports match quality rather than
conflict:
    dark    matched a taxon exactly
    amber   matched only after dropping name tokens, so it sits at a coarser
            rank than the label suggests
    grey    matched nothing and is excluded from the analysis
    red     intrudes into the group you currently have selected

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
for 10 of the 22 conflicts, so "recheck this accession" replaces "investigate
22 findings".

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
Conflict       narrow to intruded, scattered or both. Use "intruded" to hunt
               foreign sequences, "scattered" to find sampling problems.
Merge groups with identical tips
               When several nested ranks contain exactly the same sampled
               tips the tree cannot tell them apart, so they share one
               expandable row (Apiales ... Apioideae). Click the triangle to
               see what was merged. Untick to get every rank back.

Filters change what is displayed, never what was computed.

TIP MATCHING

Which tips matched nothing, which matched only after dropping name tokens
(so they sit at a coarser rank than their label suggests), and which names are
ambiguous because unrelated taxa share them.

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
tab: you get a ranked list of sequences to re-check before concatenation,
coalescent analysis or divergence dating.

The diagnostic pattern is a biologically impossible placement on a
WELL-SUPPORTED branch. A poorly supported misplacement is phylogenetic noise.
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
biology. Deep groups showing as broken often reflect sampling and rooting
rather than a real result, and the tool flags them the same way it flags a
single contaminated tip. Interpretation stays with you.

It also does not weigh conflicts by branch support. Support values are
preserved, drawn and exported, but a conflict on a 30%-support branch and one
on a 100%-support branch are reported alike, and only the second is worth
acting on. Check the support on the branches involved before drawing a
conclusion.

A polytomy is reported as a conflict even though an unresolved node does not
actually contradict a group - it simply fails to confirm it. And rooting
matters: change the root and you change which groups look monophyletic, so
root deliberately before running this.

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
