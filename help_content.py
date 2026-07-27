"""
help_content.py - In-app documentation.

Kept as plain data so the GUI stays layout code and the wording can be edited
without touching widget plumbing. Two audiences:

- PANEL_HELP: one-line hover tooltips, for someone who is looking at a control
  right now and wants to know what it does.
- HELP_SECTIONS: the Help window, for someone who wants to know what the tool
  is for and how to read what it produced.
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
        "Clade names sit on chips, support values are the small grey numbers.\n"
        "Selecting a row in the results below highlights that clade here.",
        "Reading Results",
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
   There is no need to unpack a download first.

2. Add one or more trees.
   "Add Tree(s)" takes Newick files, plain or gzipped. For many trees use
   File > Open Multiple Trees (a directory) or File > Open Trees from Archive
   (a .tar.gz / .zip). Click a file in the list to display it.

3. Run the check.
   "Run Monophyly Check" analyses the selected tree. "Run Batch Analysis"
   runs every loaded tree and writes a summary line for each.

4. Read the results, starting with the Conflicting Tips tab.

TIP NAMES

Tips are matched to NCBI by name. "Homo_sapiens", "homo sapiens" and
"Homo sapiens" all match. Labels carrying extra tokens still resolve:
"Homo_sapiens|COX1" and "Homo_sapiens_KX12345" both find Homo sapiens, and
"Homo_notaspecies" falls back to the genus.

Anything matched loosely is listed in the Tip Matching tab, because a
genus-level match is a weaker claim than a species-level one and the analysis
should not quietly pretend otherwise.

EXPORTING

File > Export Labeled Tree writes Newick with taxonomy names added to the
clades, keeping the bootstrap values your input had (Anthocerotophyta/100).
File > Export Report writes the full text report.
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
    Filled in when a single tip is at fault, blank when the conflict needs
    two or more moves.

Intruders / Outliers
    The raw evidence: foreign tips inside the group, and members of the group
    sitting outside it.

MERGED ROWS

When several nested ranks contain exactly the same sampled tips, the tree
cannot tell them apart, and reporting each separately is several ways of
saying one thing. Those share one expandable row, labelled with the coarsest
and finest name (Apiales ... Apioideae). Click the triangle to see the taxa
merged. Untick "Merge groups with identical tips" to see every rank.

TIP MATCHING

Which tips matched nothing, which matched only after dropping name tokens
(so they sit at a coarser rank than their label suggests), and which names are
ambiguous because unrelated taxa share them.

Tips that match nothing are reported but never counted as intruders. A
sequence we could not place is not evidence that a group is broken, and
blaming it would let one unmatched tip appear to break every group it happens
to sit inside.

THE TREE VIEW

Green branches agree with the taxonomy, red conflict, grey are unresolved.
Tip colour marks match quality: grey unmatched, amber matched loosely, red
when the tip intrudes into the group you have selected. Support values are the
small grey numbers; clade names sit on chips.

"Branch lengths" switches between a cladogram (x from depth, tips aligned) and
a phylogram (x from cumulative branch length).
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

WHAT IT CANNOT TELL YOU

It finds disagreement. It cannot tell you which of three things is at fault:
a wrong sequence, a wrong tree, or a taxonomy that has not caught up with the
biology. Deep groups showing as broken often reflect sampling and rooting
rather than a real result, and the tool flags them the same way it flags a
single contaminated tip. Interpretation stays with you.

It also does not yet weigh conflicts by branch support. Support values are
preserved and drawn, but a conflict on a 30%-support branch and one on a
100%-support branch are currently reported alike, and only the second is
worth acting on.

IT DOES NOT INFER TREES

PhyLabeler reads a tree you already have. Build it with RAxML, IQ-TREE,
MrBayes, ASTRAL or anything else that writes Newick.
"""


HELP_SECTIONS = [
    ("Getting Started", _GETTING_STARTED),
    ("Reading Results", _READING_RESULTS),
    ("Applications", _APPLICATIONS),
]
