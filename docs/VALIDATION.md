# Validating PhyLabeler

This document is for someone who has no reason to trust this tool and wants to
find out whether it is wrong.

The short version: run `python validate.py`. It re-derives the tool's answers
by brute force, from a different representation of the tree, and tells you
where the two disagree. It needs no network, no NCBI download, and no
dependencies beyond the Python standard library.

```
python validate.py                    # ~10 s
python validate.py --trials 1000      # more random trees
python validate.py --seed 7           # a different random sample
python validate.py --tree mytree.nwk  # also check your own tree
```

Exit status is 0 only if every check passes.

---

## The questions worth asking

### "How do I know the numbers aren't just made up?"

The headline number is **tips to move**: the fewest tips that would have to be
relocated for a taxonomic group to form a clade. That claim is checkable by
exhaustive search, and section 2 of `validate.py` does exactly that.

For every group in every random tree it enumerates **every subset of tips up
to size 3**, asks "if I moved these, would the group become a clade?", and
compares the smallest working answer against what the tool reported. The
oracle is written against the tree's *clade sets* — a `frozenset` of tip
labels per node — and shares no code with the algorithm under test, which
works from tip indices and a bottom-up counting pass.

Typical run: **750 group-comparisons per seed, 0 disagreements**, and it holds
across seeds.

If the two ever disagree, the script prints the offending tree in Newick so
you can reproduce it by hand.

### "Fine, but does it find a real error?"

Section 3 plants one. It builds a tree that matches the taxonomy exactly
(the tool reports zero conflicts — that is checked), then grafts one known
tip into a different genus, and asks whether the tool names *that* tip.

Ground truth is known before you look, so this measures sensitivity rather
than self-consistency. Genera hold 4 tips, so when one leaves, the majority
stays put and the moved tip is unambiguously the culprit.

Typical run: **120/120 transplants detected, 120/120 named correctly.**

This is the check to put on screen if you are demonstrating the tool. It is
the phylogenetics-native validation idiom: simulate a known truth, measure
recovery.

### "What if it's self-consistently wrong?"

That is the right worry, and it is why sections 1 and 2 exist.

Section 1 prints six trees small enough to verify in your head, with the
tool's answer beside each. For example:

```
tree: ((a1,b1),a2);        a1,a2 are genus G0; b1 is genus G1
-> not_monophyletic, 1 tip(s) to move, intruded, blame b1
```

You do not need to trust anything to check that one. Move `b1` and `(a1,a2)`
is a clade.

Section 2 is the real defence: an oracle that cannot be self-consistently
wrong with the implementation, because it does not share any of it.

### "Aren't the unit tests just testing that the code does what it does?"

Largely, yes — that is what unit tests are for, and it is a weaker claim than
people usually take it to be. `test_monophyly.py` (88 assertions) pins
behaviour so it does not regress. `validate.py` is the one that checks the
answer is *right*.

Worth knowing: **the oracle caught a real bug** when it was first run. Several
clades routinely tie on cost — in `((a1,b1),a2)` both "move `b1` out" and
"move `a2` across" cost one — and the tie-break preferred the reading that
split the genus and blamed one of its own members instead of naming the
intruder. Both answers were numerically correct; one was useless. The fix is
in `monophyly.py::_taxon_pass`, and the case is now example 2 in section 1.

### "Why should a group with 159 foreign tips inside it be 'one tip to move'?"

Because the count of foreign tips measures how big the enclosing clade
happens to be, not how wrong the tree is.

If a group has two sampled tips and they are not sister, the smallest clade
containing both may hold most of the tree — but moving *one* of those two
tips next to the other repairs it. That is one move, and the oracle confirms
it independently.

This matters in practice. On the bundled 232-tip test tree, `rosids` scores
**2** under this definition and **59** under the naive
`min(intruders, stray members)`. The naive number would have buried a
two-tip problem at the bottom of the list.

### "What exactly does 'monophyletic' mean here?"

A group is monophyletic when the tips assigned to it form a clade — a single
node whose descendants are exactly those tips. The tool computes this as
"exactly one maximal all-member clade", and `validate.py` checks the
equivalence `clusters == 1 ⟺ monophyletic ⟺ cost == 0` on every group of
every random tree.

Four notes on scope:

- This is **exclusivity relative to the tips you sampled**, not a claim about
  the group in nature. A genus can look monophyletic simply because you
  sampled two of its species.
- Tips that match no taxon are reported but **never counted as intruders**.
  We cannot say an unplaced sequence does not belong, and blaming it would
  let one unmatched tip appear to break every group it sits inside.
- **Weak branches do not count as contradiction.** An unresolved node fails
  to confirm a group but does not reject it. After the strict test above, each
  non-clade group is re-tested with weak branches collapsed into polytomies
  (default: internal length ≤ 1e-05; optionally support below
  `--min-support`). If it is then compatible with the tree — every lineage
  leaving its polytomy holds only members or no members — it is reported as
  *contradicted only by weak branches*, not broken, and never used to blame a
  tip. Its "tips to move" is still the as-drawn cost, so the oracle checks in
  section 2 are unaffected. `test_monophyly.py` covers the length and support
  thresholds and the no-blame rule.
- **NCBI placeholder bins are not groups.** "X incertae sedis",
  "unclassified X" and "environmental samples" are removed from tip lineages
  before any test, so they are never reported as broken.

### "Does it depend on how I ordered my tips, or which run it is?"

Section 4 checks both: the same input twice gives the same answer, and
reversing the input order of the tips does not change any group's cost.

### "Will it mangle my tree file?"

Section 5 parses, writes and re-parses 200 random trees carrying random
support values and branch lengths, and checks that topology, support and
lengths all survive unchanged.

Separately: labelling a clade never overwrites its bootstrap value. The
labelled export writes `Fabaceae/98`, and `--label-format support` writes the
input's values untouched.

### "Can I check it against an established tool?"

Yes, and you should. PhyLabeler deliberately implements the same question as
**MonoPhy** (Schwery & O'Meara 2016), so the monophyletic/non-monophyletic
call is directly comparable.

In R:

```r
install.packages("MonoPhy")
library(MonoPhy)
tree <- ape::read.tree("Rooted.tre")
# taxa: a data frame of tip label -> genus/family
result <- AssessMonophyly(tree, taxonomy = taxa)
GetSummaryMonophyly(result)
```

Or with ETE3 in Python:

```python
from ete3 import Tree
t = Tree("Rooted.tre")
ok, kind, intruders = t.check_monophyly(
    values={"Bupleurum_falcatum", ...}, target_attr="name")
```

Expect the **status** to agree. Do not expect the numbers to: neither tool
reports a "tips to move" cost, and MonoPhy's intruder list is taken at the
MRCA rather than at the cheapest clade, so it will match PhyLabeler's
`mrca_intruder_count`, not its `intruders`.

### "Can I just read the code?"

It is about 3,000 lines of pure standard-library Python with no dependencies.
The part doing the actual work is `MonophylyChecker._taxon_pass` in
`monophyly.py` — one bottom-up pass, roughly 60 lines. `docs/ONBOARDING.md`
explains the algorithm and why the cost is minimised over every node rather
than the two obvious candidates.

---

## What validation does *not* establish

Worth saying plainly, because overclaiming is how tools lose trust.

`validate.py` establishes that the tool computes what it says it computes. It
says nothing about whether a reported conflict is biologically meaningful.

A conflict means the tree and NCBI disagree. It does not tell you which is at
fault:

- the sequence may be mislabelled or contaminated
- the tree may be wrong — poor alignment, poor model, too little signal
- the taxonomy may simply not reflect current understanding
- or the conflict may be real biology: paralogy, horizontal transfer,
  introgression, incomplete lineage sorting

Weak branches are handled (see above), but only by a threshold you choose:
support is off by default because trees mix 0–1 and 0–100 scales, so set
`--min-support` (or Analysis > Weak Branch Settings) if your tree carries
bootstrap or posterior values. A conflict just above the threshold is still
reported as supported — check the branches involved before drawing a
conclusion.

Species-level conflicts are listed separately in the report because they are
so often real biology rather than error.

The diagnostic pattern for a data problem is a **biologically impossible
placement on a well-supported branch**. Everything else needs judgement.

---

## Running the whole suite

```
python validate.py         # correctness of the answers (this document)
python test_monophyly.py   # 132 behavioural assertions
python test_archive_io.py  # 48 input-handling assertions
```
