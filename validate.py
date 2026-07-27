#!/usr/bin/env python3
"""
validate.py - Independent verification of PhyLabeler's analysis.

The unit tests in test_monophyly.py check that the code does what its author
expected. That is not the same as checking that the ANSWER IS RIGHT, and it is
not what a sceptical reviewer should accept. This script does the second job:
it re-derives the answers by brute force, from a different data
representation, and reports where the two disagree.

Five independent lines of evidence:

  1. Worked examples      Trees small enough to check by eye, printed with
                          the answer, so you can verify them yourself.
  2. Brute-force oracle   Random trees. The reported "tips to move" is checked
                          against exhaustive search over every subset of tips.
                          This is the strongest check: the fast algorithm and
                          the oracle share no code.
  3. Spike-in recovery    Start from a tree that matches the taxonomy exactly,
                          move one known taxon, and see whether the tool names
                          that taxon. Ground truth is planted, so this
                          measures sensitivity, not self-consistency.
  4. Invariants           Properties that must hold on any input at all.
  5. Round-trip           Parsing and writing Newick loses nothing.

Usage:
    python validate.py                      # all checks, fixed seed
    python validate.py --trials 500         # more random trees
    python validate.py --seed 7             # a different random sample
    python validate.py --tree mytree.nwk    # also check invariants on your
                                            # own tree (needs a taxonomy)

Exit status is 0 only if every check passes.
"""

import argparse
import itertools
import random
import sys

from tree_parser import Node, parse_newick, parse_newick_file
from monophyly import MonophylyChecker


# ---------------------------------------------------------------------------
# A taxonomy stub. check_taxa() takes lineages directly and only consults the
# database for display names, so the random trials need no NCBI download.
# ---------------------------------------------------------------------------

class StubDB:
    def __init__(self, ranks=None):
        self.code_to_rank = ranks or {}

    def get_scientific_name(self, code):
        return code

    def get_lineage_codes(self, code):
        return [code]


# ---------------------------------------------------------------------------
# The oracle: a completely separate representation of the same question.
#
# A tree is reduced to the SET OF ITS CLADES - one frozenset of tip labels per
# node. Nothing here uses tip indices, spans, prefix sums, or the bottom-up
# pass that the real implementation is built on, so agreement between the two
# is meaningful rather than circular.
# ---------------------------------------------------------------------------

def clade_sets(root):
    """Every node's set of tip labels."""
    return [frozenset(node.get_tip_labels()) for node in root.post_order()]


def fixable_by_moving(clades, members, moved):
    """Would relocating exactly the tips in `moved` make `members` a clade?

    Moving is deleting plus re-grafting, and a tip can be grafted anywhere.
    So the group becomes a clade precisely when the members left behind
    already form a clade of the reduced tree: graft the moved members onto
    that clade, and the moved non-members anywhere else.

    A single remaining member counts - one tip is a clade - which is why
    "move one of two non-sister tips next to the other" costs one, not the
    number of tips that happen to lie between them.
    """
    remaining = members - moved
    if not remaining:
        return False          # nothing left to graft onto
    return any((clade - moved) == remaining for clade in clades)


def min_moves(clades, members, all_tips, max_k):
    """Fewest tips to relocate to make `members` a clade, by exhaustive search.

    Returns None when no set of size <= max_k works.
    """
    members = frozenset(members)
    for k in range(max_k + 1):
        for moved in itertools.combinations(all_tips, k):
            if fixable_by_moving(clades, members, frozenset(moved)):
                return k
    return None


# ---------------------------------------------------------------------------
# Random inputs
# ---------------------------------------------------------------------------

def random_tree(rng, labels):
    """A random rooted binary topology over the given tip labels."""
    pool = list(labels)
    rng.shuffle(pool)
    while len(pool) > 1:
        i = rng.randrange(len(pool))
        a = pool.pop(i)
        j = rng.randrange(len(pool))
        b = pool.pop(j)
        pool.append(f"({a},{b})")
    return parse_newick(pool[0] + ";")


def nested_lineages(labels, group_size):
    """Give each tip a lineage: its own code, its genus, its family, root."""
    lineages = {}
    for index, label in enumerate(labels):
        genus = f"G{index // group_size}"
        family = f"F{index // (group_size * 2)}"
        lineages[label] = [label, genus, family, "root"]
    return lineages


def matching_tree(labels, group_size):
    """A tree that agrees with nested_lineages() perfectly - zero conflicts."""
    def join(items):
        while len(items) > 1:
            items = ["(" + items[i] + "," + items[i + 1] + ")"
                     if i + 1 < len(items) else items[i]
                     for i in range(0, len(items), 2)]
        return items[0]

    genera = [join(labels[i:i + group_size])
              for i in range(0, len(labels), group_size)]
    families = [join(genera[i:i + 2]) for i in range(0, len(genera), 2)]
    return join(families) + ";"


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.failures = []

    def check(self, label, ok, detail=""):
        print(f"  {'PASS' if ok else 'FAIL'}  {label}"
              f"{'  -  ' + detail if detail else ''}")
        if not ok:
            self.failures.append(label)
        return ok

    def note(self, text):
        print(f"        {text}")


def heading(text):
    print(f"\n{'=' * 72}\n{text}\n{'=' * 72}")


# ---------------------------------------------------------------------------
# 1. Worked examples
# ---------------------------------------------------------------------------

EXAMPLES = [
    (
        "A clean genus",
        "((a1,a2),(b1,b2));",
        "G0",
        "a1 and a2 are sister and nothing else is in with them.",
        dict(status="monophyletic", cost=0, clusters=1),
    ),
    (
        "One foreign tip nested inside",
        "((a1,b1),a2);",
        "G0",
        "b1 sits between a1 and a2. Delete b1 and the genus is whole:\n"
        "        one tip to move, and b1 is named.",
        dict(status="not_monophyletic", cost=1, clusters=2,
             conflict_type="intruded", culprit="b1"),
    ),
    (
        "One member parked elsewhere",
        "((a1,(b1,b2)),a2);",
        "G0",
        "a1 and a2 are not sister. The enclosing clade holds 4 tips, but\n"
        "        moving ONE tip fixes it - the count of foreign tips inside\n"
        "        is not the difficulty.",
        dict(status="not_monophyletic", cost=1, clusters=2,
             conflict_type="scattered"),
    ),
    (
        "Two foreign tips inside",
        "(((a1,b1),a2),(b2,a3));",
        "G0",
        "Move b1 and b2 out and a1,a2,a3 are a clade. Two moves.",
        dict(status="not_monophyletic", cost=2, conflict_type="intruded",
             culprit=None),
    ),
    (
        "Both at once",
        "((a2,(b2,(a3,a4))),(a1,(b1,b3)));",
        "G0",
        "The best reading of this tree is the clade (a2,b2,a3,a4): b2 is\n"
        "        inside it and doesn't belong, and a1 belongs but is outside.\n"
        "        Neither repair alone is enough - the answer mixes them, and\n"
        "        no single tip is to blame.",
        dict(status="not_monophyletic", cost=2, conflict_type="both",
             culprit=None),
    ),
    (
        "An unplaced tip is not an intruder",
        "(((a1,zzz),a2),(b1,b2));",
        "G0",
        "zzz matched no taxon. We cannot say it does not belong, so it is\n"
        "        reported separately and the genus stays monophyletic.",
        dict(status="monophyletic", cost=0, clusters=1),
    ),
]


def run_worked_examples(report):
    heading("1. WORKED EXAMPLES  (check these by eye)")
    print("Tips a1,a2,a3 are genus G0; b1,b2 are genus G1; zzz matches nothing.\n")

    checker = MonophylyChecker(StubDB())
    for title, newick, group, explanation, expected in EXAMPLES:
        tree = parse_newick(newick)
        labels = tree.get_tip_labels()
        lineages = {}
        unresolved = []
        for label in labels:
            if label.startswith("a"):
                lineages[label] = ["G0", "root"]
            elif label.startswith("b"):
                lineages[label] = ["G1", "root"]
            else:
                lineages[label] = []
                unresolved.append(label)

        results = {r.tax_id: r for r in
                   checker.check_taxa(tree, lineages, unresolved)}
        result = results.get(group)

        print(f"\n  {title}")
        print(f"        tree: {newick}")
        print(f"        {explanation}")
        if result is None:
            report.check(f"{title}: group reported", False, "group missing")
            continue

        print(f"        -> {result.status}, "
              f"{result.conflict_cost} tip(s) to move"
              + (f", {result.conflict_type}" if result.conflict_type else "")
              + (f", blame {result.culprit}" if result.culprit else ""))

        ok = True
        for field, want in expected.items():
            got = {"status": result.status,
                   "cost": result.conflict_cost,
                   "clusters": result.cluster_count,
                   "conflict_type": result.conflict_type,
                   "culprit": result.culprit}[field]
            if got != want:
                ok = False
                report.note(f"expected {field}={want!r}, got {got!r}")
        report.check(f"{title}", ok)


# ---------------------------------------------------------------------------
# 2. Brute-force oracle
# ---------------------------------------------------------------------------

def run_oracle(report, trials, seed, n_tips=9, max_k=3):
    heading("2. BRUTE-FORCE ORACLE  (the strongest check)")
    print(f"{trials} random trees of {n_tips} tips. For every taxonomic group")
    print("in each, the reported 'tips to move' is compared against an")
    print("exhaustive search over every subset of tips up to size", max_k, end=".\n")
    print("The oracle works from the tree's clade sets and shares no code")
    print("with the algorithm under test.\n")

    rng = random.Random(seed)
    checker = MonophylyChecker(StubDB())
    labels = [f"t{i}" for i in range(n_tips)]
    lineages = nested_lineages(labels, group_size=3)

    compared = mismatches = 0
    culprits_checked = culprits_bad = 0
    beyond = 0

    for _ in range(trials):
        tree = random_tree(rng, labels)
        clades = clade_sets(tree)
        all_tips = tuple(sorted(tree.get_tip_labels()))

        for result in checker.check_taxa(tree, lineages):
            truth = min_moves(clades, set(result.members), all_tips, max_k)
            if truth is None:
                # Oracle found no repair within max_k, so the true
                # answer exceeds it; the tool must agree it is expensive.
                beyond += 1
                if result.conflict_cost <= max_k:
                    mismatches += 1
                    report.note(
                        f"tool says {result.conflict_cost}, oracle found no "
                        f"repair <= {max_k} for {result.name} in "
                        f"{tree.to_newick()}")
                continue

            compared += 1
            if result.conflict_cost != truth:
                mismatches += 1
                if mismatches <= 5:
                    report.note(
                        f"tool={result.conflict_cost} oracle={truth} "
                        f"group={result.name} tree={tree.to_newick()}")

            # When one tip is blamed, moving exactly that tip must work.
            if result.culprit is not None:
                culprits_checked += 1
                if not fixable_by_moving(
                        clades, frozenset(result.members),
                        frozenset([result.culprit])):
                    culprits_bad += 1
                    report.note(
                        f"blamed {result.culprit} but moving it does not "
                        f"repair {result.name} in {tree.to_newick()}")

    print()
    report.check("reported cost equals exhaustive search",
                 mismatches == 0,
                 f"{compared:,} group-comparisons, {mismatches} disagreements")
    report.check("named culprits really are sufficient",
                 culprits_bad == 0,
                 f"{culprits_checked:,} single-tip claims verified")
    report.note(f"{beyond:,} groups needed more than {max_k} moves "
                f"(oracle skipped, cost only bounded)")


# ---------------------------------------------------------------------------
# 3. Spike-in recovery
# ---------------------------------------------------------------------------

def transplant(tree, victim_label, host_label):
    """Detach one tip and re-graft it as sister to another. Edits in place.

    Done on the tree structure rather than by editing Newick text: with
    labels like t1 and t10, textual substitution silently hits the wrong tip.
    """
    tips = {node.label: node for node in tree.get_tip_nodes()}
    victim, host = tips[victim_label], tips[host_label]

    parent = victim.parent
    parent.remove_child(victim)
    # A node left with a single child is not a real internal node; splice it
    # out so the transplant does not also change the tree's shape.
    if len(parent.children) == 1 and parent.parent is not None:
        survivor = parent.children[0]
        grandparent = parent.parent
        position = grandparent.children.index(parent)
        parent.remove_child(survivor)
        grandparent.remove_child(parent)
        grandparent.children.insert(position, survivor)
        survivor.parent = grandparent

    holder = host.parent
    position = holder.children.index(host)
    joint = Node()
    holder.children[position] = joint
    joint.parent = holder
    joint.add_child(host)
    joint.add_child(victim)
    return tree


def run_spike_in(report, trials, seed, n_tips=16, group_size=4):
    heading("3. SPIKE-IN RECOVERY  (planted ground truth)")
    print("Start from a tree that matches the taxonomy exactly, so the tool")
    print("reports zero conflicts. Then graft ONE known tip into a different")
    print("genus and ask whether the tool names that tip. This is the standard")
    print("way to measure a detector: the answer is known before you look.\n")
    print(f"Genera hold {group_size} tips, so when one leaves the majority")
    print("stays put and the moved tip is the unambiguous culprit.\n")

    rng = random.Random(seed + 1)
    checker = MonophylyChecker(StubDB())
    labels = [f"t{i}" for i in range(n_tips)]
    lineages = nested_lineages(labels, group_size)
    genus_of = {label: lineages[label][1] for label in labels}

    clean_text = matching_tree(labels, group_size)
    clean_results = checker.check_taxa(parse_newick(clean_text), lineages)
    broken_before = [r for r in clean_results if not r.is_monophyletic]
    report.check("the clean tree reports no conflict at all",
                 not broken_before,
                 f"{len(clean_results)} groups, {len(broken_before)} broken")

    detected = named = 0
    for _ in range(trials):
        victim = rng.choice(labels)
        elsewhere = [t for t in labels if genus_of[t] != genus_of[victim]]
        host = rng.choice(elsewhere)

        tree = transplant(parse_newick(clean_text), victim, host)
        assert sorted(tree.get_tip_labels()) == sorted(labels)

        broken = [r for r in checker.check_taxa(tree, lineages)
                  if not r.is_monophyletic]
        if not broken:
            continue
        detected += 1
        if victim in {r.culprit for r in broken if r.culprit}:
            named += 1

    print()
    report.check("moving one taxon is always detected",
                 detected == trials,
                 f"{detected}/{trials} transplants produced a conflict")
    report.check("the transplanted taxon is the one named",
                 named == detected,
                 f"{named}/{detected} named correctly")


# ---------------------------------------------------------------------------
# 4. Invariants
# ---------------------------------------------------------------------------

def check_invariants(report, results, context):
    bad_cluster = bad_cost = bad_sum = bad_overlap = 0
    for r in results:
        if (r.cluster_count == 1) != r.is_monophyletic:
            bad_cluster += 1
        if (r.conflict_cost == 0) != r.is_monophyletic:
            bad_cost += 1
        if r.conflict_cost != len(r.intruders) + len(r.outliers):
            bad_sum += 1
        if set(r.intruders) & set(r.members):
            bad_overlap += 1
        if set(r.intruders) & set(r.unknown_inside):
            bad_overlap += 1

    report.check(f"{context}: one cluster iff monophyletic", bad_cluster == 0)
    report.check(f"{context}: zero cost iff monophyletic", bad_cost == 0)
    report.check(f"{context}: cost = intruders + outliers", bad_sum == 0)
    report.check(f"{context}: intruders exclude members and unplaced tips",
                 bad_overlap == 0)


def run_invariants(report, trials, seed):
    heading("4. INVARIANTS  (must hold for any input)")
    print("Properties that follow from the definitions. If any fails, the")
    print("output is internally inconsistent regardless of the biology.\n")

    rng = random.Random(seed + 2)
    checker = MonophylyChecker(StubDB())
    labels = [f"t{i}" for i in range(12)]
    lineages = nested_lineages(labels, group_size=3)

    collected = []
    for _ in range(trials):
        tree = random_tree(rng, labels)
        collected.extend(checker.check_taxa(tree, lineages))
    check_invariants(report, collected, f"{trials} random trees")

    print()
    print("  Determinism and order-independence:")
    tree = random_tree(random.Random(seed + 3), labels)
    first = checker.check_taxa(tree, lineages)
    second = checker.check_taxa(tree, lineages)
    same = ([(r.tax_id, r.conflict_cost, r.status) for r in first]
            == [(r.tax_id, r.conflict_cost, r.status) for r in second])
    report.check("same input gives the same answer twice", same)

    shuffled = dict(sorted(lineages.items(), key=lambda kv: kv[0], reverse=True))
    third = checker.check_taxa(tree, shuffled)
    same_order = ({(r.tax_id, r.conflict_cost) for r in first}
                  == {(r.tax_id, r.conflict_cost) for r in third})
    report.check("answer does not depend on tip input order", same_order)


# ---------------------------------------------------------------------------
# 5. Newick round-trip
# ---------------------------------------------------------------------------

def run_round_trip(report, seed, trials=200):
    heading("5. NEWICK ROUND-TRIP  (nothing is lost or invented)")
    print("Parse, write, parse again. Topology, branch lengths and support")
    print("values must survive unchanged.\n")

    rng = random.Random(seed + 4)
    labels = [f"t{i}" for i in range(8)]
    bad_topology = bad_support = bad_lengths = 0

    for _ in range(trials):
        tree = random_tree(rng, labels)
        for node in tree.post_order():
            if not node.is_tip:
                node.support = str(rng.randint(50, 100))
            if node.parent is not None:
                node.branch_length = f"{rng.random():.4f}"

        text = tree.to_newick(branch_lengths=True, label_format="support")
        again = parse_newick(text + ";")

        if sorted(clade_sets(tree)) != sorted(clade_sets(again)):
            bad_topology += 1
        before = [n.support for n in tree.post_order() if not n.is_tip]
        after = [n.support for n in again.post_order() if not n.is_tip]
        if before != after:
            bad_support += 1
        lb = [n.branch_length for n in tree.post_order()]
        la = [n.branch_length for n in again.post_order()]
        if lb != la:
            bad_lengths += 1

    report.check("topology survives a round trip", bad_topology == 0,
                 f"{trials} trees")
    report.check("support values survive", bad_support == 0)
    report.check("branch lengths survive", bad_lengths == 0)


# ---------------------------------------------------------------------------
# Optional: the reviewer's own tree
# ---------------------------------------------------------------------------

def run_user_tree(report, path):
    heading(f"6. YOUR OWN TREE  ({path})")
    try:
        from taxonomy_db import TaxonomyDB
    except ImportError:
        report.check("taxonomy available", False, "could not import TaxonomyDB")
        return

    db = TaxonomyDB()
    if not db.load():
        print("  No cached taxonomy. Run:  python main.py --download")
        report.check("taxonomy loaded", False, "no cache")
        return

    tree = parse_newick_file(path)
    checker = MonophylyChecker(db)
    analysis = checker.analyze(tree)
    results = analysis["taxon_results"]
    print(f"  {tree.count_tips()} tips, {len(results)} groups with 2+ tips, "
          f"{len(analysis['unresolved'])} unresolved\n")

    check_invariants(report, results, "your tree")

    print()
    print("  Spot-checking the cheapest fixes against the oracle:")
    clades = clade_sets(tree)
    all_tips = tuple(sorted(tree.get_tip_labels()))
    cheap = sorted((r for r in results if not r.is_monophyletic),
                   key=lambda r: r.conflict_cost)[:8]
    mismatches = 0
    for r in cheap:
        truth = min_moves(clades, set(r.members), all_tips, max_k=1)
        expected = r.conflict_cost if r.conflict_cost <= 1 else None
        got = truth
        agree = (got == expected)
        mismatches += not agree
        print(f"    {r.name[:34]:36} tool={r.conflict_cost:<3} "
              f"oracle={'>1' if truth is None else truth:<3} "
              f"{'ok' if agree else 'MISMATCH'}")
    report.check("your tree's one-tip fixes match the oracle", mismatches == 0)


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Independently verify PhyLabeler's analysis.")
    parser.add_argument("--trials", type=int, default=200,
                        help="random trees per section (default 200)")
    parser.add_argument("--seed", type=int, default=1,
                        help="random seed (default 1; change it and re-run)")
    parser.add_argument("--tree", type=str,
                        help="also check invariants on your own Newick tree")
    args = parser.parse_args()

    print(__doc__.split("Usage:")[0].strip())

    report = Report()
    run_worked_examples(report)
    run_oracle(report, args.trials, args.seed)
    run_spike_in(report, args.trials, args.seed)
    run_invariants(report, args.trials, args.seed)
    run_round_trip(report, args.seed)
    if args.tree:
        run_user_tree(report, args.tree)

    heading("VERDICT")
    if report.failures:
        print(f"{len(report.failures)} CHECK(S) FAILED:")
        for name in report.failures:
            print(f"  - {name}")
        print("\nDo not trust the output until these are explained.")
        return 1

    print("Every check passed.")
    print("\nWhat this does and does not establish:")
    print("  It does establish that 'tips to move' is the true minimum, that")
    print("  a named culprit really repairs its group, that a planted error")
    print("  is recovered, and that the output is internally consistent.")
    print("  It does NOT establish that a conflict is biologically real -")
    print("  only that the tree and the taxonomy disagree as reported.")
    print("\nRe-run with --seed N to sample different random trees.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
