#!/usr/bin/env python3
"""
test_monophyly.py - Tests for tip resolution, parsing, and monophyly.

Uses a small synthetic taxonomy where the right answer is known by hand, so
the assertions do not depend on the NCBI download being present:

    python test_monophyly.py
"""

import sys

from taxonomy_db import TaxonomyDB, normalize_name
from tree_parser import parse_newick, index_tips, quote_if_needed
from monophyly import MonophylyChecker, is_placeholder_taxon

# root(1)
#  +- Metazoa(100, kingdom)
#  |   +- Chordata(200, phylum) -> Homo(300, genus) -> sapiens(301), erectus(302)
#  |   +- Arthropoda(210, phylum) -> Drosophila(310, genus) -> melanogaster(311),
#  |                                                          simulans(312)
#  +- Plantae(110, kingdom) -> Rosales(220, order) -> Rosa(320, genus) -> canina(321)
NODES = [
    ("1", "1", "no rank"),
    ("100", "1", "kingdom"), ("110", "1", "kingdom"),
    ("200", "100", "phylum"), ("210", "100", "phylum"),
    ("220", "110", "order"),
    ("300", "200", "genus"), ("310", "210", "genus"), ("320", "220", "genus"),
    ("301", "300", "species"), ("302", "300", "species"),
    ("303", "300", "species"), ("304", "300", "species"),
    ("311", "310", "species"), ("312", "310", "species"),
    ("313", "310", "species"),
    ("321", "320", "species"),
]
NAMES = [
    ("1", "root", "scientific name"),
    ("100", "Metazoa", "scientific name"),
    ("110", "Plantae", "scientific name"),
    ("200", "Chordata", "scientific name"),
    ("210", "Arthropoda", "scientific name"),
    ("220", "Rosales", "scientific name"),
    ("300", "Homo", "scientific name"),
    ("310", "Drosophila", "scientific name"),
    ("320", "Rosa", "scientific name"),
    ("301", "Homo sapiens", "scientific name"),
    ("302", "Homo erectus", "scientific name"),
    ("303", "Homo habilis", "scientific name"),
    ("304", "Homo naledi", "scientific name"),
    ("311", "Drosophila melanogaster", "scientific name"),
    ("312", "Drosophila simulans", "scientific name"),
    ("313", "Drosophila yakuba", "scientific name"),
    ("321", "Rosa canina", "scientific name"),
    # A synonym that must never outrank a real scientific name.
    ("311", "Fruit fly", "genbank common name"),
    ("312", "Homo sapiens", "synonym"),
    # A genuine homonym: two different taxa legitimately named "Rosa".
    ("999", "Rosa", "scientific name"),
    # A numeric strain name, the kind a retired taxon ID could collide with.
    ("313", "4040", "equivalent name"),
]
NODES.append(("999", "210", "genus"))

_passed = _failed = 0


def check(label, got, want):
    global _passed, _failed
    if got == want:
        _passed += 1
        print(f"  ok    {label}")
    else:
        _failed += 1
        print(f"  FAIL  {label}\n          got  {got!r}\n          want {want!r}")


def build_db():
    db = TaxonomyDB.__new__(TaxonomyDB)   # skip cache-dir creation
    db.cache_dir = None
    db.name_to_code = {}
    db.code_to_names = {}
    db.code_to_scientific = {}
    db.code_to_parent = {}
    db.code_to_rank = {}
    db.homonyms = {}
    db.merged = {"4001": "301"}          # retired ID -> Homo sapiens
    db.loaded = False
    db._cache_version = 4

    names = "".join(f"{t}\t|\t{n}\t|\t\t|\t{c}\t|\n" for t, n, c in NAMES)
    nodes = "".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t, p, r in NODES)
    db._parse_names_stream(names.splitlines(keepends=True))
    db._parse_nodes_stream(nodes.splitlines(keepends=True))
    return db


def test_name_resolution(db):
    print("\nScientific names win over synonyms and case/underscore variants")
    check("Homo_sapiens", db.lookup_name("Homo_sapiens"), "301")
    check("homo sapiens (case-folded)", db.lookup_name("homo sapiens"), "301")
    check("HOMO_SAPIENS", db.lookup_name("HOMO_SAPIENS"), "301")
    check("Homo  sapiens (extra space)", db.lookup_name("Homo  sapiens"), "301")
    # 302 lists "Homo sapiens" as a synonym; it must not steal the name.
    check("synonym never overrides a scientific name",
          db.lookup_name("Homo sapiens"), "301")
    check("common name still resolves", db.lookup_name("Fruit fly"), "311")
    check("unknown name", db.lookup_name("Nonexistent thing"), None)


def test_taxid_tips(db):
    print("\nTips labelled with NCBI taxon IDs")
    match = db.resolve_tip("311")
    check("bare taxon ID", (match.tax_id, match.method), ("311", "taxid"))
    check("taxon ID counts as exact", match.is_exact, True)
    check("taxon ID takes the scientific name", match.matched_name,
          "Drosophila melanogaster")
    match = db.resolve_tip("4001")
    check("retired ID follows merged.dmp", (match.tax_id, match.method),
          ("301", "merged"))
    check("numeric label never matched as a name", db.resolve_tip("4040"),
          None)
    check("unknown ID", db.resolve_tip("123456"), None)
    check("display_label names an ID", db.display_label("311"),
          "311 (Drosophila melanogaster)")
    check("display_label leaves names alone",
          db.display_label("Homo_sapiens"), "Homo_sapiens")

    print("\nHomonyms are recorded rather than silently resolved")
    check("Rosa is flagged ambiguous", "rosa" in db.homonyms, True)
    check("both Rosa taxa listed",
          sorted(db.get_homonyms("Rosa")), ["320", "999"])
    check("first NCBI entry wins", db.lookup_name("Rosa"), "320")

    print("\nInexact tip matching reports how it matched")
    for label, tax_id, method in [
        ("Homo_sapiens", "301", "exact"),
        ("Homo_sapiens|COX1", "301", "trimmed"),
        ("Homo_sapiens_KX12345", "301", "partial"),
        ("Homo_notaspecies", "300", "genus"),
    ]:
        match = db.resolve_tip(label)
        check(f"{label} -> {tax_id} via {method}",
              (match.tax_id, match.method), (tax_id, method))
    check("garbage stays unresolved", db.resolve_tip("Xyzzy_qqq"), None)

    print("\nNormalisation")
    check("normalize underscores/case",
          normalize_name("Homo_Sapiens"), "homo sapiens")
    check("normalize collapses whitespace",
          normalize_name("  Homo   sapiens "), "homo sapiens")


def test_parsing():
    print("\nParser keeps support separate from labels")
    tree = parse_newick("((A:0.1,B:0.2)95:0.3,C:0.4);")
    internal = tree.children[0]
    check("numeric internal label becomes support", internal.support, "95")
    check("label left free for taxonomy", internal.label, "")

    internal.label = "Hominidae"
    check("export keeps both",
          tree.to_newick(True, "taxon_support"),
          "((A:0.1,B:0.2)Hominidae/95:0.3,C:0.4)")
    check("export taxon only",
          tree.to_newick(True, "taxon"),
          "((A:0.1,B:0.2)Hominidae:0.3,C:0.4)")
    check("export support only leaves the input untouched",
          tree.to_newick(True, "support"),
          "((A:0.1,B:0.2)95:0.3,C:0.4)")

    print("\nParser robustness")
    check("comments ignored",
          parse_newick("((A,B)[&rate=0.4]90,C);").count_tips(), 3)
    check("polytomy", parse_newick("(A,B,C,D);").count_tips(), 4)
    check("scientific notation branch length",
          parse_newick("(A:1e-3,B:2.5E-4);").children[0].branch_length, "1e-3")
    check("quoted label with spaces",
          parse_newick("('Homo sapiens',B);").children[0].label, "Homo sapiens")
    check("labels needing quotes are re-quoted on output",
          quote_if_needed("Homo sapiens"), "'Homo sapiens'")

    print("\nDeep trees do not overflow the stack")
    sys.setrecursionlimit(1000)
    ladder = "t0"
    for i in range(1, 5000):
        ladder = f"({ladder},t{i})"
    deep = parse_newick(ladder + ";")
    check("5000-tip ladder parses", deep.count_tips(), 5000)
    check("5000-tip ladder renders", len(deep.to_newick()) > 0, True)


def test_index_tips():
    print("\nTip indexing")
    tree = parse_newick("((A,B),(C,D));")
    tips, spans = index_tips(tree)
    check("tip order", [t.label for t in tips], ["A", "B", "C", "D"])
    check("root span", spans[id(tree)], (0, 3, 4))
    check("left clade span", spans[id(tree.children[0])], (0, 1, 2))
    check("right clade span", spans[id(tree.children[1])], (2, 3, 2))


def test_placeholder_taxa():
    print("\nNCBI placeholder bins are never tested as groups")
    check("incertae sedis", is_placeholder_taxon("Salvia incertae sedis"),
          True)
    check("unclassified", is_placeholder_taxon("unclassified Poaceae"), True)
    check("real genus", is_placeholder_taxon("Salvia"), False)

    # Homo naledi and Homo habilis sit in a "Homo incertae sedis" bin, which
    # the tree splits; the bin must not be reported as a broken group.
    db = build_db()
    db.code_to_parent.update({"305": "300", "303": "305", "304": "305"})
    db.code_to_rank["305"] = "no rank"
    db.code_to_scientific["305"] = "Homo incertae sedis"
    tree = parse_newick("(((Homo_habilis,Homo_sapiens),Homo_naledi),"
                        "(Drosophila_melanogaster,Drosophila_simulans));")
    out = MonophylyChecker(db).analyze(tree)
    names = {r.name for r in out["taxon_results"]}
    check("bin not among tested groups", "Homo incertae sedis" in names,
          False)
    check("real parent still tested", "Homo" in names, True)
    labels = {r.mrca_name for r in out["node_results"]}
    check("no node labelled with the bin", "Homo incertae sedis" in labels,
          False)


def test_weak_branches(db):
    print("\nConflicts across weak branches are not evidence")
    # Homo erectus sits with Drosophila melanogaster, but only across a
    # zero-length branch; Drosophila is split across real branches.
    newick = ("((Homo_sapiens,(Homo_erectus,Drosophila_melanogaster)"
              "{weak}):0.1,(Drosophila_simulans,Rosa_canina):0.1);")
    tree = parse_newick(newick.format(weak=":0.000001"))
    by_name = {r.name: r for r in
               MonophylyChecker(db).analyze(tree)["taxon_results"]}
    check("zero-length conflict is unsupported",
          by_name["Homo"].status, "unsupported_conflict")
    check("still costs a move as drawn", by_name["Homo"].conflict_cost, 1)
    check("real conflict stays broken",
          by_name["Drosophila"].status, "not_monophyletic")

    tree = parse_newick(newick.format(weak=":0.000001"))
    by_name = {r.name: r for r in MonophylyChecker(
        db, min_branch_length=None).analyze(tree)["taxon_results"]}
    check("threshold off: every branch counts",
          by_name["Homo"].status, "not_monophyletic")

    tree = parse_newick(newick.format(weak="50:0.2"))
    by_name = {r.name: r for r in MonophylyChecker(
        db, min_support=70).analyze(tree)["taxon_results"]}
    check("low support counts as weak",
          by_name["Homo"].status, "unsupported_conflict")
    tree = parse_newick(newick.format(weak="95:0.2"))
    by_name = {r.name: r for r in MonophylyChecker(
        db, min_support=70).analyze(tree)["taxon_results"]}
    check("high support counts", by_name["Homo"].status, "not_monophyletic")

    out = MonophylyChecker(db).analyze(
        parse_newick(newick.format(weak=":0.000001")))
    blamed = {t.label for t in out["conflicting_tips"]}
    check("weak-only conflicts blame no tip",
          blamed.isdisjoint({"Homo_sapiens", "Homo_erectus"}), True)


def test_monophyly(db):
    checker = MonophylyChecker(db)

    print("\nA tree where every group is monophyletic")
    tree = parse_newick(
        "(((Homo_sapiens,Homo_erectus),"
        "(Drosophila_melanogaster,Drosophila_simulans)),Rosa_canina);")
    out = checker.analyze(tree)
    by_name = {r.name: r for r in out["taxon_results"]}
    check("no unresolved tips", out["unresolved"], [])
    check("Homo monophyletic", by_name["Homo"].status, "monophyletic")
    check("Drosophila monophyletic",
          by_name["Drosophila"].status, "monophyletic")
    check("Metazoa monophyletic", by_name["Metazoa"].status, "monophyletic")
    check("no broken groups",
          [r.name for r in out["taxon_results"] if not r.is_monophyletic], [])

    print("\nA tree where a fly is nested inside Homo")
    tree = parse_newick(
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));")
    out = checker.analyze(tree)
    by_name = {r.name: r for r in out["taxon_results"]}
    check("Homo is broken", by_name["Homo"].status, "not_monophyletic")
    check("the fly is named as the intruder",
          by_name["Homo"].intruders, ["Drosophila_melanogaster"])
    check("Homo still has 2 members", len(by_name["Homo"].members), 2)
    check("Drosophila is broken too",
          by_name["Drosophila"].status, "not_monophyletic")

    print("\nOne row per group, not one per node")
    names = [r.name for r in out["taxon_results"]]
    check("no duplicate group rows", len(names), len(set(names)))
    check("groups with a single sampled tip are omitted",
          "Rosa" in names, False)

    print("\nLabelling preserves support")
    tree = parse_newick("((Homo_sapiens:0.1,Homo_erectus:0.2)99:0.3,"
                        "Drosophila_melanogaster:0.4);")
    out = checker.analyze(tree)
    checker.label_tree(tree, out["node_results"])
    node = tree.children[0]
    check("clade labelled", node.label, "Homo")
    check("support survives labelling", node.support, "99")
    check("both written out",
          tree.to_newick(True, "taxon_support"),
          "((Homo_sapiens:0.1,Homo_erectus:0.2)Homo/99:0.3,"
          "Drosophila_melanogaster:0.4)")

    print("\nTips that resolve to no lineage are reported, not dropped")
    tree = parse_newick("((Homo_sapiens,Homo_erectus),Xyzzy_qqq);")
    out = checker.analyze(tree)
    check("unknown tip listed as unresolved",
          out["unresolved"], ["Xyzzy_qqq"])


def test_conflict_costs(db):
    checker = MonophylyChecker(db)

    def groups_of(newick):
        out = checker.analyze(parse_newick(newick))
        return out, {r.name: r for r in out["taxon_results"]}

    print("\nClusters are clades, not runs of adjacent tips")
    # The two Homo tips sit at positions 1 and 2, adjacent in tip order, but
    # no node holds exactly them. An index-run heuristic says 1 cluster.
    _out, by_name = groups_of(
        "(((Drosophila_melanogaster,Homo_sapiens),Homo_erectus),Rosa_canina);")
    homo = by_name["Homo"]
    check("adjacent but not a clade -> 2 clusters", homo.cluster_count, 2)
    check("cluster sizes are member counts",
          sorted(len([t for t in c.get_tip_labels()
                      if t in set(homo.members)]) for c in homo.clusters),
          [1, 1])

    print("\nA nested intruder: cheap to fix, badly fragmented")
    _out, by_name = groups_of(
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));")
    homo = by_name["Homo"]
    check("intruder named", homo.intruders, ["Drosophila_melanogaster"])
    check("no members stranded", homo.outliers, [])
    check("one tip to move", homo.conflict_cost, 1)
    check("type", homo.conflict_type, "intruded")
    check("culprit", homo.culprit, "Drosophila_melanogaster")

    print("\nTwo members that simply are not sister")
    # The Proteales shape: a huge enclosing clade, but one tip out of place.
    _out, by_name = groups_of(
        "((Homo_sapiens,(Drosophila_melanogaster,Drosophila_simulans)),"
        "Homo_erectus);")
    homo = by_name["Homo"]
    check("no intruders blamed", homo.intruders, [])
    check("one member sits outside", len(homo.outliers), 1)
    check("still only one tip to move", homo.conflict_cost, 1)
    check("type", homo.conflict_type, "scattered")
    check("the enclosing clade is bigger than the cost",
          homo.mrca_intruder_count > homo.conflict_cost, True)

    print("\nCost is the true minimum, not min(intruders, stray members)")
    # The regression that would have caught shipping min(a, b): one foreign
    # tip sits inside the group's core and one member is parked outside, so
    # the cheapest repair mixes both moves. Evicting everything under the
    # MRCA costs 4; gathering every stray member costs 3; the answer is 2.
    _out, by_name = groups_of(
        "((Homo_sapiens,Homo_erectus,Rosa_canina,Homo_habilis),"
        "Drosophila_melanogaster,Drosophila_simulans,Drosophila_yakuba,"
        "Homo_naledi);")
    homo = by_name["Homo"]
    largest = max(
        len([t for t in c.get_tip_labels() if t in set(homo.members)])
        for c in homo.clusters)
    stray = homo.member_count - largest
    check("evict one and gather one", homo.conflict_cost, 2)
    check("type", homo.conflict_type, "both")
    check("no single culprit", homo.culprit, None)
    check("beats evicting everything under the MRCA",
          homo.conflict_cost < homo.mrca_intruder_count, True)
    check("beats gathering every stray member",
          homo.conflict_cost < stray, True)
    check("and beats the cheaper of the two",
          homo.conflict_cost < min(homo.mrca_intruder_count, stray), True)

    print("\nUnplaced tips are reported, never blamed")
    _out, by_name = groups_of(
        "(((Homo_sapiens,Xyzzy_qqq),Homo_erectus),"
        "(Drosophila_melanogaster,Drosophila_simulans));")
    homo = by_name["Homo"]
    check("still monophyletic", homo.status, "monophyletic")
    check("not called an intruder", homo.intruders, [])
    check("reported separately", homo.unknown_inside, ["Xyzzy_qqq"])
    check("costs nothing", homo.conflict_cost, 0)
    check("one cluster", homo.cluster_count, 1)

    print("\nInvariants across a battery of trees")
    trees = [
        "(((Homo_sapiens,Homo_erectus),(Drosophila_melanogaster,"
        "Drosophila_simulans)),Rosa_canina);",
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));",
        "(Homo_sapiens,Homo_erectus,Drosophila_melanogaster,"
        "Drosophila_simulans);",
        "((((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "Rosa_canina),Homo_habilis);",
    ]
    bad_cluster = bad_cost = bad_sum = 0
    for newick in trees:
        out, _ = groups_of(newick)
        for r in out["taxon_results"]:
            if (r.cluster_count == 1) != r.is_monophyletic:
                bad_cluster += 1
            if (r.conflict_cost == 0) != r.is_monophyletic:
                bad_cost += 1
            if r.conflict_cost != len(r.intruders) + len(r.outliers):
                bad_sum += 1
    check("clusters == 1 iff monophyletic", bad_cluster, 0)
    check("cost == 0 iff monophyletic", bad_cost, 0)
    check("cost == intruders + outliers", bad_sum, 0)

    print("\nPolytomies")
    _out, by_name = groups_of(
        "(Homo_sapiens,Homo_erectus,Drosophila_melanogaster,"
        "Drosophila_simulans);")
    check("star tree splits Homo", by_name["Homo"].cluster_count, 2)
    check("one move", by_name["Homo"].conflict_cost, 1)

    print("\nDuplicate tip labels are counted by node")
    out, by_name = groups_of(
        "(((Homo_sapiens,Homo_sapiens),Drosophila_melanogaster),"
        "Drosophila_simulans);")
    hs = by_name["Homo sapiens"]
    check("two member nodes", hs.member_count, 2)
    check("monophyletic", hs.status, "monophyletic")


def test_collapsing_and_culprits(db):
    checker = MonophylyChecker(db)
    out = checker.analyze(parse_newick(
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));"))

    print("\nGroups covering the same tips merge into one row")
    # Homo and Homo sapiens+erectus differ, but Rosa/Rosales/Plantae all cover
    # the single Rosa tip, and Drosophila/Arthropoda cover the same two flies.
    for group in out["groups"]:
        names = {t.name for t in group.taxa}
        if {"Drosophila", "Arthropoda"} <= names:
            check("Drosophila and Arthropoda merged", group.is_collapsed, True)
            check("ordered coarse to fine",
                  [t.name for t in group.taxa][0], "Arthropoda")
            check("label spans the chain", " … " in group.name, True)
            check("rank range", "–" in group.rank, True)
            break
    else:
        check("a collapsed group exists", False, True)

    print("\nEvery collapsed row shares one analysis")
    bad = 0
    for group in out["groups"]:
        for taxon in group.taxa:
            if (taxon.conflict_cost != group.result.conflict_cost
                    or set(taxon.members) != set(group.result.members)):
                bad += 1
    check("constituents agree with the representative", bad, 0)

    print("\nTips are ranked by what removing them would repair")
    ranked = [t for t in out["conflicting_tips"] if t.repairs]
    check("ranked most-repairs first",
          [len(t.repairs) for t in ranked],
          sorted((len(t.repairs) for t in ranked), reverse=True))
    check("every listed group is really repaired by that tip",
          all(g.culprit == t.label for t in ranked for g in t.repairs), True)
    check("the nested fly is credited for the group it intrudes into",
          any(t.label == "Drosophila_melanogaster"
              and any("Homo" in g.name for g in t.repairs) for t in ranked),
          True)
    check("role recorded", ranked[0].role in ("intruder", "member", "both"),
          True)

    print("\nA sprawling clade credits nobody")
    # Proteales-style: one member elsewhere pulls in many tips, but removing
    # any one of those tips repairs nothing.
    out = checker.analyze(parse_newick(
        "((Homo_sapiens,(Drosophila_melanogaster,Drosophila_simulans,"
        "Drosophila_yakuba)),Homo_erectus);"))
    ranked = {t.label: len(t.repairs) for t in out["conflicting_tips"]}
    flies = [n for n in ranked if n.startswith("Drosophila")]
    check("no fly is credited with repairing Homo",
          all(ranked[f] == 0 for f in flies), True)


def test_pruning():
    print("\nPruning removes tips without distorting what is left")
    from tree_parser import copy_tree, prune_tips

    def pruned(newick, drop):
        root, removed, missing = prune_tips(parse_newick(newick), drop)
        text = root.to_newick(True, "support") + ";" if root else None
        return text, removed, missing

    check("unary parent is spliced out",
          pruned("((a:1,b:2)90:3,(c:4,d:5)80:6);", ["a"])[0],
          "(b:5,(c:4,d:5)80:6);")
    check("emptying a clade removes it and collapses the root",
          pruned("((a:1,b:2)90:3,(c:4,d:5)80:6);", ["a", "b"])[0],
          "(c:4,d:5);")
    check("a unary root collapses",
          pruned("(((a,b),c),d);", ["d"])[0], "((a,b),c);")
    check("a polytomy stays a polytomy",
          pruned("(a,b,c,d);", ["a"])[0], "(b,c,d);")
    check("refuses to leave fewer than two tips",
          pruned("((a,b),(c,d));", ["a", "b", "c"])[0], None)
    check("labels not in the tree are reported, not fatal",
          pruned("((a,b),(c,d));", ["zzz"])[2], ["zzz"])
    check("every tip carrying a duplicated label goes",
          pruned("(((a,b),a),c);", ["a"])[1], ["a", "a"])

    original = parse_newick("((a:0.1,b:0.2)95:0.3,(c:0.4,d:0.5)88:0.6)100;")
    duplicate = copy_tree(original)
    check("a copy is identical",
          duplicate.to_newick(True, "support"),
          original.to_newick(True, "support"))
    prune_tips(duplicate, ["a"])
    check("pruning the copy leaves the original alone",
          original.count_tips(), 4)

    sys.setrecursionlimit(1000)
    ladder = "t0"
    for index in range(1, 4000):
        ladder = f"({ladder},t{index})"
    deep = copy_tree(parse_newick(ladder + ";"))
    check("a 4000-tip ladder copies without recursing", deep.count_tips(), 4000)
    check("and prunes", prune_tips(deep, ["t0"])[0].count_tips(), 3999)


def test_removal_simulation(db):
    checker = MonophylyChecker(db)

    def simulate(newick, drop):
        tree = parse_newick(newick)
        out = checker.analyze(tree)
        sim = checker.simulate_removal(
            tree, out["tip_lineages"], out["unresolved"], drop,
            out["taxon_results"])
        return tree, sim, {c.name: c for c in sim.changes}

    print("\nRemoving an intruder is confirmed by re-analysis")
    tree, sim, by_name = simulate(
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));", ["Drosophila_melanogaster"])
    check("Homo comes out monophyletic", by_name["Homo"].outcome, "resolved")
    check("fewer broken groups than before",
          sim.broken_after < sim.broken_before, True)
    check("nothing became broken", sim.newly_broken, [])
    check("the analysed tree is untouched", tree.count_tips(), 5)

    print("\nA group that loses its members is reported as dropped, not fixed")
    check("Drosophila is not called resolved",
          by_name["Drosophila"].outcome, "dropped_small")
    check("and its cost after is blank rather than zero",
          by_name["Drosophila"].cost_after, None)

    print("\nRemoving an unmatched tip provably changes nothing")
    tree, sim, by_name = simulate(
        "(((Homo_sapiens,Xyzzy_qqq),Homo_erectus),"
        "(Drosophila_melanogaster,Drosophila_simulans));", ["Xyzzy_qqq"])
    check("flagged as a no-op", sim.only_unresolved, True)
    check("no group changed cost",
          [c.name for c in sim.changes
           if c.cost_after is not None and c.cost_after != c.cost_before], [])

    print("\nRefusals")
    tree = parse_newick("((Homo_sapiens,Homo_erectus),Drosophila_melanogaster);")
    out = checker.analyze(tree)
    sim = checker.simulate_removal(tree, out["tip_lineages"],
                                   out["unresolved"], ["Not_a_tip"],
                                   out["taxon_results"])
    check("unknown labels are refused with a reason",
          bool(sim.refused) and sim.missing == ["Not_a_tip"], True)
    sim = checker.simulate_removal(
        tree, out["tip_lineages"], out["unresolved"],
        ["Homo_sapiens", "Homo_erectus"], out["taxon_results"])
    check("refuses to leave fewer than two tips", bool(sim.refused), True)

    print("\nInvariants over a battery of trees")
    trees = [
        "(((Homo_sapiens,Drosophila_melanogaster),Homo_erectus),"
        "(Drosophila_simulans,Rosa_canina));",
        "((Homo_sapiens,(Drosophila_melanogaster,Drosophila_simulans)),"
        "Homo_erectus);",
        "((Homo_sapiens,Homo_erectus,Rosa_canina,Homo_habilis),"
        "Drosophila_melanogaster,Drosophila_simulans,Drosophila_yakuba,"
        "Homo_naledi);",
    ]
    rose = broke = 0
    for newick in trees:
        tree = parse_newick(newick)
        out = checker.analyze(tree)
        for victim in sorted({t.label for t in tree.get_tip_nodes()}):
            sim = checker.simulate_removal(
                tree, out["tip_lineages"], out["unresolved"], [victim],
                out["taxon_results"])
            if sim.refused:
                continue
            broke += len(sim.newly_broken)
            rose += sum(1 for c in sim.changes
                        if c.cost_after is not None
                        and c.cost_after > c.cost_before)
    check("removing a tip never raises a group's cost", rose, 0)
    check("removing a tip never breaks an intact group", broke, 0)


def main():
    db = build_db()
    test_name_resolution(db)
    test_taxid_tips(db)
    test_parsing()
    test_placeholder_taxa()
    test_weak_branches(db)
    test_index_tips()
    test_monophyly(db)
    test_conflict_costs(db)
    test_collapsing_and_culprits(db)
    test_pruning()
    test_removal_simulation(db)
    print(f"\n{_passed} passed, {_failed} failed")
    return 1 if _failed else 0


if __name__ == "__main__":
    sys.exit(main())
