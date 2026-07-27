"""
monophyly.py - Monophyly Checking Engine

Checks whether taxonomic groups are monophyletic in a gene tree with respect
to the NCBI taxonomy. Inspired by MonoPhy (Schwery & O'Meara, 2016).

Two views of the same data:

1. Taxon-level (check_taxa) - the question a biologist actually asks:
   "is Pentapetalae monophyletic in this tree, and if not, what is inside it
   that shouldn't be?" One row per named group.

2. Node-level (check_tree) - one result per internal node, used to colour the
   tree drawing and to write taxonomy labels into the output Newick.

The node-level view alone is misleading as a report: a group that is broken
once shows up at every nested node beneath it, so a 232-tip tree produces
four consecutive "Mesangiospermae" rows describing a single conflict.

Pure Python stdlib - no external dependencies.
"""

from tree_parser import clean_newick_name, index_tips

# Ranks worth reporting on by default, coarse to fine. NCBI also emits
# "clade" and "no rank" for a great many nodes.
MAIN_RANKS = ("superkingdom", "kingdom", "phylum", "class", "order",
              "family", "genus", "species")

_RANK_ORDER = {rank: i for i, rank in enumerate(MAIN_RANKS)}


def plural(count, word):
    """Render "1 tip" / "2 tips" so reports read like prose."""
    return f"{count} {word}" if count == 1 else f"{count} {word}s"


class MonophylyResult:
    """Result of monophyly analysis for a single tree node."""

    __slots__ = ("node", "tip_labels", "mrca_code", "mrca_name", "mrca_rank",
                 "status", "candidate_names", "intruders")

    def __init__(self, node, tip_labels):
        self.node = node
        self.tip_labels = tip_labels
        self.mrca_code = None
        self.mrca_name = ""
        self.mrca_rank = ""
        self.status = "unknown"  # "monophyletic", "not_monophyletic", "unresolved"
        self.candidate_names = []
        self.intruders = []  # Taxa outside clade that share the MRCA

    def __repr__(self):
        return (f"MonophylyResult(name={self.mrca_name}, rank={self.mrca_rank}, "
                f"status={self.status}, tips={len(self.tip_labels)})")


class TaxonResult:
    """Whether one named taxonomic group is monophyletic in the tree.

    `home_node` is the clade it is cheapest to call "the group". The tips
    there split three ways, and the members split two:

      members inside    tips that belong here and are here
      intruders         resolved tips inside that belong elsewhere
      unknown_inside    tips inside that matched no taxon at all
      outliers          members of the group sitting outside home_node

    Unknown tips are held separately and never counted as intruders. A tip we
    could not place is not evidence that the group is broken, and blaming it
    would let one unmatched sequence appear to break every group it sits in.

    `clusters` are the maximal clades containing only members (and unknowns):
    exactly one cluster means monophyletic. `conflict_cost` is the fewest tips
    that would have to move to make the group monophyletic.
    """

    __slots__ = ("tax_id", "name", "rank", "members", "member_count",
                 "node", "home_node", "clusters", "intruders", "outliers",
                 "unknown_inside", "mrca_intruder_count", "status",
                 "_clade_tips")

    def __init__(self, tax_id, name, rank, members, member_count=None):
        self.tax_id = tax_id
        self.name = name
        self.rank = rank
        self.members = members
        self.member_count = len(members) if member_count is None else member_count
        self.node = None          # MRCA of the members
        self.home_node = None     # cheapest clade to call "the group"
        self.clusters = []        # maximal clades holding only members
        self.intruders = []       # foreign resolved tips inside home_node
        self.outliers = []        # members sitting outside home_node
        self.unknown_inside = []  # unplaced tips inside home_node
        self.mrca_intruder_count = 0
        self.status = "monophyletic"
        self._clade_tips = None

    @property
    def is_monophyletic(self):
        return self.status == "monophyletic"

    @property
    def clade_tips(self):
        """Tips under the group's home clade, computed on demand."""
        if self._clade_tips is None:
            self._clade_tips = (self.home_node.get_tip_labels()
                                if self.home_node else [])
        return self._clade_tips

    @property
    def conflict_cost(self):
        """Fewest tips that must move to make this group a clade.

        Evaluated over every candidate home clade, not just two of them.
        Taking `min(intruders under the MRCA, members outside the largest
        cluster)` looks equivalent but is only an upper bound, and a loose
        one: it rates rosids at 59 when evicting one tip and gathering one
        other makes it a clade.
        """
        return len(self.intruders) + len(self.outliers)

    @property
    def conflict_type(self):
        """Whether foreign tips are inside, members are outside, or both."""
        if self.is_monophyletic:
            return ""
        if not self.outliers:
            return "intruded"
        if not self.intruders:
            return "scattered"
        return "both"

    @property
    def culprit(self):
        """The one tip whose removal would repair this group, if there is one.

        A group needs exactly one move precisely when a single tip is at
        fault, so this needs no search.
        """
        if self.conflict_cost != 1:
            return None
        return (self.intruders or self.outliers)[0]

    @property
    def cluster_count(self):
        """How fragmented the members are.

        Fragmentation, NOT severity: one intruder placed deep in a ladder
        splits Fabaceae's 25 members into 8 pure clades while still being a
        one-tip fix. Never sort on this.
        """
        return len(self.clusters)

    def severity_sort_key(self):
        """Easiest to fix first, then the biggest groups."""
        return (self.conflict_cost, -self.member_count, self.name)

    def rank_sort_key(self):
        """Sort coarse ranks first, unranked groups last."""
        return (_RANK_ORDER.get(self.rank, len(MAIN_RANKS)), -self.member_count)

    def __repr__(self):
        return (f"TaxonResult({self.name}, {self.rank}, "
                f"{self.member_count} members, {self.status})")


class CollapsedGroup:
    """Taxa that share an identical set of sampled tips.

    Apiales, Apiineae, Apiaceae and Apioideae contain exactly the same four
    tips in the test tree, so the tree cannot tell them apart and reporting
    them as four findings is four ways of saying one thing.

    `taxa` is ordered coarse to fine. Because they share a member set they
    also share an MRCA, so every analysis field is identical and `result`
    (the coarsest) speaks for all of them.
    """

    __slots__ = ("taxa", "result")

    def __init__(self, taxa):
        self.taxa = taxa
        self.result = taxa[0]

    @property
    def is_collapsed(self):
        return len(self.taxa) > 1

    @property
    def name(self):
        if not self.is_collapsed:
            return self.result.name
        return f"{self.taxa[0].name} … {self.taxa[-1].name}"

    @property
    def rank(self):
        first = self.taxa[0].rank or "no rank"
        if not self.is_collapsed:
            return first
        last = self.taxa[-1].rank or "no rank"
        return first if first == last else f"{first}–{last}"

    def __getattr__(self, item):
        # Delegate the analysis fields to the representative taxon.
        return getattr(self.result, item)

    def __repr__(self):
        return f"CollapsedGroup({self.name}, {len(self.taxa)} taxa)"


class ConflictingTip:
    """A tip ranked by how many groups its removal would repair."""

    __slots__ = ("label", "repairs", "appears_in", "role")

    def __init__(self, label):
        self.label = label
        self.repairs = []      # groups this tip alone would fix
        self.appears_in = []   # groups it intrudes into, repairable or not
        self.role = ""         # "intruder", "member", or "both"

    def sort_key(self):
        return (-len(self.repairs), -len(self.appears_in), self.label)

    def __repr__(self):
        return (f"ConflictingTip({self.label}, repairs "
                f"{len(self.repairs)})")


class MonophylyChecker:
    """Check monophyly of tree clades against NCBI taxonomy."""

    def __init__(self, taxonomy_db):
        self.db = taxonomy_db

    # ---- Tip resolution ----

    def resolve_tips(self, root, progress_callback=None):
        """Resolve every tip label to a taxonomy lineage.

        Returns (tip_lineages, matches, unresolved) where tip_lineages maps
        tip label -> list of tax_id codes from tip to root, matches maps tip
        label -> TipMatch, and unresolved lists labels with no match at all.
        """
        all_tips = root.get_tip_labels()
        tip_lineages = {}
        matches = {}
        unresolved = []

        for tip in all_tips:
            match = self.db.resolve_tip(tip)
            lineage = self.db.get_lineage_codes(match.tax_id) if match else []
            if match and lineage:
                matches[tip] = match
                tip_lineages[tip] = lineage
            else:
                # A name that matched but has no lineage is just as unusable
                # as one that never matched - report both rather than
                # dropping the first kind silently.
                unresolved.append(tip)
                tip_lineages[tip] = []

        if progress_callback:
            resolved = len(all_tips) - len(unresolved)
            inexact = sum(1 for m in matches.values() if not m.is_exact)
            progress_callback(
                f"Resolved {resolved}/{len(all_tips)} tip names "
                f"({inexact} inexact, {len(unresolved)} unresolved)"
            )
        return tip_lineages, matches, unresolved

    # ---- Taxon-level analysis (the report biologists read) ----

    def check_taxa(self, root, tip_lineages, unresolved=(), min_members=2,
                   progress_callback=None):
        """Test each taxonomic group represented in the tree for monophyly.

        Groups with fewer than min_members sampled tips are skipped: a group
        with one representative is monophyletic by definition and says
        nothing. A group covering every resolved tip is skipped for the same
        reason.
        """
        tip_nodes, spans = index_tips(root)
        position = {id(node): i for i, node in enumerate(tip_nodes)}

        members_by_taxon = {}
        total_resolved = 0
        for node in tip_nodes:
            lineage = tip_lineages.get(node.label) or ()
            if lineage:
                total_resolved += 1
            for code in lineage:
                members_by_taxon.setdefault(code, []).append(node)

        clusters, homes, mrcas = self._taxon_pass(root, tip_lineages)

        unresolved_set = set(unresolved)
        results = []
        for code, member_nodes in members_by_taxon.items():
            if not min_members <= len(member_nodes) < total_resolved:
                continue

            result = TaxonResult(
                tax_id=code,
                name=self.db.get_scientific_name(code) or code,
                rank=self.db.code_to_rank.get(code, ""),
                members=[n.label for n in member_nodes],
                member_count=len(member_nodes),
            )
            result.clusters = [node for node, _size in clusters.get(code, ())]

            home, inside_members, inside_resolved = homes[code]
            result.home_node = home
            mrca_node, mrca_resolved = mrcas.get(code, (home, inside_resolved))
            result.node = mrca_node
            result.mrca_intruder_count = mrca_resolved - len(member_nodes)

            lo, hi, span_tips = spans[id(home)]
            result.outliers = [n.label for n in member_nodes
                               if not lo <= position[id(n)] <= hi]

            # The counts come free from the pass; only walk the clade when
            # there is actually something in it to name.
            if inside_resolved > inside_members or span_tips > inside_resolved:
                member_labels = {n.label for n in member_nodes}
                for tip in home.get_tip_labels():
                    if tip in member_labels:
                        continue
                    if tip in unresolved_set:
                        result.unknown_inside.append(tip)
                    else:
                        result.intruders.append(tip)

            result.status = ("monophyletic" if result.conflict_cost == 0
                             else "not_monophyletic")
            results.append(result)

        results.sort(key=TaxonResult.rank_sort_key)

        if progress_callback:
            broken = sum(1 for r in results if not r.is_monophyletic)
            progress_callback(
                f"Checked {len(results)} taxonomic groups: "
                f"{len(results) - broken} monophyletic, {broken} broken"
            )
        return results

    @staticmethod
    def _taxon_pass(root, tip_lineages):
        """One bottom-up pass that scores every taxon at every node.

        For each node this carries the number of resolved tips below it and,
        per taxon, how many of those tips are members. From that, three
        things fall out at once for every taxon in the tree:

        * clusters - a child is a maximal all-member clade when the child is
          pure and its parent is not
        * home clade - adopting node v as "the group" costs
          `foreign_inside + members_outside`, which is
          `M + resolved[v] - 2 * members[v]`; since M is fixed per taxon,
          minimising `resolved[v] - 2 * members[v]` finds the cheapest clade
          over the whole tree
        * MRCA - the first node holding every member

        Unresolved tips are counted in neither total, so they are wildcards:
        they never spoil a cluster and never become intruders. That is what
        stops one unmatched sequence appearing to break every group it sits
        in.

        Cost is O(nodes x lineage depth) for the whole tree rather than
        O(taxa x tips), which matters because a tree has far more taxa in
        play than it has ranks.
        """
        counts = {}     # id(node) -> {code: member tips below}
        resolved = {}   # id(node) -> resolved tips below
        clusters = {}   # code -> [(node, size)]
        homes = {}      # code -> (node, members inside, resolved inside)
        home_key = {}
        mrcas = {}

        totals = {}
        for node in root.get_tip_nodes():
            for code in tip_lineages.get(node.label) or ():
                totals[code] = totals.get(code, 0) + 1

        for node in root.post_order():
            nid = id(node)
            if node.is_tip:
                lineage = tip_lineages.get(node.label) or ()
                counts[nid] = dict.fromkeys(lineage, 1)
                resolved[nid] = 1 if lineage else 0
            else:
                merged = {}
                total = 0
                for child in node.children:
                    total += resolved[id(child)]
                    for code, n in counts[id(child)].items():
                        merged[code] = merged.get(code, 0) + n
                # A pure child whose parent is impure is a maximal cluster.
                for child in node.children:
                    child_resolved = resolved.pop(id(child))
                    child_counts = counts.pop(id(child))
                    if not child_resolved:
                        continue
                    for code, n in child_counts.items():
                        if n == child_resolved and merged[code] != total:
                            clusters.setdefault(code, []).append(
                                (child, n))
                counts[nid] = merged
                resolved[nid] = total

            # This loop runs once per (node, taxon) pair - about 91 times per
            # tip on a real tree - so the lookups are hoisted deliberately.
            here = resolved[nid]
            root_rank = 1 if node.parent is None else 0
            best_get = home_key.get
            totals_get = totals.get
            for code, m in counts[nid].items():
                # Prefer cheapest; then the reading that keeps the most
                # members together; then a real clade over the whole tree.
                #
                # Order matters. Several clades routinely tie on cost - in
                # ((a1,b1),a2) both "move b1 out" and "move a2 across" cost
                # one - and preferring the clade holding more members is what
                # names the foreign tip rather than splitting the group and
                # blaming one of its own.
                key = (here - 2 * m, -m, root_rank)
                previous = best_get(code)
                if previous is None or key < previous:
                    home_key[code] = key
                    homes[code] = (node, m, here)
                if m == totals_get(code) and code not in mrcas:
                    mrcas[code] = (node, here)

        # Anything still pure when we reach the root is a cluster at the root.
        root_resolved = resolved[id(root)]
        for code, m in counts[id(root)].items():
            if m == root_resolved:
                clusters.setdefault(code, []).append((root, m))

        return clusters, homes, mrcas

    # ---- Node-level analysis (drives tree colouring and labels) ----

    def check_tree(self, root, progress_callback=None):
        """Run node-level monophyly analysis on an entire tree.

        Returns (results, unresolved) for backwards compatibility. Use
        analyze() when the taxon-level report is wanted too.
        """
        tip_lineages, _matches, unresolved = self.resolve_tips(
            root, progress_callback)
        results = self._check_nodes(root, tip_lineages, progress_callback)
        return results, unresolved

    def analyze(self, root, progress_callback=None):
        """Full analysis: resolve tips, check nodes, and check taxa.

        Returns a dict with keys: node_results, taxon_results, matches,
        unresolved, tip_lineages.
        """
        tip_lineages, matches, unresolved = self.resolve_tips(
            root, progress_callback)
        node_results = self._check_nodes(root, tip_lineages, progress_callback)
        taxon_results = self.check_taxa(root, tip_lineages, unresolved,
                                        progress_callback=progress_callback)
        # Rank tips against the collapsed rows so the two tables agree: the
        # four Apiales-chain taxa are one finding, and should count once.
        groups = self.collapse_groups(taxon_results)
        return {
            "node_results": node_results,
            "taxon_results": taxon_results,
            "groups": groups,
            "conflicting_tips": rank_conflicting_tips(groups),
            "matches": matches,
            "unresolved": unresolved,
            "tip_lineages": tip_lineages,
        }

    def collapse_groups(self, taxon_results):
        """Merge taxa that cover an identical set of sampled tips.

        Constituents are ordered by taxonomic depth rather than by rank name.
        Within a shared member set the taxa are always strictly nested, so
        lineage length is a total order - and it sidesteps the fact that most
        of the ranks involved ("clade", "no rank", "subtribe", ...) have no
        position in any rank ordering.
        """
        by_members = {}
        for result in taxon_results:
            by_members.setdefault(frozenset(result.members), []).append(result)

        groups = []
        for taxa in by_members.values():
            if len(taxa) > 1:
                taxa = sorted(
                    taxa,
                    key=lambda r: (len(self.db.get_lineage_codes(r.tax_id)),
                                   r.name))
            groups.append(CollapsedGroup(taxa))
        return groups

    def _check_nodes(self, root, tip_lineages, progress_callback=None):
        """Check every internal node's bipartition against the taxonomy."""
        all_tips = root.get_tip_labels()
        bipartitions = root.get_bipartitions()

        results = []
        for idx, (bip_tips, node) in enumerate(bipartitions):
            if progress_callback and idx % 50 == 0:
                progress_callback(
                    f"Checking node {idx + 1}/{len(bipartitions)}...")

            result = self._check_bipartition(
                bip_tips, all_tips, tip_lineages, node)
            results.append(result)

            node.monophyly_status = result.status
            node.taxonomy_code = result.mrca_code
            node.taxonomy_name = result.mrca_name

        if progress_callback:
            mono = sum(1 for r in results if r.status == "monophyletic")
            non_mono = sum(1 for r in results
                           if r.status == "not_monophyletic")
            progress_callback(
                f"Done: {mono} monophyletic, {non_mono} non-monophyletic, "
                f"{len(results) - mono - non_mono} unresolved"
            )
        return results

    def _check_bipartition(self, bip_tips, all_tips, tip_lineages, node):
        """Check if a single bipartition is monophyletic."""
        result = MonophylyResult(node, bip_tips)

        # Build the membership set once, not once per tip of the tree.
        inside = set(bip_tips)
        complement = [t for t in all_tips if t not in inside]

        clade_lineages = [tip_lineages[t] for t in bip_tips if tip_lineages[t]]
        if not clade_lineages:
            result.status = "unresolved"
            return result

        mrca = self._find_mrca(clade_lineages)
        if not mrca:
            result.status = "unresolved"
            return result

        result.mrca_code = mrca
        result.mrca_name = self.db.get_scientific_name(mrca)
        result.mrca_rank = self.db.code_to_rank.get(mrca, "")
        result.candidate_names = self.db.get_names_for_code(mrca)

        for tip in complement:
            if tip_lineages[tip] and mrca in tip_lineages[tip]:
                result.intruders.append(tip)

        result.status = ("not_monophyletic" if result.intruders
                         else "monophyletic")
        return result

    def _find_mrca(self, lineage_list):
        """Find the most recent common ancestor from a list of lineages.

        Each lineage is a list of tax_id codes from tip to root.
        Returns the first code in lineage[0] that appears in ALL other lineages.
        """
        if not lineage_list:
            return None
        if len(lineage_list) == 1:
            return lineage_list[0][0] if lineage_list[0] else None

        other_sets = [set(lin) for lin in lineage_list[1:]]
        for code in lineage_list[0]:
            if all(code in s for s in other_sets):
                return code
        return None

    # ---- Labelling ----

    def label_tree(self, root, results, auto_select=True):
        """Apply taxonomy labels from node results to the tree.

        Only the node label is set; any support value the tree arrived with
        stays in node.support and is written back out on export.
        """
        for result in results:
            if result.status == "monophyletic" and result.mrca_name:
                result.node.label = clean_newick_name(result.mrca_name)

    # ---- Reporting ----

    def get_summary(self, results, taxon_results=None, matches=None):
        """Generate a text report.

        When taxon_results are supplied the report leads with the taxon-level
        view, which has one row per named group instead of one per node.
        """
        lines = []

        if taxon_results is not None:
            groups = self.collapse_groups(taxon_results)
            broken = [g for g in groups if not g.is_monophyletic]
            intact = [g for g in groups if g.is_monophyletic]
            collapsed = sum(len(g.taxa) for g in groups if g.is_collapsed)

            repairing = [t for t in rank_conflicting_tips(groups) if t.repairs]
            if repairing:
                fixable = sum(len(t.repairs) for t in repairing)
                lines.append("=== What To Check First ===")
                lines.append(
                    f"Rechecking {plural(len(repairing), 'sequence')} would "
                    f"repair {fixable} of {len(broken)} conflicting groups.")
                lines.append("")
                for tip in repairing[:10]:
                    names = ", ".join(g.name for g in tip.repairs[:6])
                    if len(tip.repairs) > 6:
                        names += f"... (+{len(tip.repairs) - 6} more)"
                    lines.append(
                        f"  {tip.label}: repairs "
                        f"{plural(len(tip.repairs), 'group')} as {tip.role}")
                    lines.append(f"      {names}")
                if len(repairing) > 10:
                    lines.append(f"  ... {len(repairing) - 10} more")
                lines.append("")

            lines.append("=== Taxonomic Groups ===")
            lines.append(f"Groups with 2+ sampled tips: {len(taxon_results)}")
            if collapsed:
                lines.append(
                    f"Rows shown: {len(groups)} "
                    f"({collapsed} taxa share a member set and are merged)")
            lines.append(f"Monophyletic: {len(intact)}")
            lines.append(f"Not monophyletic: {len(broken)}")
            lines.append("")

            if broken:
                lines.append("--- Groups Broken In This Tree ---")
                lines.append("Cost = fewest tips that would have to move.")
                lines.append("")
                for g in sorted(broken, key=lambda x: x.severity_sort_key()):
                    lines.append(
                        f"  {g.name} ({g.rank}) "
                        f"[{g.member_count} members]")
                    lines.append(
                        f"      fix: move {plural(g.conflict_cost, 'tip')} "
                        f"({g.conflict_type})")
                    if g.culprit:
                        lines.append(f"      blame: {g.culprit}")
                    if g.intruders:
                        listed = ", ".join(g.intruders[:5])
                        if len(g.intruders) > 5:
                            listed += f"... (+{len(g.intruders) - 5} more)"
                        lines.append(
                            f"      {plural(len(g.intruders), 'tip')} "
                            f"inside: {listed}")
                    if g.outliers:
                        listed = ", ".join(g.outliers[:5])
                        if len(g.outliers) > 5:
                            listed += f"... (+{len(g.outliers) - 5} more)"
                        lines.append(
                            f"      {plural(len(g.outliers), 'member')} "
                            f"outside: {listed}")
                    if g.is_collapsed:
                        lines.append(
                            "      same tips as: "
                            + ", ".join(f"{t.name} ({t.rank or 'no rank'})"
                                        for t in g.taxa))
                    if g.unknown_inside:
                        lines.append(
                            f"      unplaced tips inside: "
                            f"{len(g.unknown_inside)} (not counted as "
                            f"intruders)")
                lines.append("")

            if intact:
                lines.append("--- Monophyletic Groups ---")
                for g in sorted(intact, key=lambda x: x.rank_sort_key()):
                    note = ""
                    if g.is_collapsed:
                        note = f"  [+{len(g.taxa) - 1} nested names]"
                    lines.append(
                        f"  {g.name} ({g.rank}) "
                        f"[{g.member_count} tips]{note}")
                lines.append("")

        if matches:
            inexact = [m for m in matches.values() if not m.is_exact]
            ambiguous = [m for m in matches.values() if m.ambiguous]
            if inexact or ambiguous:
                lines.append("--- Tip Matching Warnings ---")
                for m in sorted(inexact, key=lambda x: x.label):
                    lines.append(
                        f"  {m.label} -> {m.matched_name} "
                        f"({m.rank or 'no rank'}, matched by {m.method})")
                for m in sorted(ambiguous, key=lambda x: x.label):
                    lines.append(
                        f"  {m.label} -> ambiguous name, several taxa share it")
                lines.append("")

        mono = [r for r in results if r.status == "monophyletic"]
        non_mono = [r for r in results if r.status == "not_monophyletic"]
        unresolved = [r for r in results if r.status == "unresolved"]

        lines.append("=== Node-Level Detail ===")
        lines.append(f"Total internal nodes: {len(results)}")
        lines.append(f"Monophyletic: {len(mono)}")
        lines.append(f"Non-monophyletic: {len(non_mono)}")
        lines.append(f"Unresolved: {len(unresolved)}")

        return "\n".join(lines)


def rank_conflicting_tips(taxon_results):
    """Rank tips by how many broken groups removing them would repair.

    Ranking on "how many broken groups is this tip inside" does not work: a
    group whose members are scattered drags a huge enclosing clade along, so
    all 159 tips inside Proteales would score alike. Asking instead whether
    removing the tip *repairs* the group is immune to that by construction -
    deleting one of those 159 leaves 158, so such a group credits nobody.

    It also catches tips that break a group by being a stray member rather
    than an intruder, which an intruder-only count cannot see.
    """
    tips = {}

    def entry(label):
        tip = tips.get(label)
        if tip is None:
            tip = tips[label] = ConflictingTip(label)
        return tip

    for result in taxon_results:
        if result.is_monophyletic:
            continue
        for label in result.intruders:
            entry(label).appears_in.append(result)

        culprit = result.culprit
        if culprit is None:
            continue
        tip = entry(culprit)
        tip.repairs.append(result)
        role = "intruder" if result.intruders else "member"
        tip.role = role if tip.role in ("", role) else "both"

    return sorted(tips.values(), key=ConflictingTip.sort_key)
