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

from tree_parser import clean_newick_name, copy_tree, index_tips, prune_tips

# Ranks worth reporting on by default, coarse to fine. NCBI also emits
# "clade" and "no rank" for a great many nodes.
MAIN_RANKS = ("superkingdom", "kingdom", "phylum", "class", "order",
              "family", "genus", "species")

_RANK_ORDER = {rank: i for i, rank in enumerate(MAIN_RANKS)}

# How many equally-cheap repairs to keep per group. Enough to show a
# biologist the real ambiguity, bounded so a pathological tree cannot
# blow up memory.
MAX_REPAIR_OPTIONS = 12


def plural(count, word):
    """Render "1 tip" / "2 tips" so reports read like prose."""
    return f"{count} {word}" if count == 1 else f"{count} {word}s"


class MonophylyResult:
    """Result of monophyly analysis for a single tree node.

    Deliberately counts, not lists. This exists to colour the drawing and to
    supply the taxonomy labels written into the exported Newick, and there is
    one of these per internal node - materialising the clade's tips and its
    intruders at every node is what made the node pass quadratic. The
    taxon-level `TaxonResult` is where named tips live.
    """

    __slots__ = ("node", "tip_count", "mrca_code", "mrca_name", "mrca_rank",
                 "status", "intruder_count")

    def __init__(self, node, tip_count=0):
        self.node = node
        self.tip_count = tip_count          # resolved tips below this node
        self.mrca_code = None
        self.mrca_name = ""
        self.mrca_rank = ""
        self.status = "unknown"  # "monophyletic", "not_monophyletic", "unresolved"
        self.intruder_count = 0  # tips elsewhere sharing this node's MRCA

    def __repr__(self):
        return (f"MonophylyResult(name={self.mrca_name}, rank={self.mrca_rank}, "
                f"status={self.status}, tips={self.tip_count})")


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
                 "unknown_inside", "mrca_intruder_count", "repair_options",
                 "status", "_clade_tips")

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
        self.repair_options = []  # tips that alone would repair the group
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
        """One tip whose removal would repair this group, if any would.

        Only one of possibly several - see `repair_options` for the rest and
        `is_ambiguous` to know whether this is a choice or the answer.
        """
        return self.repair_options[0] if self.repair_options else None

    @property
    def is_ambiguous(self):
        """True when more than one single tip would repair this group."""
        return len(self.repair_options) > 1

    def culprit_summary(self):
        """Short form for a table cell: the tip, plus how many others qualify."""
        if not self.repair_options:
            return ""
        others = len(self.repair_options) - 1
        if not others:
            return self.repair_options[0]
        return (f"{self.repair_options[0]}  or "
                f"{plural(others, 'other')}")

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
        first_position = {}
        for index, node in enumerate(tip_nodes):
            first_position.setdefault(node.label, index)

        # Running count of tips that matched the taxonomy, so "how many
        # resolved tips does this clade hold?" is two array lookups. Built
        # once for the tree rather than once per taxon.
        resolved_prefix = [0] * (len(tip_nodes) + 1)
        for index, node in enumerate(tip_nodes):
            resolved_prefix[index + 1] = (
                resolved_prefix[index] + (1 if tip_lineages.get(node.label)
                                          else 0))

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
            if result.conflict_cost == 1:
                result.repair_options = self._single_tip_repairs(
                    result, member_nodes, position, spans, resolved_prefix,
                    first_position)
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
    def _lowest_covering(positions, start_node, spans):
        """Lowest ancestor of start_node whose tips cover every position."""
        lo, hi = positions[0], positions[-1]
        node = start_node
        while node.parent is not None:
            first, last, _span = spans[id(node)]
            if first <= lo and last >= hi:
                break
            node = node.parent
        return node

    @staticmethod
    def _single_tip_repairs(result, member_nodes, position, spans,
                            resolved_prefix, first_position):
        """Every tip whose removal on its own would make the group a clade.

        A group can have several equally cheap repairs. In ((a1,b1),a2) both
        "move b1 out" and "move a2 across" cost one move, and reporting only
        whichever the tie-break happened to pick presents a guess as the
        answer. On a real 232-tip tree a third of the broken groups have more
        than one valid single-tip fix.

        Only meaningful at cost 1: beyond that a repair is a *set* of tips and
        enumerating the alternatives explodes combinatorially.

        Each candidate is tested directly - would deleting it leave the rest
        of the members forming a clade? - rather than by collecting ties
        during the main bottom-up pass. Tracking ties there cost 5x on the
        whole analysis, because the hot loop runs once per (node, taxon) pair
        and there are tens of thousands of those; this runs only for the
        handful of groups that cost exactly one move.
        """
        member_positions = sorted(position[id(n)] for n in member_nodes)
        if len(member_positions) < 2:
            return []
        anchor_by_position = {position[id(n)]: n for n in member_nodes}

        def repairs(label):
            dropped = first_position.get(label)
            if dropped is None:
                return False
            remaining = [p for p in member_positions if p != dropped]
            if not remaining:
                return False
            node = MonophylyChecker._lowest_covering(
                remaining, anchor_by_position[remaining[0]], spans)
            lo, hi, _span = spans[id(node)]
            inside = resolved_prefix[hi + 1] - resolved_prefix[lo]
            if lo <= dropped <= hi:
                inside -= 1          # the candidate itself is gone
            return inside == len(remaining)

        member_labels = set(result.members)
        foreign, own, seen = [], [], set()
        for label in list(result.intruders) + list(result.members):
            if label in seen:
                continue
            seen.add(label)
            if not repairs(label):
                continue
            (own if label in member_labels else foreign).append(label)
            if len(foreign) + len(own) >= MAX_REPAIR_OPTIONS:
                break

        # A foreign tip sitting inside the group is the actionable finding, so
        # it wins outright. Only when no intruder alone repairs the group -
        # typically two sampled members that simply are not sister - do we
        # fall back to naming the group's own tips, and then every one of them
        # is an equally good suspect and all get listed.
        return sorted(foreign) if foreign else sorted(own)

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

    def simulate_removal(self, root, tip_lineages, unresolved, labels,
                         baseline, min_members=2, progress_callback=None):
        """Re-analyse the tree with the named tips deleted.

        A genuine second analysis, not a prediction: the answer to "this tip
        repairs 10 groups" is confirmed or contradicted by actually removing
        it and counting again.

        The loaded tree is never touched - everything happens on a copy. Tip
        resolution is reused rather than redone, which keeps the experiment
        clean: topology is the only thing that differs between before and
        after, so nothing can be blamed on a name matching differently.
        """
        simulation = RemovalSimulation()
        simulation.tips_before = root.count_tips()
        simulation.broken_before = sum(1 for r in baseline
                                       if not r.is_monophyletic)

        if progress_callback:
            progress_callback(
                f"Simulating removal of {plural(len(labels), 'tip')}...")

        pruned, removed, missing = prune_tips(copy_tree(root), labels)
        simulation.removed = removed
        simulation.missing = missing

        if not removed:
            simulation.refused = "None of those labels are tips in this tree."
            return simulation
        if pruned is None:
            simulation.refused = (
                f"Removing {plural(len(removed), 'tip')} would leave fewer "
                f"than two behind. Nothing left to analyse.")
            return simulation

        simulation.tree = pruned
        simulation.tips_after = pruned.count_tips()

        # Unplaced tips are wildcards in the analysis, so deleting only those
        # provably changes nothing. Say so rather than let a null result read
        # as a malfunction.
        unresolved_set = set(unresolved)
        simulation.only_unresolved = all(label in unresolved_set
                                         for label in removed)

        after = self.check_taxa(pruned, tip_lineages, unresolved, min_members,
                               progress_callback)
        after_by_id = {r.tax_id: r for r in after}
        simulation.broken_after = sum(1 for r in after
                                      if not r.is_monophyletic)

        removed_counts = {}
        for label in removed:
            removed_counts[label] = removed_counts.get(label, 0) + 1
        resolved_removed = sum(count for label, count in removed_counts.items()
                               if label not in unresolved_set)
        total_resolved_after = (
            sum(1 for tip in root.get_tip_nodes()
                if tip_lineages.get(tip.label)) - resolved_removed)

        for result in baseline:
            survivor = after_by_id.get(result.tax_id)
            members_after = result.member_count - sum(
                removed_counts.get(label, 0) for label in set(result.members))

            if survivor is not None:
                cost_after = survivor.conflict_cost
                if not result.is_monophyletic:
                    if cost_after == 0:
                        outcome = "resolved"
                    elif cost_after < result.conflict_cost:
                        outcome = "improved"
                    else:
                        outcome = "unchanged"
                else:
                    outcome = "intact"
            else:
                cost_after = None
                # Worked out arithmetically rather than by asking check_taxa
                # why it skipped the group - the two skip rules are the same
                # ones it applies, and this way the numbers also fill the
                # "members before / after" column.
                if members_after < min_members:
                    outcome = "dropped_small"
                elif members_after >= total_resolved_after:
                    outcome = "dropped_covers_all"
                else:
                    outcome = "dropped_small"
            simulation.changes.append(
                GroupChange(result, outcome, members_after, cost_after))

        return simulation

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
        """Label every internal node with its taxonomic MRCA and status.

        One bottom-up pass. A node's MRCA is the lowest common ancestor, in
        the taxonomy, of its children's MRCAs; and the clade is exclusive
        exactly when the number of tips carrying that taxon anywhere in the
        tree equals the number inside the clade.

        The previous version rebuilt the clade's tip list, the complement
        list, and a set per tip lineage at every node, which is O(nodes x
        tips x lineage depth) - about 105 seconds on a 5000-tip tree. Counts
        replace all three, because nothing downstream reads the node-level
        tip or intruder lists: they feed the tree colouring and the exported
        labels, both of which need only status and name.
        """
        lineage_cache = {}

        def lineage_set(code):
            cached = lineage_cache.get(code)
            if cached is None:
                cached = lineage_cache[code] = set(
                    self.db.get_lineage_codes(code))
            return cached

        def common_ancestor(left, right):
            if left is None or right is None or left == right:
                return left if right is None else (right if left is None
                                                   else left)
            above = lineage_set(right)
            for code in self.db.get_lineage_codes(left):
                if code in above:
                    return code
            return None

        # How many tips in the whole tree carry each taxon.
        tip_nodes = root.get_tip_nodes()
        totals = {}
        for tip in tip_nodes:
            for code in tip_lineages.get(tip.label) or ():
                totals[code] = totals.get(code, 0) + 1

        mrca_of, resolved_of = {}, {}
        results = []
        for node in root.post_order():
            key = id(node)
            if node.is_tip:
                lineage = tip_lineages.get(node.label) or ()
                mrca_of[key] = lineage[0] if lineage else None
                resolved_of[key] = 1 if lineage else 0
                continue

            mrca, resolved = None, 0
            for child in node.children:
                child_mrca = mrca_of.pop(id(child))
                child_resolved = resolved_of.pop(id(child))
                resolved += child_resolved
                if child_resolved:      # all-unresolved subtrees say nothing
                    mrca = (child_mrca if mrca is None
                            else common_ancestor(mrca, child_mrca))
            mrca_of[key] = mrca
            resolved_of[key] = resolved

            if node.parent is None:
                continue                # the root is not a bipartition

            result = MonophylyResult(node, resolved)
            if mrca is None:
                result.status = "unresolved"
            else:
                result.mrca_code = mrca
                result.mrca_name = self.db.get_scientific_name(mrca)
                result.mrca_rank = self.db.code_to_rank.get(mrca, "")
                result.intruder_count = max(0, totals.get(mrca, 0) - resolved)
                result.status = ("not_monophyletic" if result.intruder_count
                                 else "monophyletic")
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
        show = self.db.display_label

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
                        f"  {show(tip.label)}: repairs "
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
                        lines.append(f"      blame: {show(g.culprit)}")
                    if g.intruders:
                        listed = ", ".join(map(show, g.intruders[:5]))
                        if len(g.intruders) > 5:
                            listed += f"... (+{len(g.intruders) - 5} more)"
                        lines.append(
                            f"      {plural(len(g.intruders), 'tip')} "
                            f"inside: {listed}")
                    if g.outliers:
                        listed = ", ".join(map(show, g.outliers[:5]))
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

            retired = [m for m in matches.values() if m.method == "merged"]
            if retired:
                lines.append("--- Outdated Taxon IDs ---")
                lines.append("NCBI merged these into newer IDs; the current "
                             "taxon was used.")
                for m in sorted(retired, key=lambda x: x.label):
                    lines.append(f"  {m.label} -> {m.tax_id} ({m.matched_name})")
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


class GroupChange:
    """What removing some tips did to one taxonomic group.

    `outcome` is deliberately more than a boolean, because "no longer in the
    report" is not the same as "fixed" and conflating them would flatter the
    result:

      resolved            was broken, is a clade in the pruned tree
      improved            costs fewer moves, still broken
      unchanged           same cost, still broken
      dropped_small       too few sampled tips remain to test it at all
      dropped_covers_all  now spans every remaining resolved tip, so it says
                          nothing - skipped exactly as a first analysis would
      intact              was already monophyletic
    """

    __slots__ = ("tax_id", "name", "rank", "outcome", "members_before",
                 "members_after", "cost_before", "cost_after",
                 "culprit_before")

    def __init__(self, result, outcome, members_after, cost_after):
        self.tax_id = result.tax_id
        self.name = result.name
        self.rank = result.rank
        self.outcome = outcome
        self.members_before = result.member_count
        self.members_after = members_after
        self.cost_before = result.conflict_cost
        self.cost_after = cost_after       # None when no longer reported
        self.culprit_before = result.culprit

    @property
    def was_broken(self):
        return self.cost_before > 0

    def __repr__(self):
        return (f"GroupChange({self.name}, {self.outcome}, "
                f"{self.cost_before}->{self.cost_after})")


class RemovalSimulation:
    """The result of asking "what if these tips were not in the tree?"."""

    __slots__ = ("removed", "missing", "tree", "changes", "tips_before",
                 "tips_after", "broken_before", "broken_after", "refused",
                 "only_unresolved")

    def __init__(self):
        self.removed = []
        self.missing = []
        self.tree = None
        self.changes = []
        self.tips_before = 0
        self.tips_after = 0
        self.broken_before = 0
        self.broken_after = 0
        self.refused = ""          # why the simulation could not run
        self.only_unresolved = False

    def of(self, outcome):
        return [c for c in self.changes if c.outcome == outcome]

    @property
    def resolved(self):
        return self.of("resolved")

    @property
    def newly_broken(self):
        """Groups that were intact and are now broken.

        Removing tips cannot do this - every clade's cost is non-increasing
        under deletion - so anything here is a bug, not a finding.
        """
        return [c for c in self.changes
                if not c.was_broken and (c.cost_after or 0) > 0]

    def counts(self):
        return {name: len(self.of(name)) for name in
                ("resolved", "improved", "unchanged", "dropped_small",
                 "dropped_covers_all", "intact")}

    def __repr__(self):
        return (f"RemovalSimulation({len(self.removed)} removed, "
                f"{self.broken_before}->{self.broken_after} broken)")


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
