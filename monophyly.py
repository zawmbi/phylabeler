"""
monophyly.py - Monophyly Checking Engine

Checks whether clades in a gene tree are monophyletic with respect to
the NCBI taxonomy. Inspired by MonoPhy (Schwery & O'Meara, 2016).

Core algorithm:
1. For each internal node's bipartition, find the most recent common
   ancestor (MRCA) in the NCBI taxonomy.
2. Verify that no taxa outside the bipartition share that MRCA
   (i.e., the MRCA is exclusive to the clade).
3. Label the node with the MRCA's taxonomy name if monophyletic.

Pure Python stdlib - no external dependencies.
"""


class MonophylyResult:
    """Result of monophyly analysis for a single node."""

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


class MonophylyChecker:
    """Check monophyly of tree clades against NCBI taxonomy.

    Uses post-order traversal and bipartition analysis (like MonoPhy).
    """

    def __init__(self, taxonomy_db):
        self.db = taxonomy_db

    def check_tree(self, root, progress_callback=None):
        """Run monophyly analysis on an entire tree.

        Args:
            root: Root Node of the parsed tree.
            progress_callback: Optional function(str) for status updates.

        Returns:
            list of MonophylyResult objects for each internal node.
        """
        all_tips = root.get_tip_labels()

        # Resolve tip names to taxonomy codes
        tip_lineages = {}
        unresolved = []
        for tip in all_tips:
            code = self.db.lookup_name(tip)
            if code:
                tip_lineages[tip] = self.db.get_lineage_codes(code)
            else:
                unresolved.append(tip)
                tip_lineages[tip] = []

        if progress_callback:
            resolved = len(all_tips) - len(unresolved)
            progress_callback(
                f"Resolved {resolved}/{len(all_tips)} tip names "
                f"({len(unresolved)} unresolved)"
            )

        # Get bipartitions via post-order traversal
        bipartitions = root.get_bipartitions()

        results = []
        for idx, (bip_tips, node) in enumerate(bipartitions):
            if progress_callback and idx % 10 == 0:
                progress_callback(
                    f"Checking node {idx + 1}/{len(bipartitions)}..."
                )

            result = self._check_bipartition(
                bip_tips, all_tips, tip_lineages, node
            )
            results.append(result)

            # Annotate the tree node
            node.monophyly_status = result.status
            node.taxonomy_code = result.mrca_code
            node.taxonomy_name = result.mrca_name

        if progress_callback:
            mono = sum(1 for r in results if r.status == "monophyletic")
            non_mono = sum(1 for r in results if r.status == "not_monophyletic")
            progress_callback(
                f"Done: {mono} monophyletic, {non_mono} non-monophyletic, "
                f"{len(results) - mono - non_mono} unresolved"
            )

        return results, unresolved

    def _check_bipartition(self, bip_tips, all_tips, tip_lineages, node):
        """Check if a single bipartition is monophyletic."""
        result = MonophylyResult(node, bip_tips)

        # Get the complement (taxa NOT in this clade)
        complement = [t for t in all_tips if t not in set(bip_tips)]

        # Collect lineages for clade members (skip unresolved)
        clade_lineages = [
            tip_lineages[t] for t in bip_tips if tip_lineages[t]
        ]

        if not clade_lineages:
            result.status = "unresolved"
            return result

        # Find MRCA: the first (most recent) code shared by ALL clade members
        mrca = self._find_mrca(clade_lineages)

        if not mrca:
            result.status = "unresolved"
            return result

        result.mrca_code = mrca
        result.mrca_name = self.db.get_scientific_name(mrca)
        result.mrca_rank = self.db.code_to_rank.get(mrca, "")
        result.candidate_names = self.db.get_names_for_code(mrca)

        # Check exclusivity: MRCA should NOT appear in complement lineages
        for tip in complement:
            if tip_lineages[tip] and mrca in tip_lineages[tip]:
                result.intruders.append(tip)

        if result.intruders:
            result.status = "not_monophyletic"
        else:
            result.status = "monophyletic"

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

        # Convert other lineages to sets for fast lookup
        other_sets = [set(lin) for lin in lineage_list[1:]]

        for code in lineage_list[0]:
            if all(code in s for s in other_sets):
                return code

        return None

    def label_tree(self, root, results, auto_select=True):
        """Apply labels from monophyly results to the tree.

        Labels monophyletic nodes with their MRCA taxonomy name.
        Non-monophyletic nodes can optionally be marked.
        """
        for result in results:
            if result.status == "monophyletic" and result.mrca_name:
                name = clean_label(result.mrca_name)
                result.node.label = name

    def get_summary(self, results):
        """Generate a summary report of monophyly analysis."""
        mono = [r for r in results if r.status == "monophyletic"]
        non_mono = [r for r in results if r.status == "not_monophyletic"]
        unresolved = [r for r in results if r.status == "unresolved"]

        lines = []
        lines.append(f"=== Monophyly Analysis Summary ===")
        lines.append(f"Total internal nodes: {len(results)}")
        lines.append(f"Monophyletic: {len(mono)}")
        lines.append(f"Non-monophyletic: {len(non_mono)}")
        lines.append(f"Unresolved: {len(unresolved)}")
        lines.append("")

        if mono:
            lines.append("--- Monophyletic Clades ---")
            for r in sorted(mono, key=lambda x: len(x.tip_labels), reverse=True):
                lines.append(
                    f"  {r.mrca_name} ({r.mrca_rank}) "
                    f"[{len(r.tip_labels)} tips]"
                )
            lines.append("")

        if non_mono:
            lines.append("--- Non-Monophyletic Clades ---")
            for r in sorted(non_mono, key=lambda x: len(x.tip_labels), reverse=True):
                intruder_str = ", ".join(r.intruders[:5])
                if len(r.intruders) > 5:
                    intruder_str += f"... (+{len(r.intruders) - 5} more)"
                lines.append(
                    f"  {r.mrca_name} ({r.mrca_rank}) "
                    f"[{len(r.tip_labels)} tips] - "
                    f"intruders: {intruder_str}"
                )

        return "\n".join(lines)


def clean_label(name):
    """Clean a taxonomy name for use in Newick labels."""
    for ch in (" ", ".", ",", "'", ")", "(", ";"):
        name = name.replace(ch, "")
    return name
