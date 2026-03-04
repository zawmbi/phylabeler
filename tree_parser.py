"""
tree_parser.py - Robust Newick Tree Parser

Improved Newick parser that handles:
- Quoted labels, branch lengths, bootstrap values
- Polytomies (multifurcating nodes)
- Scientific notation in branch lengths
- Post-order traversal for bipartitions
- Tree-to-text representation

Pure Python stdlib - no external dependencies.
"""

import re


def clean_newick_name(name):
    """Remove characters that are illegal in Newick format from a name."""
    for ch in (" ", ".", ",", "'", ")", "(", ";", "[", "]"):
        name = name.replace(ch, "")
    return name


class Node:
    """A node in a phylogenetic tree.

    Supports hierarchical tree structures with labels, branch lengths,
    parent-child relationships, and bipartition computation.
    """

    __slots__ = ("label", "branch_length", "parent", "children", "is_tip",
                 "monophyly_status", "taxonomy_code", "taxonomy_name")

    def __init__(self, label="", branch_length=None, is_tip=False):
        self.label = label
        self.branch_length = branch_length
        self.parent = None
        self.children = []
        self.is_tip = is_tip
        self.monophyly_status = None   # "monophyletic", "not_monophyletic", None
        self.taxonomy_code = None
        self.taxonomy_name = None

    def add_child(self, child):
        self.children.append(child)
        child.parent = self

    def remove_child(self, child):
        self.children.remove(child)
        child.parent = None

    # ---- Tip collection ----

    def get_tip_labels(self):
        """Recursively collect all tip labels in this clade."""
        tips = []
        self._collect_tips(tips)
        return tips

    def _collect_tips(self, tips):
        if self.is_tip:
            tips.append(self.label)
        for child in self.children:
            child._collect_tips(tips)

    def get_tip_nodes(self):
        """Recursively collect all tip Node objects."""
        nodes = []
        self._collect_tip_nodes(nodes)
        return nodes

    def _collect_tip_nodes(self, nodes):
        if self.is_tip:
            nodes.append(self)
        for child in self.children:
            child._collect_tip_nodes(nodes)

    # ---- Bipartition computation (post-order) ----

    def get_bipartitions(self):
        """Get bipartitions for all internal nodes using post-order traversal.

        Returns list of (tip_label_set, node_reference) tuples.
        Only internal (non-root) nodes are included.
        """
        bips = []
        self._post_order_bips(bips)
        return bips

    def _post_order_bips(self, bips):
        """Post-order traversal to collect bipartitions."""
        # Process children first (post-order)
        for child in self.children:
            child._post_order_bips(bips)

        # Then process this node
        if not self.is_tip and self.parent is not None:
            tips = self.get_tip_labels()
            bips.append((tips, self))

    # ---- Newick output ----

    def to_newick(self, branch_lengths=True):
        """Generate Newick string representation of the subtree."""
        ret = ""
        if self.children:
            child_strings = []
            for child in self.children:
                child_strings.append(child.to_newick(branch_lengths))
            ret = "(" + ",".join(child_strings) + ")"

        if self.label:
            ret += str(self.label)

        if branch_lengths and self.branch_length is not None:
            ret += ":" + str(self.branch_length)

        return ret

    # ---- Tree statistics ----

    def count_tips(self):
        if self.is_tip:
            return 1
        return sum(child.count_tips() for child in self.children)

    def count_internal(self):
        if self.is_tip:
            return 0
        return 1 + sum(child.count_internal() for child in self.children)

    def depth(self):
        if self.is_tip:
            return 0
        return 1 + max((child.depth() for child in self.children), default=0)

    # ---- Tree traversal ----

    def post_order(self):
        """Iterate over all nodes in post-order."""
        for child in self.children:
            yield from child.post_order()
        yield self

    def pre_order(self):
        """Iterate over all nodes in pre-order."""
        yield self
        for child in self.children:
            yield from child.pre_order()

    # ---- Label assignment ----

    def assign_labels(self, bipart_and_name, auto_select=True):
        """Assign taxonomy labels to internal nodes based on bipartition matching.

        bipart_and_name: list of (bipartition_tips, [candidate_names]) tuples
        auto_select: if True, automatically pick first name; if False, skip ambiguous
        """
        if not self.is_tip:
            my_tips = sorted(self.get_tip_labels())
            for tips, names in bipart_and_name:
                if sorted(tips) == my_tips:
                    if names and names[0] != "":
                        selected = names[0] if auto_select else (
                            names[1] if len(names) > 1 else names[0]
                        )
                        self.label = clean_newick_name(selected)
                    break

        for child in self.children:
            child.assign_labels(bipart_and_name, auto_select)

    def __repr__(self):
        if self.is_tip:
            return f"Node(tip={self.label})"
        n = self.count_tips()
        return f"Node(label={self.label}, tips={n})"


# ---- Newick Parsing ----

class NewickParser:
    """Improved Newick format parser.

    Handles quoted labels, branch lengths, bootstrap values,
    polytomies, and various edge cases in real-world Newick strings.
    """

    def __init__(self, string):
        self.s = string.strip().rstrip(";").strip()
        self.pos = 0

    def parse(self):
        """Parse the Newick string and return the root Node."""
        root = Node()
        self._parse_node(root)
        return root

    def _parse_node(self, node):
        """Recursively parse a node and its children."""
        if self.pos >= len(self.s):
            return

        if self.s[self.pos] == "(":
            self.pos += 1  # skip (
            self._parse_children(node)

            # After closing ), read node's own label/length
            if self.pos < len(self.s):
                label, length = self._read_label_and_length()
                node.label = label
                node.branch_length = length

        else:
            # This is a tip node
            label, length = self._read_label_and_length()
            node.label = label
            node.branch_length = length
            node.is_tip = True

    def _parse_children(self, parent):
        """Parse comma-separated children until closing paren."""
        while True:
            child = Node()

            if self.pos < len(self.s) and self.s[self.pos] == "(":
                # Internal child
                self.pos += 1
                self._parse_children(child)
                label, length = self._read_label_and_length()
                child.label = label
                child.branch_length = length
            else:
                # Tip child
                label, length = self._read_label_and_length()
                child.label = label
                child.branch_length = length
                child.is_tip = True

            parent.add_child(child)

            if self.pos < len(self.s) and self.s[self.pos] == ",":
                self.pos += 1  # skip comma
                continue
            elif self.pos < len(self.s) and self.s[self.pos] == ")":
                self.pos += 1  # skip )
                break
            else:
                break

    def _read_label_and_length(self):
        """Read a label and optional branch length at current position."""
        label = ""
        length = None

        if self.pos >= len(self.s):
            return label, length

        # Read quoted label
        if self.s[self.pos] == "'":
            label = self._read_quoted("'")
        elif self.s[self.pos] == '"':
            label = self._read_quoted('"')
        else:
            # Read unquoted label
            while self.pos < len(self.s) and self.s[self.pos] not in (
                ":", ",", "(", ")", ";", " ", "\t", "\n"
            ):
                label += self.s[self.pos]
                self.pos += 1

        # Read branch length
        if self.pos < len(self.s) and self.s[self.pos] == ":":
            self.pos += 1  # skip :
            length_str = ""
            while self.pos < len(self.s) and self.s[self.pos] in (
                "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                ".", "-", "e", "E", "+"
            ):
                length_str += self.s[self.pos]
                self.pos += 1
            if length_str:
                length = length_str

        return label, length

    def _read_quoted(self, quote_char):
        """Read a quoted string."""
        self.pos += 1  # skip opening quote
        result = ""
        while self.pos < len(self.s):
            if self.s[self.pos] == quote_char:
                self.pos += 1
                break
            result += self.s[self.pos]
            self.pos += 1
        return result


def parse_newick(newick_string):
    """Parse a Newick format string and return the root Node."""
    parser = NewickParser(newick_string)
    return parser.parse()


def parse_newick_file(filepath):
    """Parse a Newick file and return the root Node."""
    with open(filepath, "r") as f:
        content = f.read().strip()
    return parse_newick(content)
