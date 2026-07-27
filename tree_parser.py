"""
tree_parser.py - Robust Newick Tree Parser

Improved Newick parser that handles:
- Quoted labels, branch lengths, bootstrap/support values
- [...] comments, including NHX-style annotations
- Polytomies (multifurcating nodes)
- Scientific notation in branch lengths
- Post-order traversal for bipartitions
- Tree-to-text representation

Support values on internal nodes are kept separate from taxonomy labels, so
labelling a clade never destroys the bootstrap value the tree arrived with.

Parsing and traversal are iterative rather than recursive, so ladder-shaped
trees of any depth are safe - a recursive walk overflows Python's stack at a
few hundred tips.

Tree files may be plain or compressed (.nwk.gz, .tre.bz2, ...).

Pure Python stdlib - no external dependencies.
"""

from archive_io import open_text, strip_compression_suffix

# Extensions understood as Newick trees, before any compression suffix.
TREE_SUFFIXES = (".nwk", ".tre", ".tree", ".newick")

# Characters that terminate an unquoted Newick label.
_LABEL_STOP = set(":,();[] \t\n\r")
_LENGTH_CHARS = set("0123456789.-+eE")


def is_tree_filename(name):
    """True if a filename looks like a Newick tree, compressed or not."""
    return strip_compression_suffix(name).lower().endswith(TREE_SUFFIXES)


def clean_newick_name(name):
    """Make a name safe to write as an unquoted Newick label.

    Spaces become underscores rather than vanishing, so "Homo sapiens" stays
    readable as "Homo_sapiens" instead of collapsing to "Homosapiens".
    """
    name = name.strip().replace(" ", "_")
    for ch in (",", "'", '"', ")", "(", ";", "[", "]", ":"):
        name = name.replace(ch, "")
    return name


def _is_number(text):
    """True if text parses as a float (bootstrap values, posteriors)."""
    try:
        float(text)
    except (TypeError, ValueError):
        return False
    return True


def quote_if_needed(label):
    """Single-quote a label that Newick cannot represent bare.

    Without this a label like "Homo sapiens" would be written unquoted and
    the resulting file would not parse anywhere else.
    """
    if not label:
        return label
    if any(ch in label for ch in " \t\n()[]',:;"):
        return "'" + label.replace("'", "''") + "'"
    return label


class Node:
    """A node in a phylogenetic tree.

    Supports hierarchical tree structures with labels, branch lengths,
    parent-child relationships, and bipartition computation.
    """

    __slots__ = ("label", "support", "branch_length", "parent", "children",
                 "is_tip", "monophyly_status", "taxonomy_code",
                 "taxonomy_name")

    def __init__(self, label="", branch_length=None, is_tip=False):
        self.label = label
        self.support = None            # Bootstrap/posterior from the input
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

    # ---- Traversal (iterative: safe at any depth) ----

    def post_order(self):
        """Iterate over all nodes in post-order."""
        stack = [(self, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                yield node
            else:
                stack.append((node, True))
                for child in reversed(node.children):
                    stack.append((child, False))

    def pre_order(self):
        """Iterate over all nodes in pre-order."""
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            for child in reversed(node.children):
                stack.append(child)

    def ancestors(self):
        """Iterate from this node's parent up to the root."""
        node = self.parent
        while node is not None:
            yield node
            node = node.parent

    # ---- Tip collection ----

    def get_tip_labels(self):
        """Collect all tip labels in this clade, left to right."""
        return [n.label for n in self.get_tip_nodes()]

    def get_tip_nodes(self):
        """Collect all tip Node objects in this clade, left to right."""
        tips = []
        stack = [self]
        while stack:
            node = stack.pop()
            if node.is_tip:
                tips.append(node)
            else:
                stack.extend(reversed(node.children))
        return tips

    # ---- Bipartition computation (post-order) ----

    def get_bipartitions(self):
        """Get bipartitions for all internal nodes using post-order traversal.

        Returns list of (tip_label_list, node_reference) tuples for internal
        non-root nodes. Tip sets accumulate bottom-up in a single pass rather
        than being recollected at every node.
        """
        tips_below = {}
        bips = []
        for node in self.post_order():
            if node.is_tip:
                tips_below[id(node)] = [node.label]
                continue
            collected = []
            for child in node.children:
                collected.extend(tips_below.pop(id(child)))
            tips_below[id(node)] = collected
            if node.parent is not None:
                bips.append((collected, node))
        return bips

    # ---- Newick output ----

    def _newick_suffix(self, branch_lengths, label_format):
        """Render this node's own label/support/length, without children."""
        taxon = str(self.label) if self.label else ""
        support = str(self.support) if self.support is not None else ""

        if self.is_tip:
            out = taxon
        elif label_format == "support":
            out = support or taxon
        elif label_format == "taxon":
            out = taxon or support
        elif taxon and support:
            out = f"{taxon}/{support}"   # keep both, e.g. Anthocerotophyta/100
        else:
            out = taxon or support

        out = quote_if_needed(out)
        if branch_lengths and self.branch_length is not None:
            out += ":" + str(self.branch_length)
        return out

    def to_newick(self, branch_lengths=True, label_format="taxon_support"):
        """Generate a Newick string for the subtree.

        label_format controls what internal nodes carry:
          "taxon_support"  Anthocerotophyta/100  (default, keeps both)
          "taxon"          Anthocerotophyta
          "support"        100                   (original tree unchanged)
        """
        rendered = {}
        for node in self.post_order():
            if node.children:
                inner = ",".join(rendered.pop(id(c)) for c in node.children)
                text = "(" + inner + ")"
            else:
                text = ""
            rendered[id(node)] = text + node._newick_suffix(
                branch_lengths, label_format)
        return rendered[id(self)]

    # ---- Tree statistics ----

    def count_tips(self):
        count = 0
        stack = [self]
        while stack:
            node = stack.pop()
            if node.is_tip:
                count += 1
            else:
                stack.extend(node.children)
        return count

    def count_internal(self):
        count = 0
        stack = [self]
        while stack:
            node = stack.pop()
            if not node.is_tip:
                count += 1
                stack.extend(node.children)
        return count

    def depth(self):
        """Maximum number of edges from this node down to a tip."""
        best = 0
        stack = [(self, 0)]
        while stack:
            node, d = stack.pop()
            if node.is_tip:
                best = max(best, d)
            else:
                for child in node.children:
                    stack.append((child, d + 1))
        return best

    def __repr__(self):
        if self.is_tip:
            return f"Node(tip={self.label})"
        return f"Node(label={self.label}, tips={self.count_tips()})"


def index_tips(root):
    """Index a tree for fast clade queries.

    Returns (tip_nodes, spans) where tip_nodes is the left-to-right list of
    tips and spans maps id(node) -> (first_tip_index, last_tip_index, n_tips).
    A node's span covers a contiguous index range, so the MRCA of a set of
    tips is the lowest ancestor whose span covers their min..max index.
    """
    tip_nodes = root.get_tip_nodes()
    position = {id(tip): i for i, tip in enumerate(tip_nodes)}

    spans = {}
    for node in root.post_order():
        if node.is_tip:
            i = position[id(node)]
            spans[id(node)] = (i, i, 1)
        else:
            lo = hi = None
            total = 0
            for child in node.children:
                c_lo, c_hi, c_n = spans[id(child)]
                lo = c_lo if lo is None else min(lo, c_lo)
                hi = c_hi if hi is None else max(hi, c_hi)
                total += c_n
            spans[id(node)] = (lo or 0, hi or 0, total)
    return tip_nodes, spans


# ---- Newick Parsing ----

class NewickParser:
    """Iterative Newick format parser.

    Handles quoted labels, branch lengths, support values, [...] comments,
    polytomies, and arbitrary nesting depth.
    """

    def __init__(self, string):
        self.s = string.strip()
        self.pos = 0

    def parse(self):
        """Parse the Newick string and return the root Node."""
        root = Node()
        current = root
        s = self.s
        n = len(s)

        while self.pos < n:
            c = s[self.pos]

            if c == "(":
                child = Node()
                current.add_child(child)
                current = child
                self.pos += 1

            elif c == ",":
                # Close the current node and start a sibling under its parent.
                if current.parent is None:
                    # Malformed input: no open paren. Ignore the separator.
                    self.pos += 1
                    continue
                current = current.parent
                child = Node()
                current.add_child(child)
                current = child
                self.pos += 1

            elif c == ")":
                if current.parent is not None:
                    current = current.parent
                self.pos += 1
                self._read_label_and_length(current)

            elif c == "[":
                self._skip_comment()

            elif c == ";":
                break

            elif c in " \t\n\r":
                self.pos += 1

            else:
                self._read_label_and_length(current)

        self._finalize(root)
        return root

    @staticmethod
    def _finalize(root):
        """Mark tips and split numeric internal labels out as support."""
        for node in root.post_order():
            node.is_tip = not node.children
            if not node.is_tip and node.label and _is_number(node.label):
                node.support = node.label
                node.label = ""

    def _skip_comment(self):
        """Skip a [...] comment, including NHX annotations."""
        depth = 0
        while self.pos < len(self.s):
            c = self.s[self.pos]
            if c == "[":
                depth += 1
            elif c == "]":
                depth -= 1
                if depth == 0:
                    self.pos += 1
                    return
            self.pos += 1

    def _read_label_and_length(self, node):
        """Read a label and optional branch length into node."""
        s = self.s
        n = len(s)

        if self.pos < n and s[self.pos] in ("'", '"'):
            label = self._read_quoted(s[self.pos])
        else:
            start = self.pos
            while self.pos < n and s[self.pos] not in _LABEL_STOP:
                self.pos += 1
            label = s[start:self.pos]

        if label:
            node.label = label

        if self.pos < n and s[self.pos] == "[":
            self._skip_comment()

        if self.pos < n and s[self.pos] == ":":
            self.pos += 1
            if self.pos < n and s[self.pos] == "[":
                self._skip_comment()
            start = self.pos
            while self.pos < n and s[self.pos] in _LENGTH_CHARS:
                self.pos += 1
            length = s[start:self.pos]
            if length:
                node.branch_length = length
            if self.pos < n and s[self.pos] == "[":
                self._skip_comment()

    def _read_quoted(self, quote_char):
        """Read a quoted string, handling doubled quotes as escapes."""
        self.pos += 1  # skip opening quote
        chars = []
        while self.pos < len(self.s):
            c = self.s[self.pos]
            if c == quote_char:
                # A doubled quote is a literal quote character.
                if self.pos + 1 < len(self.s) and self.s[self.pos + 1] == quote_char:
                    chars.append(quote_char)
                    self.pos += 2
                    continue
                self.pos += 1
                break
            chars.append(c)
            self.pos += 1
        return "".join(chars)


def parse_newick(newick_string):
    """Parse a Newick format string and return the root Node."""
    return NewickParser(newick_string).parse()


def parse_newick_file(filepath):
    """Parse a Newick file and return the root Node.

    Transparently reads gzip/bz2/xz-compressed trees. Files holding several
    trees (one per line) yield the first one.
    """
    with open_text(filepath) as f:
        content = f.read().strip()
    return parse_newick(content)
