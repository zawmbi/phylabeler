"""
taxonomy_db.py - NCBI Taxonomy Database with Caching

Downloads, parses, and caches NCBI taxonomy data (names.dmp, nodes.dmp).
Uses pickle for caching so the GenBank tree doesn't need to be regenerated
every time. Pure Python stdlib - no external dependencies.

Source files may be supplied in any of the forms NCBI ships or that users
end up with locally: extracted names.dmp/nodes.dmp, gzipped copies of them,
or the un-extracted taxdump.tar.gz / taxdump.tgz / taxdmp.zip archive.
"""

import os
import pickle
import urllib.request
import time

from archive_io import ArchiveError, is_archive, open_text, stream_members

DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".phylabeler_cache")
NCBI_TAXONOMY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# The only two members of an NCBI taxonomy dump this tool needs.
TAXDUMP_MEMBERS = ("names.dmp", "nodes.dmp")

# Tokens gene-tree tips commonly carry alongside the organism name.
_TIP_SEPARATORS = ("|", "/", "@", "#")


def normalize_name(name):
    """Fold a name to the form used as a lookup key.

    Tip labels arrive as "Homo_sapiens", "homo sapiens", "Homo  sapiens"; all
    three must find tax_id 9606. Normalising at both insert and lookup time
    makes matching O(1) instead of a scan over millions of names.
    """
    return " ".join(name.replace("_", " ").split()).lower()


class TipMatch:
    """How a tree tip label was resolved against the taxonomy.

    `method` records the strength of the match so a genus-level fallback is
    never mistaken for a confirmed species identification:
      "exact"    the whole label matched a taxon name
      "trimmed"  matched after stripping a gene/accession suffix
      "partial"  matched after dropping trailing name tokens
      "genus"    only the first token matched
    """

    __slots__ = ("label", "tax_id", "matched_name", "method", "rank",
                 "ambiguous")

    def __init__(self, label, tax_id, matched_name, method, rank,
                 ambiguous=False):
        self.label = label
        self.tax_id = tax_id
        self.matched_name = matched_name
        self.method = method
        self.rank = rank
        self.ambiguous = ambiguous

    @property
    def is_exact(self):
        return self.method == "exact"

    def __repr__(self):
        return (f"TipMatch({self.label} -> {self.tax_id} "
                f"via {self.method})")


class TaxonomyDB:
    """Cached NCBI taxonomy database.

    Parses names.dmp and nodes.dmp into efficient lookup structures,
    then caches them as pickle files so subsequent runs are fast.
    """

    def __init__(self, cache_dir=None):
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)

        # Core data structures
        self.name_to_code = {}       # normalized name -> tax_id
        self.code_to_names = {}      # tax_id -> [list of names]
        self.code_to_scientific = {}  # tax_id -> scientific name
        self.code_to_parent = {}     # tax_id -> parent tax_id
        self.code_to_rank = {}       # tax_id -> rank string
        self.homonyms = {}           # normalized name -> [competing tax_ids]
        self.loaded = False
        self._cache_version = 3      # Bump to invalidate old caches

    # ---- Public API ----

    def load(self, names_file=None, nodes_file=None, archive=None,
             progress_callback=None):
        """Load taxonomy data from cache, source files, or an NCBI archive.

        With no arguments, loads a previously built cache. `names_file` and
        `nodes_file` are the extracted dmp files and may be gzipped.
        `archive` is an un-extracted dump (taxdump.tar.gz, .tgz, taxdmp.zip);
        passing one as `names_file` with no `nodes_file` works too, so a
        single file picker can accept either form.
        """
        if archive is None and names_file and not nodes_file:
            if is_archive(names_file):
                archive, names_file = names_file, None

        if archive:
            return self.load_archive(archive, progress_callback)

        if names_file is None and nodes_file is None:
            if self._load_cache(progress_callback):
                return True
            # The cache may merely be stale after a format change. Rebuild it
            # from whatever the user already downloaded rather than making
            # them fetch 60 MB from NCBI again.
            return self._rebuild_from_cache_dir(progress_callback)

        if names_file and nodes_file:
            self._parse_names(names_file, progress_callback)
            self._parse_nodes(nodes_file, progress_callback)
            self._save_cache(progress_callback)
            self.loaded = True
            return True

        return False

    def _rebuild_from_cache_dir(self, progress_callback=None):
        """Re-parse taxonomy from source files left in the cache directory."""
        names_path = os.path.join(self.cache_dir, "names.dmp")
        nodes_path = os.path.join(self.cache_dir, "nodes.dmp")
        if os.path.exists(names_path) and os.path.exists(nodes_path):
            if progress_callback:
                progress_callback("Rebuilding cache from names.dmp/nodes.dmp...")
            self._parse_names(names_path, progress_callback)
            self._parse_nodes(nodes_path, progress_callback)
            self._save_cache(progress_callback)
            self.loaded = True
            return True

        for name in ("taxdump.tar.gz", "taxdump.tgz", "taxdmp.zip"):
            archive = os.path.join(self.cache_dir, name)
            if os.path.exists(archive):
                if progress_callback:
                    progress_callback(f"Rebuilding cache from {name}...")
                return self.load_archive(archive, progress_callback)

        return False

    def load_archive(self, archive_path, progress_callback=None):
        """Parse names.dmp and nodes.dmp straight out of an NCBI dump archive.

        Accepts taxdump.tar.gz, taxdump.tgz, taxdump.tar and taxdmp.zip. The
        archive is read in a single streaming pass and nothing is unpacked to
        disk - only the parsed pickle cache is written.
        """
        if progress_callback:
            progress_callback(
                f"Reading {os.path.basename(archive_path)}...")

        seen = set()
        for basename, stream in stream_members(archive_path,
                                               wanted=TAXDUMP_MEMBERS):
            if basename in seen:
                continue
            seen.add(basename)
            if basename == "names.dmp":
                self._parse_names_stream(stream, progress_callback)
            else:
                self._parse_nodes_stream(stream, progress_callback)

        missing = [m for m in TAXDUMP_MEMBERS if m not in seen]
        if missing:
            raise ArchiveError(
                f"{archive_path} does not contain "
                f"{' and '.join(missing)}. Expected an NCBI taxonomy dump "
                f"such as taxdump.tar.gz or taxdmp.zip."
            )

        self._save_cache(progress_callback)
        self.loaded = True
        return True

    def download_and_load(self, progress_callback=None):
        """Download NCBI taxonomy dump and load it.

        The archive is parsed in place, so the ~450 MB of dmp files it
        contains are never written to the cache directory.
        """
        tar_path = os.path.join(self.cache_dir, "taxdump.tar.gz")
        self._download_taxonomy(tar_path, progress_callback)
        return self.load_archive(tar_path, progress_callback)

    def get_lineage(self, tax_id):
        """Get full lineage from a tax_id up to root.

        Returns: list of (tax_id, rank, name) tuples from tip to root.
        """
        if not self.loaded:
            return []

        lineage = []
        current = str(tax_id)
        visited = set()

        while current != "1" and current in self.code_to_parent:
            if current in visited:
                break
            visited.add(current)

            rank = self.code_to_rank.get(current, "")
            name = self.code_to_scientific.get(current, "")
            lineage.append((current, rank, name))
            current = self.code_to_parent[current]

        # Add root
        if current == "1":
            lineage.append(("1", "no rank", "root"))

        return lineage

    def get_lineage_codes(self, tax_id):
        """Get just the tax_id codes in the lineage."""
        return [x[0] for x in self.get_lineage(tax_id)]

    def lookup_name(self, name):
        """Look up a taxon name and return its tax_id, or None."""
        match = self.resolve_tip(name)
        return match.tax_id if match else None

    def resolve_tip(self, label):
        """Resolve a tree tip label to a taxon, reporting how it matched.

        Gene-tree tips are rarely bare organism names: they carry accessions
        ("Anthoceros_agrestis_Bonn"), gene ids ("Homo_sapiens|COX1") and
        strain suffixes. Rather than giving up, progressively drop trailing
        tokens - but say so, because a genus-level match is a weaker claim
        than a species-level one.

        Returns a TipMatch, or None if nothing matched.
        """
        if not label:
            return None

        # Whole label, normalized.
        match = self._match(label, "exact")
        if match:
            return match

        # Strip anything after a gene/accession separator.
        head = label
        for sep in _TIP_SEPARATORS:
            if sep in head:
                head = head.split(sep, 1)[0]
        if head != label:
            match = self._match(head, "trimmed")
            if match:
                return match

        # Drop trailing tokens: "Anthoceros agrestis Bonn" -> "Anthoceros
        # agrestis" -> "Anthoceros". Stop at the genus.
        tokens = normalize_name(head).split()
        for end in range(len(tokens) - 1, 0, -1):
            candidate = " ".join(tokens[:end])
            code = self.name_to_code.get(candidate)
            if code:
                rank = self.code_to_rank.get(code, "")
                method = "genus" if end == 1 else "partial"
                return TipMatch(label, code, candidate, method, rank,
                                candidate in self.homonyms)
        return None

    def _match(self, text, method):
        key = normalize_name(text)
        code = self.name_to_code.get(key)
        if not code:
            return None
        return TipMatch(label=text, tax_id=code, matched_name=key,
                        method=method, rank=self.code_to_rank.get(code, ""),
                        ambiguous=key in self.homonyms)

    def get_homonyms(self, name):
        """Return competing tax_ids for an ambiguous name, newest listing first."""
        return self.homonyms.get(normalize_name(name), [])

    def get_names_for_code(self, tax_id):
        """Get all names associated with a tax_id."""
        return self.code_to_names.get(str(tax_id), [])

    def get_scientific_name(self, tax_id):
        """Get the scientific name for a tax_id."""
        return self.code_to_scientific.get(str(tax_id), "")

    def get_cache_info(self):
        """Return info about the cached taxonomy data."""
        cache_file = os.path.join(self.cache_dir, "taxonomy.pkl")
        if os.path.exists(cache_file):
            mtime = os.path.getmtime(cache_file)
            date_str = time.strftime("%Y-%m-%d %H:%M", time.localtime(mtime))
            size_mb = os.path.getsize(cache_file) / (1024 * 1024)
            return {
                "exists": True,
                "date": date_str,
                "size_mb": round(size_mb, 1),
                "num_taxa": len(self.code_to_parent) if self.loaded else "unknown",
            }
        return {"exists": False}

    # ---- Internal methods ----

    def _parse_names(self, names_file, progress_callback=None):
        """Parse an NCBI names.dmp file, plain or gzipped."""
        with open_text(names_file) as f:
            self._parse_names_stream(f, progress_callback)

    def _parse_nodes(self, nodes_file, progress_callback=None):
        """Parse an NCBI nodes.dmp file, plain or gzipped."""
        with open_text(nodes_file) as f:
            self._parse_nodes_stream(f, progress_callback)

    def _parse_names_stream(self, stream, progress_callback=None):
        """Parse names.dmp content from any open text stream.

        names.dmp holds ~4.9M rows, only 2.9M of which are scientific names;
        the rest are synonyms, authorities ("Bacteria Latreille et al. 1825")
        and strain designations. Letting those overwrite each other in file
        order makes "Bacteria" resolve to a stick-insect genus instead of the
        domain, so scientific names win and true homonyms are recorded.
        """
        self.name_to_code = {}
        self.code_to_names = {}
        self.code_to_scientific = {}
        self.homonyms = {}

        if progress_callback:
            progress_callback("Parsing names.dmp...")

        count = 0
        for line in stream:
            count += 1
            if count % 100000 == 0 and progress_callback:
                progress_callback(f"Parsing names.dmp... {count:,} entries")

            parts = line.split("|")
            if len(parts) < 4:
                continue

            tax_id = parts[0].strip()
            name = parts[1].strip()
            name_class = parts[3].strip()
            is_scientific = name_class == "scientific name"

            if tax_id not in self.code_to_names:
                self.code_to_names[tax_id] = []
            self.code_to_names[tax_id].append(name)

            if is_scientific:
                self.code_to_scientific[tax_id] = name

            key = normalize_name(name)
            if not key:
                continue

            held_by = self.name_to_code.get(key)
            if held_by is None:
                self.name_to_code[key] = tax_id
            elif held_by == tax_id:
                pass
            elif is_scientific:
                # names.dmp lists each taxon's scientific name first, so an
                # incumbent's scientific name is already known by now.
                if self.code_to_scientific.get(held_by) != name:
                    self.name_to_code[key] = tax_id       # promote over a synonym
                else:
                    self.homonyms.setdefault(key, [held_by]).append(tax_id)

        if progress_callback:
            progress_callback(
                f"Parsed {count:,} name entries "
                f"({len(self.homonyms):,} ambiguous names)")

    def _parse_nodes_stream(self, stream, progress_callback=None):
        """Parse nodes.dmp content from any open text stream."""
        self.code_to_parent = {}
        self.code_to_rank = {}

        if progress_callback:
            progress_callback("Parsing nodes.dmp...")

        count = 0
        for line in stream:
            count += 1
            if count % 100000 == 0 and progress_callback:
                progress_callback(f"Parsing nodes.dmp... {count:,} entries")

            parts = line.split("|")
            if len(parts) < 3:
                continue

            # Remove whitespace from fields
            tax_id = parts[0].strip()
            parent_id = parts[1].strip()
            rank = parts[2].strip()

            self.code_to_parent[tax_id] = parent_id
            self.code_to_rank[tax_id] = rank

        self.loaded = True
        if progress_callback:
            progress_callback(f"Parsed {count:,} node entries")

    def _save_cache(self, progress_callback=None):
        """Save parsed taxonomy to pickle cache."""
        cache_file = os.path.join(self.cache_dir, "taxonomy.pkl")
        if progress_callback:
            progress_callback("Saving taxonomy cache...")

        data = {
            "version": self._cache_version,
            "name_to_code": self.name_to_code,
            "code_to_names": self.code_to_names,
            "code_to_scientific": self.code_to_scientific,
            "code_to_parent": self.code_to_parent,
            "code_to_rank": self.code_to_rank,
            "homonyms": self.homonyms,
        }
        with open(cache_file, "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

        if progress_callback:
            size_mb = os.path.getsize(cache_file) / (1024 * 1024)
            progress_callback(f"Cache saved ({size_mb:.1f} MB)")

    def _load_cache(self, progress_callback=None):
        """Load taxonomy from pickle cache."""
        cache_file = os.path.join(self.cache_dir, "taxonomy.pkl")
        if not os.path.exists(cache_file):
            return False

        if progress_callback:
            progress_callback("Loading taxonomy from cache...")

        try:
            with open(cache_file, "rb") as f:
                data = pickle.load(f)

            if data.get("version") != self._cache_version:
                if progress_callback:
                    progress_callback("Cache version mismatch, need to rebuild")
                return False

            self.name_to_code = data["name_to_code"]
            self.code_to_names = data["code_to_names"]
            self.code_to_scientific = data["code_to_scientific"]
            self.code_to_parent = data["code_to_parent"]
            self.code_to_rank = data["code_to_rank"]
            self.homonyms = data.get("homonyms", {})
            self.loaded = True

            if progress_callback:
                progress_callback(
                    f"Loaded {len(self.code_to_parent):,} taxa from cache"
                )
            return True

        except Exception as e:
            if progress_callback:
                progress_callback(f"Cache load failed: {e}")
            return False

    def _download_taxonomy(self, dest_path, progress_callback=None):
        """Download NCBI taxonomy dump."""
        if progress_callback:
            progress_callback("Downloading NCBI taxonomy dump...")

        def reporthook(block_num, block_size, total_size):
            if progress_callback and total_size > 0:
                downloaded = block_num * block_size
                pct = min(100, int(downloaded * 100 / total_size))
                mb = downloaded / (1024 * 1024)
                progress_callback(f"Downloading... {mb:.1f} MB ({pct}%)")

        urllib.request.urlretrieve(NCBI_TAXONOMY_URL, dest_path, reporthook)

        if progress_callback:
            progress_callback("Download complete")
