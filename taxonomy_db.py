"""
taxonomy_db.py - NCBI Taxonomy Database with Caching

Downloads, parses, and caches NCBI taxonomy data (names.dmp, nodes.dmp).
Uses pickle for caching so the GenBank tree doesn't need to be regenerated
every time. Pure Python stdlib - no external dependencies.
"""

import os
import sys
import pickle
import urllib.request
import tarfile
import io
import time

DEFAULT_CACHE_DIR = os.path.join(os.path.expanduser("~"), ".phylabeler_cache")
NCBI_TAXONOMY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


class TaxonomyDB:
    """Cached NCBI taxonomy database.

    Parses names.dmp and nodes.dmp into efficient lookup structures,
    then caches them as pickle files so subsequent runs are fast.
    """

    def __init__(self, cache_dir=None):
        self.cache_dir = cache_dir or DEFAULT_CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)

        # Core data structures
        self.name_to_code = {}       # name -> tax_id
        self.code_to_names = {}      # tax_id -> [list of names]
        self.code_to_scientific = {}  # tax_id -> scientific name
        self.code_to_parent = {}     # tax_id -> parent tax_id
        self.code_to_rank = {}       # tax_id -> rank string
        self.loaded = False
        self._cache_version = 2      # Bump to invalidate old caches

    # ---- Public API ----

    def load(self, names_file=None, nodes_file=None, progress_callback=None):
        """Load taxonomy data from cache or source files.

        If a valid cache exists and no source files are given, loads from cache.
        If source files are provided, parses them and updates the cache.
        """
        if names_file is None and nodes_file is None:
            if self._load_cache(progress_callback):
                return True
            return False

        if names_file and nodes_file:
            self._parse_names(names_file, progress_callback)
            self._parse_nodes(nodes_file, progress_callback)
            self._save_cache(progress_callback)
            self.loaded = True
            return True

        return False

    def download_and_load(self, progress_callback=None):
        """Download NCBI taxonomy dump and load it."""
        tar_path = os.path.join(self.cache_dir, "taxdump.tar.gz")
        names_path = os.path.join(self.cache_dir, "names.dmp")
        nodes_path = os.path.join(self.cache_dir, "nodes.dmp")

        self._download_taxonomy(tar_path, progress_callback)
        self._extract_taxonomy(tar_path, progress_callback)
        self.load(names_file=names_path, nodes_file=nodes_path,
                  progress_callback=progress_callback)

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
        """Look up a species name and return its tax_id, or None."""
        # Exact match
        if name in self.name_to_code:
            return self.name_to_code[name]
        # Try replacing underscores with spaces
        spaced = name.replace("_", " ")
        if spaced in self.name_to_code:
            return self.name_to_code[spaced]
        # Case-insensitive search
        lower = spaced.lower()
        for n, code in self.name_to_code.items():
            if n.lower() == lower:
                return code
        return None

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
        """Parse NCBI names.dmp file."""
        self.name_to_code = {}
        self.code_to_names = {}
        self.code_to_scientific = {}

        if progress_callback:
            progress_callback("Parsing names.dmp...")

        with open(names_file, "r", encoding="utf-8", errors="replace") as f:
            count = 0
            for line in f:
                count += 1
                if count % 100000 == 0 and progress_callback:
                    progress_callback(f"Parsing names.dmp... {count:,} entries")

                parts = line.split("|")
                if len(parts) < 4:
                    continue

                tax_id = parts[0].strip()
                name = parts[1].strip()
                name_class = parts[3].strip()

                self.name_to_code[name] = tax_id

                if tax_id not in self.code_to_names:
                    self.code_to_names[tax_id] = []
                self.code_to_names[tax_id].append(name)

                if name_class == "scientific name":
                    self.code_to_scientific[tax_id] = name

        if progress_callback:
            progress_callback(f"Parsed {count:,} name entries")

    def _parse_nodes(self, nodes_file, progress_callback=None):
        """Parse NCBI nodes.dmp file."""
        self.code_to_parent = {}
        self.code_to_rank = {}

        if progress_callback:
            progress_callback("Parsing nodes.dmp...")

        with open(nodes_file, "r", encoding="utf-8", errors="replace") as f:
            count = 0
            for line in f:
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

    def _extract_taxonomy(self, tar_path, progress_callback=None):
        """Extract names.dmp and nodes.dmp from taxdump.tar.gz."""
        if progress_callback:
            progress_callback("Extracting taxonomy files...")

        with tarfile.open(tar_path, "r:gz") as tar:
            for member in tar.getmembers():
                if member.name in ("names.dmp", "nodes.dmp"):
                    tar.extract(member, self.cache_dir)

        if progress_callback:
            progress_callback("Extraction complete")
