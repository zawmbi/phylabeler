#!/usr/bin/env python3
"""
test_archive_io.py - Tests for compressed and archived input handling.

Covers every form the NCBI taxonomy and tree files can arrive in: plain,
gzipped, tar/tgz/zip, and nested or misnamed archives. Builds its own
fixtures in a temp directory, so it needs no data files and no test runner:

    python test_archive_io.py
"""

import gzip
import os
import shutil
import sys
import tarfile
import tempfile
import zipfile

import archive_io
from archive_io import ArchiveError, is_archive, probe
from taxonomy_db import TaxonomyDB
from tree_parser import is_tree_filename, parse_newick_file

# A tiny synthetic taxonomy:
#   root > Metazoa > Chordata  > Homo        > {sapiens, neanderthalensis}
#                  > Arthropoda > Drosophila > melanogaster
NODES = [("1", "1", "no rank"), ("33208", "1", "kingdom"),
         ("7711", "33208", "phylum"), ("9605", "7711", "genus"),
         ("9606", "9605", "species"), ("63221", "9605", "species"),
         ("6656", "33208", "phylum"), ("7215", "6656", "genus"),
         ("7227", "7215", "species")]
NAMES = [("1", "root"), ("33208", "Metazoa"), ("7711", "Chordata"),
         ("9605", "Homo"), ("9606", "Homo sapiens"),
         ("63221", "Homo sapiens neanderthalensis"), ("6656", "Arthropoda"),
         ("7215", "Drosophila"), ("7227", "Drosophila melanogaster")]

NAMES_DMP = "".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n"
                    for t, n in NAMES)
NODES_DMP = "".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t, p, r in NODES)
TREE = ("((Homo_sapiens:0.1,Homo_sapiens_neanderthalensis:0.1)80:0.3,"
        "Drosophila_melanogaster:0.9);")
HUMAN_LINEAGE = ["9606", "9605", "7711", "33208", "1"]

_passed = _failed = 0


def check(label, got, want):
    global _passed, _failed
    if got == want:
        _passed += 1
        print(f"  ok    {label}")
    else:
        _failed += 1
        print(f"  FAIL  {label}\n          got  {got!r}\n          want {want!r}")


def build_fixtures(root):
    """Write every input packaging variant into root."""
    def path(*parts):
        return os.path.join(root, *parts)

    os.makedirs(path("plain"))
    with open(path("plain", "names.dmp"), "w") as f:
        f.write(NAMES_DMP)
    with open(path("plain", "nodes.dmp"), "w") as f:
        f.write(NODES_DMP)
    with gzip.open(path("names.dmp.gz"), "wt") as f:
        f.write(NAMES_DMP)
    with gzip.open(path("nodes.dmp.gz"), "wt") as f:
        f.write(NODES_DMP)

    # Real dumps carry more than the two files we care about.
    members = {"names.dmp": NAMES_DMP, "nodes.dmp": NODES_DMP,
               "division.dmp": "0\t|\tUNA\t|\n"}

    def write_tar(name, mode, prefix=""):
        with tarfile.open(path(name), mode) as tar:
            for member, data in members.items():
                scratch = path("_scratch")
                with open(scratch, "w") as f:
                    f.write(data)
                tar.add(scratch, arcname=prefix + member)
                os.remove(scratch)

    write_tar("taxdump.tar.gz", "w:gz")
    write_tar("taxdump.tgz", "w:gz")
    write_tar("taxdump.tar", "w")
    write_tar("new_taxdump.tar.gz", "w:gz", prefix="new_taxdump/")
    # A tar.gz saved under a bare .gz name: the extension lies, magic doesn't.
    shutil.copy(path("taxdump.tar.gz"), path("renamed.gz"))

    with zipfile.ZipFile(path("taxdmp.zip"), "w", zipfile.ZIP_DEFLATED) as zf:
        for member, data in members.items():
            zf.writestr(member, data)

    with tarfile.open(path("wrong.tar.gz"), "w:gz") as tar:
        scratch = path("_scratch")
        with open(scratch, "w") as f:
            f.write("not a taxonomy dump\n")
        tar.add(scratch, arcname="readme.txt")
        os.remove(scratch)

    os.makedirs(path("trees"))
    with open(path("trees", "gene1.nwk"), "w") as f:
        f.write(TREE + "\n")
    with gzip.open(path("trees", "gene2.nwk.gz"), "wt") as f:
        f.write(TREE + "\n")
    with open(path("trees", "notes.txt"), "w") as f:
        f.write("ignore me\n")
    with tarfile.open(path("trees.tgz"), "w:gz") as tar:
        tar.add(path("trees"), arcname="trees")

    return path


def test_format_detection(path):
    print("\nFormat detection is by magic bytes, not extension")
    check("taxdump.tar.gz", probe(path("taxdump.tar.gz")), ("gzip", "tar"))
    check("taxdump.tgz", probe(path("taxdump.tgz")), ("gzip", "tar"))
    check("taxdump.tar", probe(path("taxdump.tar")), (None, "tar"))
    check("taxdmp.zip", probe(path("taxdmp.zip")), (None, "zip"))
    check("names.dmp.gz", probe(path("names.dmp.gz")), ("gzip", None))
    check("names.dmp", probe(path("plain", "names.dmp")), (None, None))
    check("tar.gz renamed to .gz", probe(path("renamed.gz")), ("gzip", "tar"))
    check("is_archive(names.dmp.gz)", is_archive(path("names.dmp.gz")), False)
    check("is_archive(renamed.gz)", is_archive(path("renamed.gz")), True)


def test_taxonomy_loading(path):
    print("\nTaxonomy loads from every packaging")

    def load(**kwargs):
        with tempfile.TemporaryDirectory() as cache:
            db = TaxonomyDB(cache_dir=cache)
            db.load(**kwargs)
            return db

    cases = [
        ("plain dmp pair", dict(names_file=path("plain", "names.dmp"),
                                nodes_file=path("plain", "nodes.dmp"))),
        ("gzipped dmp pair", dict(names_file=path("names.dmp.gz"),
                                  nodes_file=path("nodes.dmp.gz"))),
        ("mixed plain and gz", dict(names_file=path("plain", "names.dmp"),
                                    nodes_file=path("nodes.dmp.gz"))),
        ("archive taxdump.tar.gz", dict(archive=path("taxdump.tar.gz"))),
        ("archive taxdump.tgz", dict(archive=path("taxdump.tgz"))),
        ("archive taxdump.tar", dict(archive=path("taxdump.tar"))),
        ("archive taxdmp.zip", dict(archive=path("taxdmp.zip"))),
        ("archive nested new_taxdump",
         dict(archive=path("new_taxdump.tar.gz"))),
        ("archive given as names_file",
         dict(names_file=path("taxdump.tar.gz"))),
        ("renamed .gz given as names_file", dict(names_file=path("renamed.gz"))),
    ]
    for label, kwargs in cases:
        db = load(**kwargs)
        check(f"{label}: taxa", len(db.code_to_parent), len(NODES))
        check(f"{label}: human lineage",
              db.get_lineage_codes(db.lookup_name("Homo_sapiens")),
              HUMAN_LINEAGE)


def test_cache_round_trip(path):
    print("\nArchive load caches without unpacking dmp files")
    with tempfile.TemporaryDirectory() as cache:
        TaxonomyDB(cache_dir=cache).load(archive=path("taxdump.tar.gz"))
        check("cache written",
              os.path.exists(os.path.join(cache, "taxonomy.pkl")), True)
        check("no dmp files left behind",
              [n for n in os.listdir(cache) if n.endswith(".dmp")], [])

        reloaded = TaxonomyDB(cache_dir=cache)
        check("reload from cache", reloaded.load(), True)
        check("cached taxa", len(reloaded.code_to_parent), len(NODES))


def test_error_paths(path):
    print("\nBad input fails with a useful message")
    try:
        with tempfile.TemporaryDirectory() as cache:
            TaxonomyDB(cache_dir=cache).load(archive=path("wrong.tar.gz"))
        check("archive without dmp files raises", False, True)
    except ArchiveError as e:
        check("archive without dmp files raises ArchiveError", True, True)
        check("message names the missing files",
              "names.dmp" in str(e) and "nodes.dmp" in str(e), True)

    try:
        archive_io.open_text(path("taxdump.tar.gz"))
        check("open_text rejects archives", False, True)
    except ArchiveError:
        check("open_text rejects archives", True, True)


def test_tree_inputs(path):
    print("\nTree files, compressed and archived")
    check("gene1.nwk", is_tree_filename("gene1.nwk"), True)
    check("gene2.nwk.gz", is_tree_filename("gene2.nwk.gz"), True)
    check("x.tre.bz2", is_tree_filename("x.tre.bz2"), True)
    check("x.tre.gz", is_tree_filename("x.tre.gz"), True)
    check("macOS ._x.tre.gz is not a tree", is_tree_filename("._x.tre.gz"),
          False)
    check("notes.txt", is_tree_filename("notes.txt"), False)
    check("archive.tar.gz is not a tree", is_tree_filename("archive.tar.gz"),
          False)

    plain = parse_newick_file(path("trees", "gene1.nwk"))
    gzipped = parse_newick_file(path("trees", "gene2.nwk.gz"))
    check("gzipped tree parses identically",
          gzipped.to_newick(True), plain.to_newick(True))
    check("tip count", gzipped.count_tips(), 3)

    with tempfile.TemporaryDirectory() as work:
        files = archive_io.extract_matching(path("trees.tgz"), work,
                                            is_tree_filename)
        check("only tree members extracted",
              sorted(os.path.basename(f) for f in files),
              ["gene1.nwk", "gene2.nwk.gz"])
        check("extracted flat into workdir",
              all(os.path.dirname(f) == work for f in files), True)
        check("extracted trees parse",
              [parse_newick_file(f).count_tips() for f in sorted(files)],
              [3, 3])


def main():
    with tempfile.TemporaryDirectory(prefix="phylabeler_test_") as root:
        path = build_fixtures(root)
        test_format_detection(path)
        test_taxonomy_loading(path)
        test_cache_round_trip(path)
        test_error_paths(path)
        test_tree_inputs(path)

    print(f"\n{_passed} passed, {_failed} failed")
    return 1 if _failed else 0


if __name__ == "__main__":
    sys.exit(main())
