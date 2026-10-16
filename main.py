#!/usr/bin/env python3
"""
PhyLabeler - Phylogenetic Tree Labeler

A tool for bioinformaticians to check gene trees against NCBI taxonomy
and identify monophyletic / non-monophyletic clades.

Usage:
    python main.py                           # Launch GUI
    python main.py --gui                     # Launch GUI (explicit)
    python main.py --tree FILE --names FILE --nodes FILE   # CLI mode
    python main.py --tree FILE --taxdump taxdump.tar.gz    # Un-extracted dump
    python main.py --tree FILE --cache       # CLI mode with cached taxonomy
    python main.py --download                # Download NCBI taxonomy and cache
    python main.py --batch DIR --cache       # Batch process all trees in directory
    python main.py --batch trees.tgz --cache # Batch process trees in an archive

Tree files and dmp files may be gzipped; the NCBI taxonomy may be given as
the un-extracted taxdump.tar.gz / .tgz / taxdmp.zip.

Based on LabelPhy by J.F. Walker (https://github.com/jfwalker/LabelPhy)
PeerJ Computer Science 56 (2016)

Pure Python - no external dependencies.
"""

import argparse
import os
import shutil
import sys
import tempfile

from archive_io import (ArchiveError, extract_matching, is_archive,
                        strip_compression_suffix)
from taxonomy_db import TaxonomyDB
from tree_parser import is_tree_filename, parse_newick_file
from monophyly import (CONFLICT_CAUSES_NOTE, DEFAULT_MIN_BRANCH_LENGTH,
                       MonophylyChecker, is_species_level)


def cli_progress(msg):
    """Print progress to stderr."""
    sys.stderr.write(f"\r{msg}                    ")
    sys.stderr.flush()


def run_cli(args):
    """Run the command-line workflow."""
    db = TaxonomyDB(cache_dir=args.cache_dir)

    # Handle taxonomy download
    if args.download:
        print("Downloading NCBI taxonomy dump...")
        db.download_and_load(progress_callback=cli_progress)
        print(f"\nTaxonomy cached: {len(db.code_to_parent):,} taxa")
        if not args.tree and not args.batch:
            return

    # Load taxonomy
    try:
        if args.taxdump:
            db.load(archive=args.taxdump, progress_callback=cli_progress)
        elif args.cache:
            if not db.load(progress_callback=cli_progress):
                print("\nNo cached taxonomy found. Use --download or provide "
                      "--names and --nodes files.", file=sys.stderr)
                sys.exit(1)
        elif args.names and (args.nodes or is_archive(args.names)):
            db.load(names_file=args.names, nodes_file=args.nodes,
                    progress_callback=cli_progress)
        elif not args.download:
            # Try cache as fallback
            if not db.load(progress_callback=cli_progress):
                print("No taxonomy data available. Options:\n"
                      "  --download         Download from NCBI\n"
                      "  --cache            Use cached data\n"
                      "  --names + --nodes  Provide local files\n"
                      "  --taxdump          Provide an un-extracted NCBI dump",
                      file=sys.stderr)
                sys.exit(1)
    except (ArchiveError, OSError) as e:
        print(f"\nCould not load taxonomy: {e}", file=sys.stderr)
        sys.exit(1)

    sys.stderr.write("\n")
    print(f"Taxonomy loaded: {len(db.code_to_parent):,} taxa",
          file=sys.stderr)

    checker = MonophylyChecker(db, min_branch_length=args.min_branch_length,
                               min_support=args.min_support)

    # Single tree mode
    if args.tree:
        _process_single(args, checker)

    # Batch mode
    elif args.batch:
        _process_batch(args, checker)


def _process_single(args, checker):
    """Process a single tree file."""
    tree = parse_newick_file(args.tree)
    print(f"Tree: {args.tree} ({tree.count_tips()} tips)",
          file=sys.stderr)

    analysis = checker.analyze(tree, progress_callback=cli_progress)
    sys.stderr.write("\n")

    results = analysis["node_results"]
    taxon_results = analysis["taxon_results"]
    unresolved = analysis["unresolved"]

    checker.label_tree(tree, results)

    outprefix = args.outfile

    tree_out = outprefix + ".labeled.nwk"
    with open(tree_out, "w", encoding="utf-8") as f:
        f.write(tree.to_newick(branch_lengths=True,
                               label_format=args.label_format) + ";\n")
    print(f"Labeled tree written to: {tree_out}", file=sys.stderr)

    report_out = outprefix + ".report.txt"
    report = checker.get_summary(results, taxon_results, analysis["matches"])
    if unresolved:
        report += "\n\n--- Unresolved Tips ---\n"
        for tip in sorted(unresolved):
            report += f"  {tip}\n"
    with open(report_out, "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(f"Report written to: {report_out}", file=sys.stderr)

    groups = analysis["groups"]
    broken = [g for g in groups if g.is_broken]
    weak_only = sum(1 for g in groups if g.is_unsupported_conflict)
    species = sum(1 for g in broken if is_species_level(g.rank))
    print(f"\nResults: {len(taxon_results)} taxonomic groups with 2+ tips "
          f"({len(groups)} distinct), {len(broken)} not monophyletic "
          f"({species} at species level or below), {weak_only} contradicted "
          f"only by weak branches, {len(unresolved)} unresolved tips")

    # Lead with the tips behind most conflicts: on a real tree one misplaced
    # sample often accounts for several broken groups.
    show = checker.db.display_label
    repairing = [t for t in analysis["conflicting_tips"] if t.repairs]
    if repairing:
        fixable = sum(len(t.repairs) for t in repairing)
        print(f"\nTips behind the most conflicts "
              f"({fixable} of {len(broken)} supported conflicts):")
        for tip in repairing[:5]:
            names = ", ".join(g.name for g in tip.repairs[:3])
            if len(tip.repairs) > 3:
                names += f", +{len(tip.repairs) - 3} more"
            print(f"  {show(tip.label):<32} repairs {len(tip.repairs):>2} "
                  f"({tip.role}): {names}")
        if len(repairing) > 5:
            print(f"  ... {len(repairing) - 5} more single-group tips")
        print("  " + CONFLICT_CAUSES_NOTE)

    above = [g for g in broken if not is_species_level(g.rank)]
    if above:
        print("\nConflicting groups above species level (cheapest fix first):")
        for g in sorted(above, key=lambda x: x.severity_sort_key())[:15]:
            blame = show(g.culprit) if g.culprit else "no single tip"
            print(f"  move {g.conflict_cost:<4} {g.name} "
                  f"({g.rank}, {g.member_count} members)"
                  f"  <- {blame}")
        if len(above) > 15:
            print(f"  ... {len(above) - 15} more; full detail in {report_out}")


def _collect_batch_trees(source, workdir):
    """Return tree file paths from a directory or an archive of trees.

    Archive members are unpacked into workdir so the rest of the batch code
    can treat every input as an ordinary file path.
    """
    if os.path.isdir(source):
        return [os.path.join(source, name)
                for name in os.listdir(source)
                if is_tree_filename(name)]

    if is_archive(source):
        files = extract_matching(source, workdir, is_tree_filename)
        print(f"Extracted {len(files)} tree file(s) from {source}",
              file=sys.stderr)
        return files

    print(f"{source} is neither a directory nor a tree archive",
          file=sys.stderr)
    sys.exit(1)


def _process_batch(args, checker):
    """Process all tree files in a directory or archive."""
    workdir = tempfile.mkdtemp(prefix="phylabeler_batch_")
    try:
        _run_batch(args, checker, _collect_batch_trees(args.batch, workdir))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _run_batch(args, checker, tree_files):
    """Analyze a resolved list of tree files and write a batch summary."""
    if not tree_files:
        print(f"No tree files found in {args.batch}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tree_files)} tree files", file=sys.stderr)

    outdir = args.outfile if os.path.isdir(args.outfile) else os.path.dirname(
        args.outfile) or "."
    os.makedirs(outdir, exist_ok=True)

    batch_summary = []
    for i, fp in enumerate(sorted(tree_files)):
        name = strip_compression_suffix(os.path.basename(fp))
        basename = os.path.splitext(name)[0]
        print(f"\n[{i + 1}/{len(tree_files)}] {basename}",
              file=sys.stderr)

        try:
            tree = parse_newick_file(fp)
            analysis = checker.analyze(tree)
            results = analysis["node_results"]
            taxa = analysis["taxon_results"]
            checker.label_tree(tree, results)

            # Write labeled tree
            out_tree = os.path.join(outdir, basename + ".labeled.nwk")
            with open(out_tree, "w", encoding="utf-8") as f:
                f.write(tree.to_newick(branch_lengths=True,
                                       label_format=args.label_format) + ";\n")

            # Per-tree report, so batch mode gives the same detail as
            # single-tree mode instead of only a one-line summary.
            out_report = os.path.join(outdir, basename + ".report.txt")
            with open(out_report, "w", encoding="utf-8") as f:
                f.write(checker.get_summary(
                    results, taxa, analysis["matches"]) + "\n")

            broken = [g for g in analysis["groups"] if g.is_broken]
            repairing = [t for t in analysis["conflicting_tips"] if t.repairs]
            # Across many gene trees, "which sequence keeps causing this?"
            # is the question worth answering per tree.
            worst = (f", worst tip {checker.db.display_label(repairing[0].label)} "
                     f"({len(repairing[0].repairs)})" if repairing else "")
            batch_summary.append(
                f"{basename}: {tree.count_tips()} tips, "
                f"{len(analysis['groups'])} groups, {len(broken)} broken, "
                f"{len(analysis['unresolved'])} unresolved{worst}"
            )
        except Exception as e:
            batch_summary.append(f"{basename}: ERROR - {e}")
            print(f"  Error: {e}", file=sys.stderr)

    # Write batch summary
    summary_path = os.path.join(outdir, "batch_summary.txt")
    with open(summary_path, "w") as f:
        f.write("PhyLabeler Batch Analysis Summary\n")
        f.write("=" * 40 + "\n\n")
        for line in batch_summary:
            f.write(line + "\n")

    print(f"\nBatch complete. Summary: {summary_path}")
    for line in batch_summary:
        print(f"  {line}")


def run_gui():
    """Launch the Tkinter GUI."""
    from gui import main as gui_main
    gui_main()


def main():
    parser = argparse.ArgumentParser(
        description="PhyLabeler - Check gene trees against NCBI taxonomy",
        epilog="Run without arguments to launch the GUI."
    )

    # Mode selection
    parser.add_argument("--gui", action="store_true",
                        help="Launch graphical interface (default if no args)")

    # Input
    parser.add_argument("--tree", type=str,
                        help="Input tree file (Newick, optionally gzipped)")
    parser.add_argument("--batch", type=str,
                        help="Directory or archive (.tar.gz/.tgz/.zip) of "
                             "tree files for batch processing")

    # Taxonomy sources
    parser.add_argument("--names", type=str,
                        help="NCBI names.dmp file (may be gzipped, or be an "
                             "un-extracted taxdump archive)")
    parser.add_argument("--nodes", type=str,
                        help="NCBI nodes.dmp file (may be gzipped)")
    parser.add_argument("--taxdump", type=str,
                        help="Un-extracted NCBI taxonomy archive "
                             "(taxdump.tar.gz, taxdump.tgz, taxdmp.zip)")
    parser.add_argument("--cache", action="store_true",
                        help="Use cached taxonomy database")
    parser.add_argument("--download", action="store_true",
                        help="Download NCBI taxonomy and cache it")
    parser.add_argument("--cache-dir", type=str, dest="cache_dir",
                        default=None,
                        help="Directory for taxonomy cache "
                             "(default: ~/.phylabeler_cache)")

    # Output
    parser.add_argument("--outfile", type=str, default="out",
                        help="Output prefix (or directory for batch mode)")
    parser.add_argument("--label-format", type=str, dest="label_format",
                        default="taxon_support",
                        choices=["taxon_support", "taxon", "support"],
                        help="What internal nodes carry in the labeled tree: "
                             "taxon_support writes Anthocerotophyta/100 and "
                             "keeps the input's bootstrap values (default)")

    # Branch support
    parser.add_argument("--min-branch-length", type=float,
                        dest="min_branch_length",
                        default=DEFAULT_MIN_BRANCH_LENGTH,
                        help="Internal branches this short or shorter are "
                             "treated as unresolved, so they never count "
                             "against a group (default: %(default)g; "
                             "-1 disables)")
    parser.add_argument("--min-support", type=float, dest="min_support",
                        default=None,
                        help="Treat branches with support below this as "
                             "unresolved, on the tree's own scale (e.g. 70 "
                             "for bootstrap, 0.95 for posteriors)")

    args = parser.parse_args()
    if args.min_branch_length is not None and args.min_branch_length < 0:
        args.min_branch_length = None

    # Decide CLI vs GUI
    has_cli_args = (args.tree or args.batch or args.download or
                    args.names or args.nodes or args.taxdump)

    if args.gui or (not has_cli_args and len(sys.argv) == 1):
        run_gui()
    else:
        run_cli(args)


if __name__ == "__main__":
    main()
