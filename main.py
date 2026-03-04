#!/usr/bin/env python3
"""
PhyLabeler - Phylogenetic Tree Labeler

A tool for bioinformaticians to check gene trees against NCBI taxonomy
and identify monophyletic / non-monophyletic clades.

Usage:
    python main.py                           # Launch GUI
    python main.py --gui                     # Launch GUI (explicit)
    python main.py --tree FILE --names FILE --nodes FILE   # CLI mode
    python main.py --tree FILE --cache       # CLI mode with cached taxonomy
    python main.py --download                # Download NCBI taxonomy and cache
    python main.py --batch DIR --cache       # Batch process all trees in directory

Based on LabelPhy by J.F. Walker (https://github.com/jfwalker/LabelPhy)
PeerJ Computer Science 56 (2016)

Pure Python - no external dependencies.
"""

import argparse
import os
import sys
import glob

from taxonomy_db import TaxonomyDB
from tree_parser import parse_newick_file
from monophyly import MonophylyChecker


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
    if args.cache:
        if not db.load(progress_callback=cli_progress):
            print("\nNo cached taxonomy found. Use --download or provide "
                  "--names and --nodes files.", file=sys.stderr)
            sys.exit(1)
    elif args.names and args.nodes:
        db.load(names_file=args.names, nodes_file=args.nodes,
                progress_callback=cli_progress)
    elif not args.download:
        # Try cache as fallback
        if not db.load(progress_callback=cli_progress):
            print("No taxonomy data available. Options:\n"
                  "  --download         Download from NCBI\n"
                  "  --cache            Use cached data\n"
                  "  --names + --nodes  Provide local files",
                  file=sys.stderr)
            sys.exit(1)

    sys.stderr.write("\n")
    print(f"Taxonomy loaded: {len(db.code_to_parent):,} taxa",
          file=sys.stderr)

    checker = MonophylyChecker(db)

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

    results, unresolved = checker.check_tree(
        tree, progress_callback=cli_progress)
    sys.stderr.write("\n")

    # Label the tree
    checker.label_tree(tree, results)

    # Output
    outprefix = args.outfile

    # Write labeled tree
    tree_out = outprefix + ".labeled.nwk"
    with open(tree_out, "w") as f:
        f.write(tree.to_newick(branch_lengths=True) + ";\n")
    print(f"Labeled tree written to: {tree_out}", file=sys.stderr)

    # Write report
    report_out = outprefix + ".report.txt"
    report = checker.get_summary(results)
    if unresolved:
        report += "\n\n--- Unresolved Tips ---\n"
        for tip in sorted(unresolved):
            report += f"  {tip}\n"
    with open(report_out, "w") as f:
        f.write(report + "\n")
    print(f"Report written to: {report_out}", file=sys.stderr)

    # Print summary to stdout
    mono_count = sum(1 for r in results if r.status == "monophyletic")
    nonmono_count = sum(1 for r in results if r.status == "not_monophyletic")
    print(f"\nResults: {mono_count} monophyletic, "
          f"{nonmono_count} non-monophyletic, "
          f"{len(unresolved)} unresolved tips")

    # Print non-monophyletic clades (PI: "if doesn't match is preferred")
    non_mono = [r for r in results if r.status == "not_monophyletic"]
    if non_mono:
        print("\nNon-monophyletic clades:")
        for r in sorted(non_mono, key=lambda x: len(x.tip_labels),
                         reverse=True):
            intruders = ", ".join(r.intruders[:3])
            if len(r.intruders) > 3:
                intruders += f"... (+{len(r.intruders) - 3})"
            print(f"  {r.mrca_name} ({r.mrca_rank}) "
                  f"[{len(r.tip_labels)} tips] "
                  f"intruders: {intruders}")


def _process_batch(args, checker):
    """Process all tree files in a directory."""
    tree_dir = args.batch
    patterns = ["*.nwk", "*.tre", "*.tree", "*.newick"]
    tree_files = []
    for pat in patterns:
        tree_files.extend(glob.glob(os.path.join(tree_dir, pat)))

    if not tree_files:
        print(f"No tree files found in {tree_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tree_files)} tree files", file=sys.stderr)

    outdir = args.outfile if os.path.isdir(args.outfile) else os.path.dirname(
        args.outfile) or "."
    os.makedirs(outdir, exist_ok=True)

    batch_summary = []
    for i, fp in enumerate(sorted(tree_files)):
        basename = os.path.splitext(os.path.basename(fp))[0]
        print(f"\n[{i + 1}/{len(tree_files)}] {basename}",
              file=sys.stderr)

        try:
            tree = parse_newick_file(fp)
            results, unresolved = checker.check_tree(tree)
            checker.label_tree(tree, results)

            # Write labeled tree
            out_tree = os.path.join(outdir, basename + ".labeled.nwk")
            with open(out_tree, "w") as f:
                f.write(tree.to_newick(branch_lengths=True) + ";\n")

            mono = sum(1 for r in results if r.status == "monophyletic")
            nonmono = sum(1 for r in results
                          if r.status == "not_monophyletic")
            batch_summary.append(
                f"{basename}: {tree.count_tips()} tips, "
                f"{mono} mono, {nonmono} non-mono, "
                f"{len(unresolved)} unresolved"
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
                        help="Input tree file (Newick format)")
    parser.add_argument("--batch", type=str,
                        help="Directory of tree files for batch processing")

    # Taxonomy sources
    parser.add_argument("--names", type=str,
                        help="NCBI names.dmp file")
    parser.add_argument("--nodes", type=str,
                        help="NCBI nodes.dmp file")
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

    args = parser.parse_args()

    # Decide CLI vs GUI
    has_cli_args = (args.tree or args.batch or args.download or
                    args.names or args.nodes)

    if args.gui or (not has_cli_args and len(sys.argv) == 1):
        run_gui()
    else:
        run_cli(args)


if __name__ == "__main__":
    main()
