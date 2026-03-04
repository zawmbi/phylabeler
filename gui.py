"""
gui.py - Tkinter GUI for PhyLabeler

A graphical interface for bioinformaticians to:
- Load and visualize phylogenetic trees
- Manage NCBI taxonomy database (download, cache, update)
- Run monophyly analysis
- Browse results with monophyletic/non-monophyletic highlighting
- Batch process multiple gene trees
- Export labeled trees and reports

Pure Python stdlib (tkinter) - no external dependencies.
"""

import os
import sys
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

from taxonomy_db import TaxonomyDB
from tree_parser import parse_newick, parse_newick_file, Node
from monophyly import MonophylyChecker


# ---- Color scheme ----
COLORS = {
    "mono": "#2d8a4e",        # green for monophyletic
    "non_mono": "#d63031",    # red for non-monophyletic
    "unresolved": "#6c757d",  # gray for unresolved
    "tip": "#2c3e50",         # dark for tips
    "branch": "#555555",      # branch lines
    "bg": "#ffffff",          # canvas background
    "highlight": "#ffeaa7",   # selection highlight
    "header_bg": "#34495e",   # header background
    "header_fg": "#ecf0f1",   # header text
}


class PhyLabelerApp:
    """Main application window."""

    def __init__(self, root):
        self.root = root
        self.root.title("PhyLabeler - Phylogenetic Tree Labeler")
        self.root.geometry("1200x800")
        self.root.minsize(900, 600)

        # State
        self.db = TaxonomyDB()
        self.current_tree = None
        self.results = None
        self.unresolved_tips = []
        self.tree_files = []        # For batch mode
        self.batch_results = {}     # filename -> (results, unresolved)

        self._build_ui()
        self._check_cache()

    def _build_ui(self):
        """Build the complete UI layout."""
        # Menu bar
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open Tree...", command=self._open_tree,
                              accelerator="Ctrl+O")
        file_menu.add_command(label="Open Multiple Trees...",
                              command=self._open_batch)
        file_menu.add_separator()
        file_menu.add_command(label="Export Labeled Tree...",
                              command=self._export_tree)
        file_menu.add_command(label="Export Report...",
                              command=self._export_report)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self.root.quit,
                              accelerator="Ctrl+Q")
        menubar.add_cascade(label="File", menu=file_menu)

        db_menu = tk.Menu(menubar, tearoff=0)
        db_menu.add_command(label="Load from Files...",
                            command=self._load_taxonomy_files)
        db_menu.add_command(label="Download from NCBI...",
                            command=self._download_taxonomy)
        db_menu.add_separator()
        db_menu.add_command(label="Cache Info...", command=self._show_cache_info)
        menubar.add_cascade(label="Taxonomy DB", menu=db_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.bind("<Control-o>", lambda e: self._open_tree())
        self.root.bind("<Control-q>", lambda e: self.root.quit())

        # Main layout: left panel + right panel
        main_pane = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # ---- Left panel: controls + file list ----
        left_frame = ttk.Frame(main_pane, width=320)
        main_pane.add(left_frame, weight=1)

        # Taxonomy DB status
        db_frame = ttk.LabelFrame(left_frame, text="Taxonomy Database")
        db_frame.pack(fill=tk.X, padx=5, pady=5)

        self.db_status_var = tk.StringVar(value="Not loaded")
        ttk.Label(db_frame, textvariable=self.db_status_var).pack(
            padx=5, pady=2, anchor="w")

        db_btn_frame = ttk.Frame(db_frame)
        db_btn_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(db_btn_frame, text="Load Files",
                    command=self._load_taxonomy_files).pack(side=tk.LEFT, padx=2)
        ttk.Button(db_btn_frame, text="Download NCBI",
                    command=self._download_taxonomy).pack(side=tk.LEFT, padx=2)
        ttk.Button(db_btn_frame, text="Load Cache",
                    command=self._load_cache).pack(side=tk.LEFT, padx=2)

        # Tree files list
        tree_frame = ttk.LabelFrame(left_frame, text="Tree Files")
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        tree_btn_frame = ttk.Frame(tree_frame)
        tree_btn_frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(tree_btn_frame, text="Add Tree(s)",
                    command=self._open_tree).pack(side=tk.LEFT, padx=2)
        ttk.Button(tree_btn_frame, text="Clear",
                    command=self._clear_trees).pack(side=tk.LEFT, padx=2)

        self.tree_listbox = tk.Listbox(tree_frame, selectmode=tk.BROWSE)
        self.tree_listbox.pack(fill=tk.BOTH, expand=True, padx=5, pady=2)
        self.tree_listbox.bind("<<ListboxSelect>>", self._on_tree_select)

        # Analysis controls
        analysis_frame = ttk.LabelFrame(left_frame, text="Analysis")
        analysis_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(analysis_frame, text="Run Monophyly Check",
                    command=self._run_analysis).pack(fill=tk.X, padx=5, pady=2)
        ttk.Button(analysis_frame, text="Run Batch Analysis",
                    command=self._run_batch).pack(fill=tk.X, padx=5, pady=2)

        # ---- Right panel: tree view + results ----
        right_frame = ttk.Frame(main_pane)
        main_pane.add(right_frame, weight=3)

        right_pane = ttk.PanedWindow(right_frame, orient=tk.VERTICAL)
        right_pane.pack(fill=tk.BOTH, expand=True)

        # Tree visualization
        tree_view_frame = ttk.LabelFrame(right_pane, text="Tree Visualization")
        right_pane.add(tree_view_frame, weight=3)

        # Canvas with scrollbars for tree drawing
        canvas_frame = ttk.Frame(tree_view_frame)
        canvas_frame.pack(fill=tk.BOTH, expand=True)

        self.tree_canvas = tk.Canvas(canvas_frame, bg=COLORS["bg"],
                                     highlightthickness=0)
        h_scroll = ttk.Scrollbar(canvas_frame, orient=tk.HORIZONTAL,
                                  command=self.tree_canvas.xview)
        v_scroll = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL,
                                  command=self.tree_canvas.yview)
        self.tree_canvas.configure(xscrollcommand=h_scroll.set,
                                   yscrollcommand=v_scroll.set)

        self.tree_canvas.grid(row=0, column=0, sticky="nsew")
        v_scroll.grid(row=0, column=1, sticky="ns")
        h_scroll.grid(row=1, column=0, sticky="ew")
        canvas_frame.grid_rowconfigure(0, weight=1)
        canvas_frame.grid_columnconfigure(0, weight=1)

        # Zoom controls
        zoom_frame = ttk.Frame(tree_view_frame)
        zoom_frame.pack(fill=tk.X, padx=5, pady=2)
        self.zoom_var = tk.DoubleVar(value=1.0)
        ttk.Label(zoom_frame, text="Zoom:").pack(side=tk.LEFT)
        ttk.Scale(zoom_frame, from_=0.2, to=3.0, variable=self.zoom_var,
                  orient=tk.HORIZONTAL, command=lambda _: self._redraw_tree()
                  ).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        self.show_branch_lengths_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(zoom_frame, text="Branch lengths",
                         variable=self.show_branch_lengths_var,
                         command=self._redraw_tree).pack(side=tk.RIGHT)

        self.show_labels_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(zoom_frame, text="Labels",
                         variable=self.show_labels_var,
                         command=self._redraw_tree).pack(side=tk.RIGHT)

        # Results panel
        results_frame = ttk.LabelFrame(right_pane, text="Results")
        right_pane.add(results_frame, weight=2)

        # Notebook for different result views
        self.results_notebook = ttk.Notebook(results_frame)
        self.results_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Summary tab
        summary_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(summary_frame, text="Summary")
        self.summary_text = scrolledtext.ScrolledText(
            summary_frame, wrap=tk.WORD, font=("Courier", 10))
        self.summary_text.pack(fill=tk.BOTH, expand=True)

        # Non-monophyletic tab
        nonmono_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(nonmono_frame, text="Non-Monophyletic")
        self.nonmono_tree = ttk.Treeview(
            nonmono_frame,
            columns=("name", "rank", "tips", "intruders"),
            show="headings"
        )
        self.nonmono_tree.heading("name", text="Clade Name")
        self.nonmono_tree.heading("rank", text="Rank")
        self.nonmono_tree.heading("tips", text="# Tips")
        self.nonmono_tree.heading("intruders", text="Intruders")
        self.nonmono_tree.column("name", width=200)
        self.nonmono_tree.column("rank", width=100)
        self.nonmono_tree.column("tips", width=60)
        self.nonmono_tree.column("intruders", width=300)
        nonmono_scroll = ttk.Scrollbar(nonmono_frame, orient=tk.VERTICAL,
                                        command=self.nonmono_tree.yview)
        self.nonmono_tree.configure(yscrollcommand=nonmono_scroll.set)
        self.nonmono_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        nonmono_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Monophyletic tab
        mono_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(mono_frame, text="Monophyletic")
        self.mono_tree = ttk.Treeview(
            mono_frame,
            columns=("name", "rank", "tips"),
            show="headings"
        )
        self.mono_tree.heading("name", text="Clade Name")
        self.mono_tree.heading("rank", text="Rank")
        self.mono_tree.heading("tips", text="# Tips")
        self.mono_tree.column("name", width=250)
        self.mono_tree.column("rank", width=120)
        self.mono_tree.column("tips", width=80)
        mono_scroll = ttk.Scrollbar(mono_frame, orient=tk.VERTICAL,
                                     command=self.mono_tree.yview)
        self.mono_tree.configure(yscrollcommand=mono_scroll.set)
        self.mono_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        mono_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        # Unresolved tab
        unresolved_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(unresolved_frame, text="Unresolved Tips")
        self.unresolved_text = scrolledtext.ScrolledText(
            unresolved_frame, wrap=tk.WORD, font=("Courier", 10))
        self.unresolved_text.pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(self.root, textvariable=self.status_var,
                                relief=tk.SUNKEN, anchor="w")
        status_bar.pack(fill=tk.X, side=tk.BOTTOM, padx=5)

    # ---- Taxonomy DB operations ----

    def _check_cache(self):
        """Check if a cached taxonomy DB exists and load it."""
        info = self.db.get_cache_info()
        if info["exists"]:
            self.db_status_var.set(f"Cache found ({info['date']})")
            self._run_threaded(self._load_cache_threaded,
                               "Loading cached taxonomy...")
        else:
            self.db_status_var.set("No cache - load files or download from NCBI")

    def _load_cache(self):
        self._run_threaded(self._load_cache_threaded, "Loading taxonomy cache...")

    def _load_cache_threaded(self):
        success = self.db.load(progress_callback=self._update_status)
        if success:
            self.root.after(0, lambda: self.db_status_var.set(
                f"Loaded: {len(self.db.code_to_parent):,} taxa"))
        else:
            self.root.after(0, lambda: self.db_status_var.set(
                "No cache found"))

    def _load_taxonomy_files(self):
        """Load taxonomy from local names.dmp and nodes.dmp files."""
        names_file = filedialog.askopenfilename(
            title="Select names.dmp",
            filetypes=[("DMP files", "*.dmp"), ("All files", "*.*")]
        )
        if not names_file:
            return

        nodes_file = filedialog.askopenfilename(
            title="Select nodes.dmp",
            filetypes=[("DMP files", "*.dmp"), ("All files", "*.*")]
        )
        if not nodes_file:
            return

        def load():
            self.db.load(names_file=names_file, nodes_file=nodes_file,
                         progress_callback=self._update_status)
            self.root.after(0, lambda: self.db_status_var.set(
                f"Loaded: {len(self.db.code_to_parent):,} taxa"))

        self._run_threaded(load, "Loading taxonomy files...")

    def _download_taxonomy(self):
        """Download NCBI taxonomy dump."""
        if not messagebox.askyesno(
            "Download Taxonomy",
            "This will download the NCBI taxonomy dump (~60 MB compressed).\n"
            "It may take a few minutes. Continue?"
        ):
            return

        def download():
            self.db.download_and_load(progress_callback=self._update_status)
            self.root.after(0, lambda: self.db_status_var.set(
                f"Loaded: {len(self.db.code_to_parent):,} taxa"))

        self._run_threaded(download, "Downloading NCBI taxonomy...")

    def _show_cache_info(self):
        info = self.db.get_cache_info()
        if info["exists"]:
            msg = (f"Cache date: {info['date']}\n"
                   f"Cache size: {info['size_mb']} MB\n"
                   f"Taxa loaded: {info['num_taxa']}")
        else:
            msg = "No taxonomy cache found.\nLoad from files or download from NCBI."
        messagebox.showinfo("Cache Info", msg)

    # ---- Tree operations ----

    def _open_tree(self):
        """Open one or more tree files."""
        filepaths = filedialog.askopenfilenames(
            title="Select Tree File(s)",
            filetypes=[
                ("Newick files", "*.nwk *.tre *.tree *.newick"),
                ("All files", "*.*"),
            ]
        )
        for fp in filepaths:
            if fp not in self.tree_files:
                self.tree_files.append(fp)
                self.tree_listbox.insert(tk.END, os.path.basename(fp))

        if filepaths:
            self.tree_listbox.selection_clear(0, tk.END)
            idx = len(self.tree_files) - 1
            self.tree_listbox.selection_set(idx)
            self._load_tree(self.tree_files[idx])

    def _open_batch(self):
        """Open a directory of tree files."""
        dirpath = filedialog.askdirectory(title="Select Directory with Trees")
        if not dirpath:
            return
        for fname in sorted(os.listdir(dirpath)):
            if fname.endswith((".nwk", ".tre", ".tree", ".newick")):
                fp = os.path.join(dirpath, fname)
                if fp not in self.tree_files:
                    self.tree_files.append(fp)
                    self.tree_listbox.insert(tk.END, fname)

    def _clear_trees(self):
        self.tree_files.clear()
        self.tree_listbox.delete(0, tk.END)
        self.current_tree = None
        self.results = None
        self.batch_results.clear()
        self.tree_canvas.delete("all")
        self._clear_results()

    def _on_tree_select(self, event):
        sel = self.tree_listbox.curselection()
        if sel:
            self._load_tree(self.tree_files[sel[0]])

    def _load_tree(self, filepath):
        """Parse and display a tree file."""
        try:
            self.current_tree = parse_newick_file(filepath)
            tips = self.current_tree.count_tips()
            internal = self.current_tree.count_internal()
            self._update_status(
                f"Loaded: {os.path.basename(filepath)} "
                f"({tips} tips, {internal} internal nodes)"
            )
            self._draw_tree()

            # Show previous results if available
            basename = os.path.basename(filepath)
            if basename in self.batch_results:
                self.results, self.unresolved_tips = self.batch_results[basename]
                self._show_results()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse tree:\n{e}")

    # ---- Tree Drawing ----

    def _draw_tree(self):
        """Draw the current tree on the canvas."""
        if not self.current_tree:
            return

        self.tree_canvas.delete("all")
        zoom = self.zoom_var.get()

        # Calculate layout
        tip_spacing = 18 * zoom
        x_scale = 150 * zoom
        margin_left = 20
        margin_top = 20

        tips = self.current_tree.get_tip_nodes()
        n_tips = len(tips)

        # Assign y positions to tips
        y_positions = {}
        for i, tip in enumerate(tips):
            y_positions[id(tip)] = margin_top + i * tip_spacing

        # Calculate x positions based on depth (post-order)
        max_depth = self.current_tree.depth()
        if max_depth == 0:
            max_depth = 1

        def assign_positions(node, depth=0):
            x = margin_left + depth * (x_scale / max_depth)
            if node.is_tip:
                y = y_positions[id(node)]
            else:
                child_ys = []
                for child in node.children:
                    assign_positions(child, depth + 1)
                    child_ys.append(y_positions[id(child)])
                y = sum(child_ys) / len(child_ys) if child_ys else margin_top
            y_positions[id(node)] = y
            return x, y

        def get_x(node, depth=0):
            return margin_left + depth * (x_scale / max_depth)

        def get_depth(node, current=0, depths=None):
            if depths is None:
                depths = {}
            depths[id(node)] = current
            for child in node.children:
                get_depth(child, current + 1, depths)
            return depths

        depths = get_depth(self.current_tree)
        assign_positions(self.current_tree)

        def draw_node(node):
            d = depths[id(node)]
            x = margin_left + d * (x_scale / max_depth)
            y = y_positions[id(node)]

            for child in node.children:
                cd = depths[id(child)]
                cx = margin_left + cd * (x_scale / max_depth)
                cy = y_positions[id(child)]

                # Determine line color based on monophyly status
                color = COLORS["branch"]
                if not child.is_tip and child.monophyly_status:
                    if child.monophyly_status == "monophyletic":
                        color = COLORS["mono"]
                    elif child.monophyly_status == "not_monophyletic":
                        color = COLORS["non_mono"]

                line_width = 1.5 * zoom
                # Horizontal then vertical (rectangular phylogram)
                self.tree_canvas.create_line(
                    x, y, x, cy, fill=color, width=line_width)
                self.tree_canvas.create_line(
                    x, cy, cx, cy, fill=color, width=line_width)

                draw_node(child)

            # Draw labels
            if self.show_labels_var.get():
                font_size = max(7, int(9 * zoom))
                if node.is_tip:
                    self.tree_canvas.create_text(
                        x + 5, y, text=node.label, anchor="w",
                        fill=COLORS["tip"],
                        font=("Helvetica", font_size))
                elif node.label:
                    color = COLORS["unresolved"]
                    if node.monophyly_status == "monophyletic":
                        color = COLORS["mono"]
                    elif node.monophyly_status == "not_monophyletic":
                        color = COLORS["non_mono"]
                    self.tree_canvas.create_text(
                        x - 3, y - 8 * zoom, text=node.label, anchor="e",
                        fill=color,
                        font=("Helvetica", max(6, font_size - 1), "bold"))

        draw_node(self.current_tree)

        # Set scroll region
        total_height = margin_top + n_tips * tip_spacing + 50
        total_width = margin_left + x_scale + 300 * zoom
        self.tree_canvas.configure(
            scrollregion=(0, 0, total_width, total_height))

    def _redraw_tree(self):
        self._draw_tree()

    # ---- Analysis ----

    def _run_analysis(self):
        """Run monophyly analysis on the current tree."""
        if not self.db.loaded:
            messagebox.showwarning(
                "No Taxonomy", "Load taxonomy database first.")
            return
        if not self.current_tree:
            messagebox.showwarning("No Tree", "Load a tree file first.")
            return

        def analyze():
            checker = MonophylyChecker(self.db)
            self.results, self.unresolved_tips = checker.check_tree(
                self.current_tree, progress_callback=self._update_status
            )
            checker.label_tree(self.current_tree, self.results)

            # Store for batch results
            sel = self.tree_listbox.curselection()
            if sel:
                basename = os.path.basename(self.tree_files[sel[0]])
                self.batch_results[basename] = (
                    self.results, self.unresolved_tips)

            self.root.after(0, self._show_results)
            self.root.after(0, self._draw_tree)

        self._run_threaded(analyze, "Running monophyly analysis...")

    def _run_batch(self):
        """Run monophyly analysis on all loaded trees."""
        if not self.db.loaded:
            messagebox.showwarning(
                "No Taxonomy", "Load taxonomy database first.")
            return
        if not self.tree_files:
            messagebox.showwarning("No Trees", "Load tree files first.")
            return

        def analyze_batch():
            checker = MonophylyChecker(self.db)
            total = len(self.tree_files)
            for idx, fp in enumerate(self.tree_files):
                basename = os.path.basename(fp)
                self._update_status(
                    f"Analyzing {basename} ({idx + 1}/{total})...")
                try:
                    tree = parse_newick_file(fp)
                    results, unresolved = checker.check_tree(tree)
                    checker.label_tree(tree, results)
                    self.batch_results[basename] = (results, unresolved)
                except Exception as e:
                    self._update_status(f"Error in {basename}: {e}")

            self._update_status(
                f"Batch analysis complete: {len(self.batch_results)} trees")

            # Refresh current view
            sel = self.tree_listbox.curselection()
            if sel:
                self.root.after(0, lambda: self._load_tree(
                    self.tree_files[sel[0]]))

        self._run_threaded(analyze_batch, "Running batch analysis...")

    def _show_results(self):
        """Update the results panels."""
        self._clear_results()

        if not self.results:
            return

        checker = MonophylyChecker(self.db)
        summary = checker.get_summary(self.results)
        self.summary_text.insert(tk.END, summary)

        # Non-monophyletic table
        non_mono = [r for r in self.results if r.status == "not_monophyletic"]
        non_mono.sort(key=lambda x: len(x.tip_labels), reverse=True)
        for r in non_mono:
            intruder_str = ", ".join(r.intruders[:5])
            if len(r.intruders) > 5:
                intruder_str += f"... (+{len(r.intruders) - 5})"
            self.nonmono_tree.insert("", tk.END, values=(
                r.mrca_name, r.mrca_rank, len(r.tip_labels), intruder_str
            ))

        # Monophyletic table
        mono = [r for r in self.results if r.status == "monophyletic"]
        mono.sort(key=lambda x: len(x.tip_labels), reverse=True)
        for r in mono:
            self.mono_tree.insert("", tk.END, values=(
                r.mrca_name, r.mrca_rank, len(r.tip_labels)
            ))

        # Unresolved tips
        if self.unresolved_tips:
            self.unresolved_text.insert(
                tk.END,
                f"The following {len(self.unresolved_tips)} tip(s) could not "
                f"be matched to NCBI taxonomy:\n\n"
            )
            for tip in sorted(self.unresolved_tips):
                self.unresolved_text.insert(tk.END, f"  {tip}\n")

        # Switch to non-monophyletic tab if there are results
        # (PI said "if doesn't match is preferred over does match")
        if non_mono:
            self.results_notebook.select(1)

    def _clear_results(self):
        self.summary_text.delete("1.0", tk.END)
        for item in self.nonmono_tree.get_children():
            self.nonmono_tree.delete(item)
        for item in self.mono_tree.get_children():
            self.mono_tree.delete(item)
        self.unresolved_text.delete("1.0", tk.END)

    # ---- Export ----

    def _export_tree(self):
        """Export the labeled tree in Newick format."""
        if not self.current_tree:
            messagebox.showwarning("No Tree", "No tree loaded.")
            return

        filepath = filedialog.asksaveasfilename(
            title="Export Labeled Tree",
            defaultextension=".nwk",
            filetypes=[
                ("Newick files", "*.nwk *.tre"),
                ("All files", "*.*"),
            ]
        )
        if filepath:
            bl = self.show_branch_lengths_var.get()
            newick = self.current_tree.to_newick(branch_lengths=bl) + ";"
            with open(filepath, "w") as f:
                f.write(newick + "\n")
            self._update_status(f"Exported tree to {filepath}")

    def _export_report(self):
        """Export the analysis report."""
        if not self.results:
            messagebox.showwarning("No Results", "Run analysis first.")
            return

        filepath = filedialog.asksaveasfilename(
            title="Export Report",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ]
        )
        if filepath:
            checker = MonophylyChecker(self.db)
            report = checker.get_summary(self.results)

            if self.unresolved_tips:
                report += "\n\n--- Unresolved Tips ---\n"
                for tip in sorted(self.unresolved_tips):
                    report += f"  {tip}\n"

            with open(filepath, "w") as f:
                f.write(report + "\n")
            self._update_status(f"Exported report to {filepath}")

    # ---- About ----

    def _show_about(self):
        messagebox.showinfo(
            "About PhyLabeler",
            "PhyLabeler - Phylogenetic Tree Labeler\n\n"
            "Checks gene trees against NCBI taxonomy\n"
            "and identifies monophyletic/non-monophyletic clades.\n\n"
            "Based on LabelPhy by J.F. Walker\n"
            "https://github.com/jfwalker/LabelPhy\n"
            "PeerJ CS 56 (2016)\n\n"
            "No external dependencies required.\n"
            "Uses Python standard library only."
        )

    # ---- Helpers ----

    def _run_threaded(self, target, status_msg="Working..."):
        """Run a function in a background thread."""
        self._update_status(status_msg)

        def wrapper():
            try:
                target()
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror(
                    "Error", str(e)))
            finally:
                self.root.after(0, lambda: self._update_status("Ready"))

        thread = threading.Thread(target=wrapper, daemon=True)
        thread.start()

    def _update_status(self, msg):
        """Thread-safe status bar update."""
        self.root.after(0, lambda: self.status_var.set(msg))


def main():
    root = tk.Tk()
    app = PhyLabelerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
