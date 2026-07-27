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
import shutil
import tempfile
import threading
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from tkinter import font as tkfont

from archive_io import extract_matching, is_archive
from taxonomy_db import TaxonomyDB
from tree_parser import is_tree_filename, parse_newick_file
from monophyly import (MAIN_RANKS, CollapsedGroup, ConflictingTip,
                       MonophylyChecker, plural as _plural)
from help_content import HELP_SECTIONS, PANEL_HELP

# File-dialog filters, kept next to each other so CLI and GUI stay in sync.
TREE_FILETYPES = [
    ("Newick files", "*.nwk *.tre *.tree *.newick"),
    ("Compressed Newick", "*.nwk.gz *.tre.gz *.tree.gz *.newick.gz"),
    ("Tree archives", "*.tar.gz *.tgz *.tar *.zip"),
    ("All files", "*.*"),
]

TAXONOMY_FILETYPES = [
    ("Taxonomy dump or archive",
     "*.dmp *.dmp.gz *.tar.gz *.tgz *.tar *.zip"),
    ("DMP files", "*.dmp *.dmp.gz"),
    ("NCBI archives", "*.tar.gz *.tgz *.tar *.zip"),
    ("All files", "*.*"),
]


# ---- Tree layout ----
TIP_SPACING = 16          # px between adjacent tips at zoom 1.0
LEVEL_STEP = 26           # px per depth level in cladogram mode
LEVEL_STEP_MIN = 10       # never squeeze levels below this
PHYLOGRAM_WIDTH = 620     # px the longest root-to-tip path spans
MARGIN_LEFT = 46          # room for support values left of the root
MARGIN_TOP = 20
LABEL_GUTTER = 320        # room for tip names right of the last node


# ---- Text sizing ----
# A discrete ladder rather than free scaling, so each click is a visible step
# and the ends are hard limits.
HELP_FONT_SIZES = (9, 10, 11, 12, 14, 16, 18, 20, 24, 28)
HELP_DEFAULT_FONT_SIZE = 11

RESULTS_FONT_SIZES = (9, 10, 11, 12, 14, 16, 18, 20)
RESULTS_DEFAULT_FONT_SIZE = 10

# Style name for the results tables, so scaling them cannot disturb any other
# ttk.Treeview the app might grow later.
RESULTS_TREE_STYLE = "Results.Treeview"


# ---- Color scheme ----
COLORS = {
    "mono": "#2d8a4e",        # green for monophyletic
    "non_mono": "#d63031",    # red for non-monophyletic
    "unresolved": "#6c757d",  # gray for unresolved
    "inexact": "#b8860b",     # amber for tips matched below species level
    "support": "#8a8a8a",     # bootstrap/posterior values
    "clade_band": "#fff3bf",  # band behind a selected clade
    "chip_edge": "#dddddd",   # border of the clade-name chip
    "tip": "#2c3e50",         # dark for tips
    "branch": "#555555",      # branch lines
    "bg": "#ffffff",          # canvas background
    "highlight": "#ffeaa7",   # selection highlight
    "header_bg": "#34495e",   # header background
    "header_fg": "#ecf0f1",   # header text
    "app_bg": "#2f2f2f",      # window background
    "panel_bg": "#3a3a3a",    # frame background
    "panel_edge": "#565656",  # frame border
    "text": "#f1f1f1",        # primary text
    "muted_text": "#d0d0d0",  # secondary text
    "button_bg": "#f2f2f2",   # light buttons for contrast
    "button_fg": "#1f1f1f",   # readable button text
    "tab_bg": "#4a4a4a",      # notebook tab background
    "tab_active": "#f2f2f2",  # selected tab background
    "tab_active_fg": "#202020",
    "field_bg": "#1f1f1f",    # list/text background
    "field_fg": "#f5f5f5",
    "info_icon": "#7fb3ff",   # the (i) affordance
    "tooltip_bg": "#f7f7f2",  # hover popups read as paper, not as chrome
    "tooltip_fg": "#1f1f1f",
}


def mousewheel_units(event):
    """Normalize a mouse wheel event to scroll steps (positive = down)."""
    if getattr(event, "delta", 0):
        if abs(event.delta) >= 120:
            return -int(event.delta / 120)
        return -1 if event.delta > 0 else 1
    if getattr(event, "num", None) == 4:
        return -1
    if getattr(event, "num", None) == 5:
        return 1
    return 0


class TextScaler:
    """An A- / A+ / Reset control for a group of widgets.

    Shared by the help window and the results panel so the two behave
    identically instead of drifting apart. Sizes step along a ladder, which
    makes every click a visible change and gives the ends hard limits.

    The chosen size lives here rather than on the widgets, so it survives a
    window being closed and rebuilt.
    """

    def __init__(self, sizes, default, family="Courier"):
        self.sizes = sizes
        self.default = default
        self.family = family
        self.size = default
        self._texts = []
        self._tree_style = None
        self._label = None
        self._smaller = None
        self._bigger = None

    # ---- registration ----

    def add_text(self, widget):
        """Scale a Text/ScrolledText widget."""
        self._texts.append(widget)
        return widget

    def scale_treeview(self, style_name):
        """Also scale the rows and headings of tables using this ttk style."""
        self._tree_style = style_name

    def release_widgets(self):
        """Forget widgets whose window has been destroyed; keep the size."""
        self._texts = []
        self._label = self._smaller = self._bigger = None

    # ---- controls ----

    def build_controls(self, parent, hint=None, label="Text size:"):
        """Build the A- / A+ / Reset row. Returns the containing frame."""
        bar = ttk.Frame(parent)
        ttk.Label(bar, text=label).pack(side=tk.LEFT)
        self._smaller = ttk.Button(bar, text="A−", width=3,
                                   style="App.TButton",
                                   command=lambda: self.zoom(-1))
        self._smaller.pack(side=tk.LEFT, padx=(6, 2))
        self._bigger = ttk.Button(bar, text="A+", width=3,
                                  style="App.TButton",
                                  command=lambda: self.zoom(1))
        self._bigger.pack(side=tk.LEFT, padx=2)
        ttk.Button(bar, text="Reset", style="App.TButton",
                   command=self.reset).pack(side=tk.LEFT, padx=(6, 2))
        self._label = ttk.Label(bar, text="")
        self._label.pack(side=tk.LEFT, padx=10)
        if hint:
            ttk.Label(bar, text=hint,
                      foreground=COLORS["muted_text"]).pack(side=tk.RIGHT)
        return bar

    def bind_keys(self, toplevel):
        """Ctrl/Cmd with plus, minus and zero."""
        for seq in ("<Control-plus>", "<Control-equal>", "<Control-KP_Add>",
                    "<Command-plus>", "<Command-equal>"):
            toplevel.bind(seq, lambda e: self.zoom(1))
        for seq in ("<Control-minus>", "<Control-KP_Subtract>",
                    "<Command-minus>"):
            toplevel.bind(seq, lambda e: self.zoom(-1))
        for seq in ("<Control-0>", "<Command-0>"):
            toplevel.bind(seq, lambda e: self.reset())

    def bind_wheel(self, widget):
        """Ctrl/Cmd + wheel resizes; the bare wheel still scrolls."""
        for seq in ("<Control-MouseWheel>", "<Command-MouseWheel>"):
            widget.bind(seq, self._on_wheel, add="+")

    def _on_wheel(self, event):
        step = mousewheel_units(event)
        if step:
            self.zoom(-step)    # wheel up is negative, and should enlarge
        return "break"

    # ---- actions ----

    def zoom(self, step):
        try:
            index = self.sizes.index(self.size)
        except ValueError:
            index = min(range(len(self.sizes)),
                        key=lambda i: abs(self.sizes[i] - self.size))
        self.size = self.sizes[max(0, min(len(self.sizes) - 1, index + step))]
        self.apply()

    def reset(self):
        self.size = self.default
        self.apply()

    def apply(self):
        """Push the current size to every registered widget and control."""
        for text in self._texts:
            text.configure(font=(self.family, self.size))

        if self._tree_style:
            metrics = tkfont.Font(family=self.family, size=self.size)
            style = ttk.Style()
            style.configure(self._tree_style, font=(self.family, self.size),
                            rowheight=metrics.metrics("linespace") + 10)
            style.configure(f"{self._tree_style}.Heading",
                            font=(self.family, self.size))

        if self._label is not None:
            self._label.configure(text=f"{self.size} pt")
        if self._smaller is not None:
            self._smaller.configure(
                state=tk.DISABLED if self.size <= self.sizes[0] else tk.NORMAL)
            self._bigger.configure(
                state=tk.DISABLED if self.size >= self.sizes[-1] else tk.NORMAL)


class Tooltip:
    """Hover text for a widget.

    Plain Tk: a borderless Toplevel shown after a short delay, so brushing
    past a control does not flash a popup.
    """

    DELAY_MS = 350

    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.window = None
        self._job = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._job = self.widget.after(self.DELAY_MS, self._show)

    def _cancel(self):
        if self._job is not None:
            self.widget.after_cancel(self._job)
            self._job = None

    def _show(self):
        if self.window is not None:
            return
        x = self.widget.winfo_rootx()
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self.window = tk.Toplevel(self.widget)
        self.window.wm_overrideredirect(True)
        self.window.wm_geometry(f"+{x}+{y}")
        tk.Label(
            self.window, text=self.text, justify=tk.LEFT,
            background=COLORS["tooltip_bg"], foreground=COLORS["tooltip_fg"],
            relief=tk.SOLID, borderwidth=1, padx=8, pady=6,
            font=("Helvetica", 11),
        ).pack()

    def _hide(self, _event=None):
        self._cancel()
        if self.window is not None:
            self.window.destroy()
            self.window = None


class PhyLabelerApp:
    """Main application window."""

    def __init__(self, root):
        self.root = root
        self.root.title("PhyLabeler - Phylogenetic Tree Labeler")
        self.root.geometry("1200x800")
        self.root.minsize(900, 600)
        self._configure_styles()

        # State
        self.db = TaxonomyDB()
        self.current_tree = None
        self.results = None
        self.unresolved_tips = []
        self.tree_files = []        # For batch mode
        self.batch_results = {}     # filename -> analysis dict
        self.taxon_results = []     # One entry per named taxonomic group
        self.groups = []            # Taxa merged when they cover the same tips
        self.conflicting_tips = []  # Tips ranked by groups they would repair
        self.tip_matches = {}       # tip label -> TipMatch
        self.highlight_tips = set()      # Tips of the selected clade
        self.highlight_intruders = set()  # Tips that break the selected group
        self._temp_dirs = []        # Scratch dirs for archive extraction
        self._redraw_job = None     # Pending debounced canvas redraw
        self._busy = False          # An analysis is running
        # Keyed by (widget name, row id): Treeview row ids are only unique
        # within one widget, so every table hands out I001, I002, ... and a
        # flat dict would have the tables silently overwriting each other.
        self._row_results = {}

        self._help_window = None
        self.help_scaler = TextScaler(HELP_FONT_SIZES, HELP_DEFAULT_FONT_SIZE)
        self.results_scaler = TextScaler(RESULTS_FONT_SIZES,
                                         RESULTS_DEFAULT_FONT_SIZE)
        self.results_scaler.scale_treeview(RESULTS_TREE_STYLE)
        self._build_ui()
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
        # Size the split once Tk knows how wide the controls actually want to
        # be. Left at the default the third taxonomy button was sliced down to
        # a blank sliver, which reads as a broken control rather than as a
        # pane that needs dragging.
        self.root.after_idle(self._fit_left_pane)
        self._check_cache()

    def _fit_left_pane(self):
        """Widen the split so no control in the left panel is cut off."""
        try:
            self.left_frame.update_idletasks()
            needed = self.left_frame.winfo_reqwidth() + 16
            if self.main_pane.sashpos(0) < needed:
                self.main_pane.sashpos(0, needed)
        except tk.TclError:
            pass    # window torn down before it was ever mapped

    def _on_close(self):
        """Remove any archive scratch directories, then close the window."""
        for path in self._temp_dirs:
            shutil.rmtree(path, ignore_errors=True)
        self._temp_dirs.clear()
        self.root.destroy()

    def _configure_styles(self):
        """Configure a consistent ttk theme so controls render correctly."""
        self.root.configure(bg=COLORS["app_bg"])

        style = ttk.Style(self.root)
        if "clam" in style.theme_names():
            style.theme_use("clam")

        style.configure(
            ".",
            background=COLORS["app_bg"],
            foreground=COLORS["text"],
            fieldbackground=COLORS["field_bg"],
        )
        style.configure("TFrame", background=COLORS["app_bg"])
        style.configure(
            "TLabelframe",
            background=COLORS["panel_bg"],
            bordercolor=COLORS["panel_edge"],
            relief="solid",
        )
        style.configure(
            "TLabelframe.Label",
            background=COLORS["app_bg"],
            foreground=COLORS["text"],
        )
        style.configure("TLabel", background=COLORS["app_bg"],
                        foreground=COLORS["text"])
        style.configure(
            "TPanedwindow",
            background=COLORS["app_bg"],
            sashthickness=6,
        )
        style.configure(
            "App.TButton",
            background=COLORS["button_bg"],
            foreground=COLORS["button_fg"],
            bordercolor=COLORS["panel_edge"],
            focuscolor="none",
            padding=(10, 6),
        )
        style.map(
            "App.TButton",
            background=[
                ("active", "#ffffff"),
                ("pressed", "#dcdcdc"),
                ("disabled", "#c8c8c8"),
            ],
            foreground=[
                ("disabled", "#666666"),
            ],
        )
        style.configure(
            "TCheckbutton",
            background=COLORS["app_bg"],
            foreground=COLORS["text"],
        )
        style.map(
            "TCheckbutton",
            background=[("active", COLORS["app_bg"])],
            foreground=[("disabled", COLORS["muted_text"])],
        )
        style.configure(
            "Horizontal.TScale",
            background=COLORS["app_bg"],
            troughcolor=COLORS["panel_edge"],
        )
        style.configure(
            "TNotebook",
            background=COLORS["panel_bg"],
            bordercolor=COLORS["panel_edge"],
            tabmargins=(4, 4, 4, 0),
        )
        style.configure(
            "TNotebook.Tab",
            background=COLORS["tab_bg"],
            foreground=COLORS["text"],
            padding=(12, 6),
            lightcolor=COLORS["panel_edge"],
            bordercolor=COLORS["panel_edge"],
        )
        style.map(
            "TNotebook.Tab",
            background=[
                ("selected", COLORS["tab_active"]),
                ("active", "#5c5c5c"),
            ],
            foreground=[
                ("selected", COLORS["tab_active_fg"]),
                ("active", COLORS["text"]),
            ],
        )
        style.configure(
            "Treeview",
            background=COLORS["field_bg"],
            foreground=COLORS["field_fg"],
            fieldbackground=COLORS["field_bg"],
            bordercolor=COLORS["panel_edge"],
        )
        style.configure(
            "Treeview.Heading",
            background=COLORS["panel_bg"],
            foreground=COLORS["text"],
            relief="flat",
        )
        style.map(
            "Treeview",
            background=[("selected", COLORS["highlight"])],
            foreground=[("selected", COLORS["button_fg"])],
        )
        # The results tables carry their own style so the text-size control
        # can scale them without touching any other table.
        style.configure(RESULTS_TREE_STYLE,
                        background=COLORS["field_bg"],
                        foreground=COLORS["field_fg"],
                        fieldbackground=COLORS["field_bg"],
                        bordercolor=COLORS["panel_edge"])
        style.map(
            RESULTS_TREE_STYLE,
            background=[("selected", COLORS["highlight"])],
            foreground=[("selected", COLORS["button_fg"])],
        )
        style.configure(f"{RESULTS_TREE_STYLE}.Heading",
                        background=COLORS["panel_bg"],
                        foreground=COLORS["text"],
                        relief="flat")

        # Read-only comboboxes render their text through the *selection*
        # colours in clam, which left pale grey on pale grey - the filters
        # were effectively unreadable. Style every state explicitly.
        style.configure(
            "TCombobox",
            fieldbackground=COLORS["field_bg"],
            background=COLORS["panel_bg"],
            foreground=COLORS["field_fg"],
            arrowcolor=COLORS["text"],
            bordercolor=COLORS["panel_edge"],
            lightcolor=COLORS["panel_edge"],
            darkcolor=COLORS["panel_edge"],
            selectbackground=COLORS["field_bg"],
            selectforeground=COLORS["field_fg"],
            padding=(6, 3),
        )
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", COLORS["field_bg"]),
                             ("disabled", COLORS["panel_bg"])],
            foreground=[("readonly", COLORS["field_fg"]),
                        ("disabled", COLORS["muted_text"])],
            selectbackground=[("readonly", COLORS["field_bg"])],
            selectforeground=[("readonly", COLORS["field_fg"])],
            background=[("readonly", COLORS["panel_bg"]),
                        ("active", COLORS["panel_bg"])],
            arrowcolor=[("readonly", COLORS["text"]),
                        ("disabled", COLORS["muted_text"])],
        )
        # The drop-down itself is a plain Tk listbox, reachable only through
        # the option database.
        self.root.option_add("*TCombobox*Listbox.background",
                             COLORS["field_bg"])
        self.root.option_add("*TCombobox*Listbox.foreground",
                             COLORS["field_fg"])
        self.root.option_add("*TCombobox*Listbox.selectBackground",
                             COLORS["highlight"])
        self.root.option_add("*TCombobox*Listbox.selectForeground",
                             COLORS["button_fg"])

        style.configure(
            "TScrollbar",
            background=COLORS["panel_bg"],
            troughcolor=COLORS["app_bg"],
            bordercolor=COLORS["panel_edge"],
            arrowcolor=COLORS["text"],
        )

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
        file_menu.add_command(label="Open Trees from Archive...",
                              command=self._open_tree_archive)
        file_menu.add_separator()
        file_menu.add_command(label="Export Labeled Tree...",
                              command=self._export_tree)
        file_menu.add_command(label="Export Report...",
                              command=self._export_report)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self._on_close,
                              accelerator="Ctrl+Q")
        menubar.add_cascade(label="File", menu=file_menu)

        db_menu = tk.Menu(menubar, tearoff=0)
        db_menu.add_command(label="Load from Files...",
                            command=self._load_taxonomy_files)
        db_menu.add_command(label="Load Archive (taxdump.tar.gz)...",
                            command=self._open_taxonomy_archive)
        db_menu.add_command(label="Download from NCBI...",
                            command=self._download_taxonomy)
        db_menu.add_separator()
        db_menu.add_command(label="Cache Info...", command=self._show_cache_info)
        menubar.add_cascade(label="Taxonomy DB", menu=db_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        for title, _body in HELP_SECTIONS:
            help_menu.add_command(
                label=title, command=lambda t=title: self._show_help(t))
        help_menu.add_command(
            label="Project Guide",
            command=lambda: self._show_help("Project Guide"))
        help_menu.add_separator()
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.bind("<Control-o>", lambda e: self._open_tree())
        self.root.bind("<Control-q>", lambda e: self._on_close())
        self.root.bind("<F1>", lambda e: self._show_help())

        # Main layout: left panel + right panel
        main_pane = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.main_pane = main_pane

        # ---- Left panel: controls + file list ----
        left_frame = ttk.Frame(main_pane, width=320)
        main_pane.add(left_frame, weight=1)
        self.left_frame = left_frame

        # Taxonomy DB status
        db_frame = self._panel(left_frame, "Taxonomy Database", "taxonomy")
        db_frame.pack(fill=tk.X, padx=5, pady=5)

        self.db_status_var = tk.StringVar(value="Not loaded")
        ttk.Label(db_frame, textvariable=self.db_status_var).pack(
            padx=5, pady=2, anchor="w")

        db_btn_frame = ttk.Frame(db_frame)
        db_btn_frame.pack(fill=tk.X, padx=5, pady=2)
        for column, (label, command) in enumerate((
                ("Load Files", self._load_taxonomy_files),
                ("Download NCBI", self._download_taxonomy),
                ("Load Cache", self._load_cache))):
            ttk.Button(db_btn_frame, text=label, command=command,
                       style="App.TButton").grid(
                           row=0, column=column, sticky="ew", padx=2)
            db_btn_frame.grid_columnconfigure(column, weight=1)

        # Tree files list
        tree_frame = self._panel(left_frame, "Tree Files", "trees")
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        tree_btn_frame = ttk.Frame(tree_frame)
        tree_btn_frame.pack(fill=tk.X, padx=5, pady=2)
        for column, (label, command) in enumerate((
                ("Add Tree(s)", self._open_tree),
                ("Clear", self._clear_trees))):
            ttk.Button(tree_btn_frame, text=label, command=command,
                       style="App.TButton").grid(
                           row=0, column=column, sticky="ew", padx=2)
            tree_btn_frame.grid_columnconfigure(column, weight=1)

        self.tree_listbox = tk.Listbox(tree_frame, selectmode=tk.BROWSE)
        self.tree_listbox.pack(fill=tk.BOTH, expand=True, padx=5, pady=2)
        self.tree_listbox.bind("<<ListboxSelect>>", self._on_tree_select)
        self.tree_listbox.configure(
            bg=COLORS["field_bg"],
            fg=COLORS["field_fg"],
            selectbackground=COLORS["highlight"],
            selectforeground=COLORS["button_fg"],
            highlightbackground=COLORS["panel_edge"],
            highlightcolor=COLORS["panel_edge"],
        )

        # Analysis controls
        analysis_frame = self._panel(left_frame, "Analysis", "analysis")
        analysis_frame.pack(fill=tk.X, padx=5, pady=5)

        self.run_button = ttk.Button(
            analysis_frame, text="Run Monophyly Check",
            command=self._run_analysis, style="App.TButton")
        self.run_button.pack(fill=tk.X, padx=5, pady=2)
        self.batch_button = ttk.Button(
            analysis_frame, text="Run Batch Analysis",
            command=self._run_batch, style="App.TButton")
        self.batch_button.pack(fill=tk.X, padx=5, pady=2)

        rank_row = ttk.Frame(analysis_frame)
        rank_row.pack(fill=tk.X, padx=5, pady=(6, 2))
        ttk.Label(rank_row, text="Show rank:").pack(side=tk.LEFT)
        self.rank_filter_var = tk.StringVar(value="all ranks")
        # "clade" and "no rank" carry most of the deep conflicts (rosids,
        # fabids, asterids, ...), so leaving them out of the filter made the
        # worst findings in the tree unreachable.
        self.rank_filter = ttk.Combobox(
            rank_row, textvariable=self.rank_filter_var, state="readonly",
            values=["all ranks"] + list(MAIN_RANKS) + ["clade", "no rank"],
            width=11)
        self.rank_filter.pack(side=tk.LEFT, padx=4)
        self.rank_filter.bind(
            "<<ComboboxSelected>>", lambda e: self._populate_group_tables())

        conflict_row = ttk.Frame(analysis_frame)
        conflict_row.pack(fill=tk.X, padx=5, pady=(0, 4))
        ttk.Label(conflict_row, text="Conflict:").pack(side=tk.LEFT)
        self.conflict_filter_var = tk.StringVar(value="all conflicts")
        conflict_filter = ttk.Combobox(
            conflict_row, textvariable=self.conflict_filter_var,
            state="readonly", width=11,
            values=["all conflicts", "intruded", "scattered", "both"])
        conflict_filter.pack(side=tk.LEFT, padx=4)
        conflict_filter.bind(
            "<<ComboboxSelected>>", lambda e: self._populate_group_tables())

        self.collapse_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            analysis_frame, text="Merge groups with identical tips",
            variable=self.collapse_var,
            command=self._populate_group_tables).pack(
                anchor="w", padx=5, pady=(0, 4))

        # Tip search
        search_frame = self._panel(left_frame, "Find Tip", "search")
        search_frame.pack(fill=tk.X, padx=5, pady=5)
        self.search_var = tk.StringVar()
        search_entry = ttk.Entry(search_frame, textvariable=self.search_var)
        search_entry.pack(fill=tk.X, padx=5, pady=4)
        search_entry.bind("<Return>", lambda e: self._search_tips())
        self.search_var.trace_add("write", lambda *a: self._search_tips())

        # ---- Right panel: tree view + results ----
        right_frame = ttk.Frame(main_pane)
        main_pane.add(right_frame, weight=3)

        right_pane = ttk.PanedWindow(right_frame, orient=tk.VERTICAL)
        right_pane.pack(fill=tk.BOTH, expand=True)

        # Tree visualization
        tree_view_frame = self._panel(right_pane, "Tree Visualization",
                                     "tree_view")
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

        self.show_branch_lengths_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(zoom_frame, text="Branch lengths",
                         variable=self.show_branch_lengths_var,
                         command=self._redraw_tree).pack(side=tk.RIGHT)

        self.show_support_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(zoom_frame, text="Support",
                         variable=self.show_support_var,
                         command=self._redraw_tree).pack(side=tk.RIGHT)

        self.show_labels_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(zoom_frame, text="Labels",
                         variable=self.show_labels_var,
                         command=self._redraw_tree).pack(side=tk.RIGHT)

        # Results panel
        results_frame = self._panel(right_pane, "Results", "results")
        right_pane.add(results_frame, weight=2)

        # Text size applies to the whole results panel - the two text tabs
        # and the row/heading font of the three tables - so switching tabs
        # never lands on a control that does nothing.
        self.results_scaler.build_controls(
            results_frame, hint="(Ctrl/Cmd +, −, 0)").pack(
                fill=tk.X, padx=5, pady=(4, 0))

        # Notebook for different result views
        self.results_notebook = ttk.Notebook(results_frame)
        self.results_notebook.pack(fill=tk.BOTH, expand=True, padx=5,
                                   pady=(2, 5))

        # Summary tab
        summary_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(summary_frame, text="Summary")
        self.summary_text = scrolledtext.ScrolledText(
            summary_frame, wrap=tk.WORD, font=("Courier", 10))
        self.summary_text.pack(fill=tk.BOTH, expand=True)
        self._style_text_widget(self.summary_text)
        self.results_scaler.add_text(self.summary_text)
        self.results_scaler.bind_wheel(self.summary_text)

        # Broken-groups tab: one row per taxonomic group, not per tree node,
        # and taxa covering the same tips share a single expandable row.
        self.nonmono_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.nonmono_frame, text="Broken Groups")
        self.nonmono_tree = self._build_result_table(
            self.nonmono_frame,
            columns=[("name", "Group", 200), ("rank", "Rank", 120),
                     ("members", "Members", 70),
                     ("cost", "Tips to move", 90), ("type", "Conflict", 80),
                     ("culprit", "One tip explains it", 175),
                     ("intruders", "Intruders", 70), ("outliers", "Outliers", 70),
                     ("tips", "Intruding tips", 380)],
            expandable=True,
        )

        # Monophyletic tab
        self.mono_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.mono_frame, text="Monophyletic")
        self.mono_tree = self._build_result_table(
            self.mono_frame,
            columns=[("name", "Group", 260), ("rank", "Rank", 130),
                     ("members", "Members", 80), ("note", "Notes", 240)],
            expandable=True,
        )

        # Conflicting tips: which sequences account for the most breakage.
        self.tips_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.tips_frame, text="Conflicting Tips")
        self.conflict_tree = self._build_result_table(
            self.tips_frame,
            columns=[("tip", "Tip", 240), ("repairs", "Groups repaired", 120),
                     ("role", "Role", 90),
                     ("appears", "Intrudes into", 100),
                     ("groups", "Groups it would repair", 480)],
        )
        self.conflict_banner = ttk.Label(self.tips_frame, text="")
        self.conflict_banner.grid(row=2, column=0, columnspan=2,
                                  sticky="w", padx=6, pady=(2, 4))

        # Tip matching tab
        unresolved_frame = ttk.Frame(self.results_notebook)
        self.results_notebook.add(unresolved_frame, text="Tip Matching")
        self.unresolved_text = scrolledtext.ScrolledText(
            unresolved_frame, wrap=tk.WORD, font=("Courier", 10))
        self.unresolved_text.pack(fill=tk.BOTH, expand=True)
        self._style_text_widget(self.unresolved_text)
        self.results_scaler.add_text(self.unresolved_text)
        self.results_scaler.bind_wheel(self.unresolved_text)

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(self.root, textvariable=self.status_var,
                                relief=tk.SUNKEN, anchor="w")
        status_bar.pack(fill=tk.X, side=tk.BOTTOM, padx=5)

        self._bind_scroll_wheels()
        self.results_scaler.bind_keys(self.root)
        self.results_scaler.apply()

    def _panel(self, parent, title, help_key):
        """A LabelFrame whose title carries a hoverable (i).

        The icon is built as the frame's `labelwidget`, so the marker sits in
        the panel's own title rather than stealing a row of the layout.
        """
        header = ttk.Frame(parent)
        ttk.Label(header, text=title).pack(side=tk.LEFT)

        tooltip_text, section = PANEL_HELP[help_key]
        icon = ttk.Label(header, text=" ⓘ", foreground=COLORS["info_icon"],
                         cursor="hand2")
        icon.pack(side=tk.LEFT)
        Tooltip(icon, f"{tooltip_text}\n\nClick for more.")
        icon.bind("<Button-1>", lambda e, s=section: self._show_help(s))

        return ttk.LabelFrame(parent, labelwidget=header)

    # ---- Help ----

    def _show_help(self, section=None):
        """Open the help window, optionally at a named section."""
        if getattr(self, "_help_window", None) is not None:
            try:
                self._help_window.deiconify()
                self._help_window.lift()
                if section:
                    self._select_help_section(section)
                return
            except tk.TclError:
                self._help_window = None

        window = tk.Toplevel(self.root)
        self._help_window = window
        window.title("PhyLabeler Help")
        window.geometry("880x680")
        window.minsize(620, 460)
        window.configure(bg=COLORS["app_bg"])
        window.protocol("WM_DELETE_WINDOW", self._close_help)

        self.help_scaler.release_widgets()
        self.help_scaler.build_controls(
            window, hint="(Ctrl/Cmd +, −, 0, or Ctrl+scroll)").pack(
                fill=tk.X, padx=8, pady=(8, 0))

        notebook = ttk.Notebook(window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=8, pady=(4, 8))
        self._help_notebook = notebook
        self._help_tabs = {}

        pages = list(HELP_SECTIONS) + [
            ("Project Guide", self._load_project_guide())]
        for title, body in pages:
            frame = ttk.Frame(notebook)
            notebook.add(frame, text=title)
            self._help_tabs[title] = frame
            text = scrolledtext.ScrolledText(frame, wrap=tk.WORD)
            text.pack(fill=tk.BOTH, expand=True)
            self._style_text_widget(text)
            self._bind_mousewheel(text, yview=True)
            self.help_scaler.bind_wheel(text)
            text.insert(tk.END, body)
            text.configure(state=tk.DISABLED)
            self.help_scaler.add_text(text)

        self.help_scaler.bind_keys(window)
        self.help_scaler.apply()

        if section:
            self._select_help_section(section)

    def _select_help_section(self, section):
        frame = self._help_tabs.get(section)
        if frame is not None:
            self._help_notebook.select(frame)

    def _close_help(self):
        if getattr(self, "_help_window", None) is not None:
            self._help_window.destroy()
            self._help_window = None
        # The chosen size lives on the scaler, so reopening honours it.
        self.help_scaler.release_widgets()

    def _build_result_table(self, parent, columns, expandable=False):
        """Create a results table with both scrollbars and click-to-navigate.

        The old tables had no horizontal scrollbar, so the intruders column -
        the one that says *why* a group is broken - was simply cut off at the
        window edge with no way to read it.

        `expandable` keeps the tree column visible so rows that merge several
        taxa can be opened to show what they merged.
        """
        table = ttk.Treeview(
            parent, columns=[c[0] for c in columns],
            show="tree headings" if expandable else "headings",
            selectmode="browse", style=RESULTS_TREE_STYLE)
        if expandable:
            table.column("#0", width=26, minwidth=26, stretch=False)
            table.heading("#0", text="")
        for key, heading, width in columns:
            table.heading(key, text=heading)
            table.column(key, width=width, minwidth=50, stretch=False)

        v_scroll = ttk.Scrollbar(parent, orient=tk.VERTICAL,
                                 command=table.yview)
        h_scroll = ttk.Scrollbar(parent, orient=tk.HORIZONTAL,
                                 command=table.xview)
        table.configure(yscrollcommand=v_scroll.set,
                        xscrollcommand=h_scroll.set)

        table.grid(row=0, column=0, sticky="nsew")
        v_scroll.grid(row=0, column=1, sticky="ns")
        h_scroll.grid(row=1, column=0, sticky="ew")
        parent.grid_rowconfigure(0, weight=1)
        parent.grid_columnconfigure(0, weight=1)

        table.bind("<<TreeviewSelect>>", self._on_result_select)
        return table

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
        """Load taxonomy from local dmp files or an un-extracted NCBI dump.

        Picking taxdump.tar.gz / taxdump.tgz / taxdmp.zip in the first dialog
        skips the second one, since the archive already holds both files.
        """
        names_file = filedialog.askopenfilename(
            title="Select names.dmp or an NCBI taxonomy archive",
            filetypes=TAXONOMY_FILETYPES
        )
        if not names_file:
            return

        if is_archive(names_file):
            self._load_taxonomy_archive(names_file)
            return

        nodes_file = filedialog.askopenfilename(
            title="Select nodes.dmp",
            filetypes=TAXONOMY_FILETYPES
        )
        if not nodes_file:
            return

        def load():
            self.db.load(names_file=names_file, nodes_file=nodes_file,
                         progress_callback=self._update_status)
            self.root.after(0, lambda: self.db_status_var.set(
                f"Loaded: {len(self.db.code_to_parent):,} taxa"))

        self._run_threaded(load, "Loading taxonomy files...")

    def _open_taxonomy_archive(self):
        """Pick an un-extracted NCBI taxonomy archive and load it."""
        archive = filedialog.askopenfilename(
            title="Select taxdump.tar.gz / taxdmp.zip",
            filetypes=[
                ("NCBI archives", "*.tar.gz *.tgz *.tar *.zip"),
                ("All files", "*.*"),
            ]
        )
        if archive:
            self._load_taxonomy_archive(archive)

    def _load_taxonomy_archive(self, archive):
        """Parse names.dmp and nodes.dmp out of an archive in the background."""
        def load():
            self.db.load(archive=archive,
                         progress_callback=self._update_status)
            self.root.after(0, lambda: self.db_status_var.set(
                f"Loaded: {len(self.db.code_to_parent):,} taxa"))

        self._run_threaded(
            load, f"Reading {os.path.basename(archive)}...")

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
        """Open one or more tree files, or an archive of them.

        Selecting a .tar.gz/.tgz/.zip unpacks its Newick members into a
        temporary directory and adds them all to the list.
        """
        filepaths = filedialog.askopenfilenames(
            title="Select Tree File(s) or Archive",
            filetypes=TREE_FILETYPES
        )
        if not filepaths:
            return

        added = []
        for fp in filepaths:
            if is_archive(fp):
                added.extend(self._extract_tree_archive(fp))
            else:
                added.append(fp)

        self._add_tree_files(added, select_last=True)

    def _open_batch(self):
        """Open a directory of tree files."""
        dirpath = filedialog.askdirectory(title="Select Directory with Trees")
        if not dirpath:
            return
        self._add_tree_files(
            [os.path.join(dirpath, name)
             for name in sorted(os.listdir(dirpath))
             if is_tree_filename(name)]
        )

    def _open_tree_archive(self):
        """Open an archive of tree files (.tar.gz, .tgz, .tar, .zip)."""
        archive = filedialog.askopenfilename(
            title="Select Archive of Trees",
            filetypes=[
                ("Tree archives", "*.tar.gz *.tgz *.tar *.zip"),
                ("All files", "*.*"),
            ]
        )
        if not archive:
            return

        extracted = self._extract_tree_archive(archive)
        if not extracted:
            messagebox.showwarning(
                "No Trees",
                f"No Newick files found in {os.path.basename(archive)}."
            )
            return
        self._add_tree_files(extracted, select_last=True)

    def _extract_tree_archive(self, archive):
        """Unpack the Newick members of an archive into a temp directory."""
        workdir = tempfile.mkdtemp(prefix="phylabeler_")
        self._temp_dirs.append(workdir)
        try:
            return extract_matching(archive, workdir, is_tree_filename)
        except Exception as e:
            messagebox.showerror(
                "Error", f"Failed to read {os.path.basename(archive)}:\n{e}")
            return []

    def _add_tree_files(self, filepaths, select_last=False):
        """Add tree paths to the file list, skipping duplicates."""
        for fp in filepaths:
            if fp not in self.tree_files:
                self.tree_files.append(fp)
                self.tree_listbox.insert(tk.END, os.path.basename(fp))

        if select_last and self.tree_files:
            self.tree_listbox.selection_clear(0, tk.END)
            idx = len(self.tree_files) - 1
            self.tree_listbox.selection_set(idx)
            self._load_tree(self.tree_files[idx])

    def _clear_trees(self):
        self.tree_files.clear()
        self.tree_listbox.delete(0, tk.END)
        self.current_tree = None
        self.batch_results.clear()
        self.tree_canvas.delete("all")
        self._clear_analysis_state()
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
            # Results carry node references belonging to the tree object that
            # was parsed when they were computed. Selecting the file again
            # reparses it, so recompute against the tree now on screen -
            # cheap, and it keeps the colouring and the tables in step.
            if os.path.basename(filepath) in self.batch_results and self.db.loaded:
                checker = MonophylyChecker(self.db)
                analysis = checker.analyze(self.current_tree)
                checker.label_tree(self.current_tree, analysis["node_results"])
                self._apply_analysis(analysis)
                self._show_results()
            else:
                self._clear_analysis_state()
                self._clear_results()

            self._draw_tree()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse tree:\n{e}")

    def _clear_analysis_state(self):
        self.results = None
        self.taxon_results = []
        self.groups = []
        self.conflicting_tips = []
        self.tip_matches = {}
        self.unresolved_tips = []
        self.highlight_tips = set()
        self.highlight_intruders = set()

    # ---- Tree Drawing ----

    def _layout_tree(self, zoom):
        """Compute canvas coordinates for every node.

        Two layouts, chosen by the "Branch lengths" checkbox:

        - Cladogram (default): x comes from depth, one fixed step per level,
          and tips are aligned in a column so their names read as a list.
          The old code divided a fixed width by the tree depth, which
          squeezed a 30-level tree into 5 px per level.
        - Phylogram: x comes from cumulative branch length, so the drawing
          carries the evolutionary distances the file actually contains.

        Returns (positions, width, height) with positions keyed by id(node).
        """
        tips = self.current_tree.get_tip_nodes()
        tip_spacing = max(4.0, TIP_SPACING * zoom)
        use_lengths = self.show_branch_lengths_var.get()

        y_of = {id(tip): MARGIN_TOP + i * tip_spacing
                for i, tip in enumerate(tips)}

        # Depth and cumulative branch length, top down.
        depth_of = {}
        dist_of = {}
        for node in self.current_tree.pre_order():
            if node.parent is None:
                depth_of[id(node)] = 0
                dist_of[id(node)] = 0.0
            else:
                depth_of[id(node)] = depth_of[id(node.parent)] + 1
                try:
                    length = float(node.branch_length or 0.0)
                except (TypeError, ValueError):
                    length = 0.0
                dist_of[id(node)] = dist_of[id(node.parent)] + max(0.0, length)

        max_depth = max(depth_of.values()) or 1
        max_dist = max(dist_of.values()) or 1.0

        if use_lengths:
            span = PHYLOGRAM_WIDTH * zoom
            x_of = {key: MARGIN_LEFT + (d / max_dist) * span
                    for key, d in dist_of.items()}
        else:
            step = max(LEVEL_STEP_MIN, LEVEL_STEP * zoom)
            x_of = {key: MARGIN_LEFT + d * step for key, d in depth_of.items()}
            # Align tips in a column at the right edge of the cladogram.
            tip_x = MARGIN_LEFT + max_depth * step
            for tip in tips:
                x_of[id(tip)] = tip_x

        # Internal nodes sit at the midpoint of their children, bottom up.
        for node in self.current_tree.post_order():
            if node.children:
                ys = [y_of[id(c)] for c in node.children]
                y_of[id(node)] = (min(ys) + max(ys)) / 2.0

        positions = {key: (x_of[key], y_of[key]) for key in x_of}
        width = max(x for x, _ in positions.values()) + LABEL_GUTTER * zoom
        height = MARGIN_TOP + max(1, len(tips)) * tip_spacing + MARGIN_TOP
        return positions, width, height

    def _branch_color(self, node):
        """Colour a branch by the monophyly status of the clade below it."""
        if node.is_tip or not node.monophyly_status:
            return COLORS["branch"]
        if node.monophyly_status == "monophyletic":
            return COLORS["mono"]
        if node.monophyly_status == "not_monophyletic":
            return COLORS["non_mono"]
        return COLORS["branch"]

    def _draw_tree(self):
        """Draw the current tree on the canvas."""
        canvas = self.tree_canvas
        canvas.delete("all")
        if not self.current_tree:
            return

        zoom = self.zoom_var.get()
        positions, width, height = self._layout_tree(zoom)
        tip_spacing = max(4.0, TIP_SPACING * zoom)
        show_labels = self.show_labels_var.get()
        font_size = max(6, int(9 * zoom))
        # Below this the text is unreadable and just paints the canvas solid.
        room_for_tip_labels = show_labels and tip_spacing >= 7

        highlight = self.highlight_tips
        intruders = self.highlight_intruders

        # Band behind the selected clade, drawn first so branches sit on top.
        if highlight:
            ys = [positions[id(n)][1] for n in self.current_tree.get_tip_nodes()
                  if n.label in highlight]
            if ys:
                canvas.create_rectangle(
                    0, min(ys) - tip_spacing / 2, width,
                    max(ys) + tip_spacing / 2,
                    fill=COLORS["clade_band"], outline="")

        # Branches.
        for node in self.current_tree.pre_order():
            if not node.children:
                continue
            x, y = positions[id(node)]
            for child in node.children:
                cx, cy = positions[id(child)]
                color = self._branch_color(child)
                emphasis = highlight and child.label in highlight
                line_width = (2.6 if emphasis else 1.5) * max(0.6, zoom)
                canvas.create_line(x, y, x, cy, fill=color, width=line_width)
                canvas.create_line(x, cy, cx, cy, fill=color, width=line_width)

        # Tip labels.
        if room_for_tip_labels:
            for tip in self.current_tree.get_tip_nodes():
                x, y = positions[id(tip)]
                match = self.tip_matches.get(tip.label)
                if tip.label in intruders:
                    color, weight = COLORS["non_mono"], "bold"
                elif match is None:
                    color, weight = COLORS["unresolved"], "normal"
                elif not match.is_exact:
                    color, weight = COLORS["inexact"], "normal"
                else:
                    color, weight = COLORS["tip"], "normal"
                canvas.create_text(
                    x + 6, y, text=tip.label, anchor="w", fill=color,
                    font=("Helvetica", font_size, weight))

        if not show_labels:
            canvas.configure(scrollregion=(0, 0, width, height))
            return

        # Support values, small and grey, above and left of their node.
        if self.show_support_var.get() and tip_spacing >= 9:
            for node in self.current_tree.pre_order():
                if node.is_tip or not node.support:
                    continue
                x, y = positions[id(node)]
                canvas.create_text(
                    x - 3, y - 6 * max(0.7, zoom), text=str(node.support),
                    anchor="e", fill=COLORS["support"],
                    font=("Helvetica", max(6, font_size - 2)))

        # Clade names, drawn to the RIGHT of the node on an opaque chip.
        # Drawing them to the left clipped every shallow clade off the canvas.
        for node in self.current_tree.pre_order():
            if node.is_tip or not node.label:
                continue
            x, y = positions[id(node)]
            color = self._branch_color(node)
            if color == COLORS["branch"]:
                color = COLORS["unresolved"]
            text_id = canvas.create_text(
                x + 5, y, text=node.label, anchor="w", fill=color,
                font=("Helvetica", max(6, font_size - 1), "bold"))
            x0, y0, x1, y1 = canvas.bbox(text_id)
            chip = canvas.create_rectangle(
                x0 - 2, y0 - 1, x1 + 2, y1 + 1,
                fill=COLORS["bg"], outline=COLORS["chip_edge"])
            canvas.tag_lower(chip, text_id)

        canvas.configure(scrollregion=(0, 0, width, height))

    def _redraw_tree(self):
        """Coalesce rapid redraws (the zoom slider fires per pixel)."""
        if self._redraw_job is not None:
            self.root.after_cancel(self._redraw_job)
        self._redraw_job = self.root.after(60, self._do_redraw)

    def _do_redraw(self):
        self._redraw_job = None
        self._draw_tree()

    def _scroll_canvas_to(self, y_center, height):
        """Scroll the tree canvas so a given y is roughly centred."""
        if height <= 0:
            return
        view = self.tree_canvas.winfo_height() or 1
        target = (y_center - view / 2) / height
        self.tree_canvas.yview_moveto(max(0.0, min(1.0, target)))

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
            analysis = checker.analyze(
                self.current_tree, progress_callback=self._update_status)
            checker.label_tree(self.current_tree, analysis["node_results"])
            self._apply_analysis(analysis)

            # Store for batch results
            sel = self.tree_listbox.curselection()
            if sel:
                basename = os.path.basename(self.tree_files[sel[0]])
                self.batch_results[basename] = analysis

            self.root.after(0, self._show_results)
            self.root.after(0, self._draw_tree)

        self._run_threaded(analyze, "Running monophyly analysis...")

    def _apply_analysis(self, analysis):
        """Adopt an analysis result as the current view state."""
        self.results = analysis["node_results"]
        self.taxon_results = analysis["taxon_results"]
        self.groups = analysis["groups"]
        self.conflicting_tips = analysis["conflicting_tips"]
        self.tip_matches = analysis["matches"]
        self.unresolved_tips = analysis["unresolved"]
        self.highlight_tips = set()
        self.highlight_intruders = set()

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
                    analysis = checker.analyze(tree)
                    checker.label_tree(tree, analysis["node_results"])
                    self.batch_results[basename] = analysis
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
        self.summary_text.insert(tk.END, checker.get_summary(
            self.results, self.taxon_results, self.tip_matches))

        self._populate_group_tables()
        self._populate_tip_matching()

        if any(not g.is_monophyletic for g in self.groups):
            self.results_notebook.select(self.nonmono_frame)

    def _remember_row(self, table, row, payload):
        """Associate a table row with what it stands for."""
        self._row_results[(str(table), row)] = payload

    def _matches_rank_filter(self, group):
        """True if any taxon merged into this row has the wanted rank."""
        wanted = self.rank_filter_var.get()
        if wanted == "all ranks":
            return True
        return any((t.rank or "no rank") == wanted for t in group.taxa)

    def _populate_group_tables(self):
        """Fill the group tables, honouring the rank filter.

        Rows are ordered by how few tips would have to move, so the groups a
        biologist can actually act on come first. Sorting by intruder count
        instead put a group whose two tips merely are not sister at the top
        purely because its enclosing clade happened to be large.
        """
        for table in (self.nonmono_tree, self.mono_tree, self.conflict_tree):
            for item in table.get_children():
                table.delete(item)
        self._row_results.clear()
        if not self.taxon_results:
            return

        rows = self.groups if self.collapse_var.get() else [
            CollapsedGroup([t]) for t in self.taxon_results]
        wanted_conflict = self.conflict_filter_var.get()

        broken = [g for g in rows if not g.is_monophyletic]
        intact = [g for g in rows if g.is_monophyletic]
        shown_broken = shown_intact = 0

        for group in sorted(broken, key=lambda g: g.severity_sort_key()):
            if not self._matches_rank_filter(group):
                continue
            if (wanted_conflict != "all conflicts"
                    and group.conflict_type != wanted_conflict):
                continue
            listed = ", ".join(group.intruders[:12])
            if len(group.intruders) > 12:
                listed += f"  ... (+{len(group.intruders) - 12} more)"
            row = self.nonmono_tree.insert("", tk.END, values=(
                group.name, group.rank, group.member_count,
                group.conflict_cost, group.conflict_type,
                group.culprit or "", len(group.intruders),
                len(group.outliers), listed))
            self._remember_row(self.nonmono_tree, row, group)
            self._add_constituent_rows(self.nonmono_tree, row, group, 9)
            shown_broken += 1

        for group in sorted(intact, key=lambda g: g.rank_sort_key()):
            if not self._matches_rank_filter(group):
                continue
            notes = []
            if group.is_collapsed:
                notes.append(f"{len(group.taxa)} nested names")
            if group.unknown_inside:
                notes.append(f"{len(group.unknown_inside)} unplaced tips inside")
            row = self.mono_tree.insert("", tk.END, values=(
                group.name, group.rank, group.member_count, "; ".join(notes)))
            self._remember_row(self.mono_tree, row, group)
            self._add_constituent_rows(self.mono_tree, row, group, 4)
            shown_intact += 1

        self._populate_conflicting_tips()

        self.results_notebook.tab(
            self.nonmono_frame, text=f"Broken Groups ({shown_broken})")
        self.results_notebook.tab(
            self.mono_frame, text=f"Monophyletic ({shown_intact})")

    def _add_constituent_rows(self, table, parent_row, group, column_count):
        """Expand a merged row to show the taxa it stands for."""
        if not group.is_collapsed:
            return
        for taxon in group.taxa:
            values = [taxon.name, taxon.rank or "no rank"] + [""] * (column_count - 2)
            child = table.insert(parent_row, tk.END, values=values)
            self._remember_row(table, child, group)

    def _populate_conflicting_tips(self):
        """List the tips whose removal would repair the most groups."""
        repairing = [t for t in self.conflicting_tips if t.repairs]
        for tip in repairing:
            names = ", ".join(g.name for g in tip.repairs[:8])
            if len(tip.repairs) > 8:
                names += f"  ... (+{len(tip.repairs) - 8} more)"
            row = self.conflict_tree.insert("", tk.END, values=(
                tip.label, len(tip.repairs), tip.role,
                len(tip.appears_in), names))
            self._remember_row(self.conflict_tree, row, tip)

        broken = sum(1 for g in self.groups if not g.is_monophyletic)
        fixable = sum(len(t.repairs) for t in repairing)
        if repairing:
            self.conflict_banner.configure(
                text=f"Rechecking {_plural(len(repairing), 'sequence')} "
                     f"would repair {fixable} of {broken} conflicting "
                     f"groups.")
        elif broken:
            self.conflict_banner.configure(
                text="No single tip explains any conflict here - every one "
                     "needs 2 or more tips to move.")
        else:
            self.conflict_banner.configure(text="")

        self.results_notebook.tab(
            self.tips_frame, text=f"Conflicting Tips ({len(repairing)})")

    def _populate_tip_matching(self):
        """Report how each tip was matched, not just which ones failed."""
        text = self.unresolved_text
        inexact = sorted((m for m in self.tip_matches.values()
                          if not m.is_exact), key=lambda m: m.label)
        ambiguous = sorted((m for m in self.tip_matches.values()
                            if m.ambiguous), key=lambda m: m.label)

        total = len(self.tip_matches) + len(self.unresolved_tips)
        text.insert(tk.END, f"{len(self.tip_matches)} of {total} tips matched "
                            f"NCBI taxonomy.\n\n")

        if self.unresolved_tips:
            text.insert(tk.END, f"--- No match ({len(self.unresolved_tips)}) ---\n")
            text.insert(tk.END, "These tips are excluded from the analysis.\n\n")
            for tip in sorted(self.unresolved_tips):
                text.insert(tk.END, f"  {tip}\n")
            text.insert(tk.END, "\n")

        if inexact:
            text.insert(tk.END, f"--- Matched loosely ({len(inexact)}) ---\n")
            text.insert(tk.END, "Name tokens were dropped to find a match, so "
                                "these sit at a coarser rank than the label "
                                "suggests.\n\n")
            for m in inexact:
                text.insert(tk.END, f"  {m.label}\n      -> {m.matched_name} "
                                    f"({m.rank or 'no rank'}, via {m.method})\n")
            text.insert(tk.END, "\n")

        if ambiguous:
            text.insert(tk.END, f"--- Ambiguous names ({len(ambiguous)}) ---\n")
            text.insert(tk.END, "Several unrelated taxa share these names; the "
                                "first NCBI entry was used.\n\n")
            for m in ambiguous:
                others = self.db.get_homonyms(m.matched_name)
                names = ", ".join(
                    f"{self.db.get_scientific_name(c)} "
                    f"({self.db.code_to_rank.get(c, '?')})" for c in others[:4])
                text.insert(tk.END, f"  {m.label}: {names}\n")
            text.insert(tk.END, "\n")

        if not self.unresolved_tips and not inexact and not ambiguous:
            text.insert(tk.END, "Every tip matched exactly. Nothing to review.\n")

    def _on_result_select(self, event):
        """Highlight what the selected row refers to and scroll the tree."""
        selection = event.widget.selection()
        if not selection or not self.current_tree:
            return
        entry = self._row_results.get((str(event.widget), selection[0]))
        if entry is None:
            return

        if isinstance(entry, ConflictingTip):
            self.highlight_tips = {entry.label}
            self.highlight_intruders = {entry.label}
            status = (f"{entry.label}: removing it would repair "
                      f"{_plural(len(entry.repairs), 'group')}; it intrudes "
                      f"into {len(entry.appears_in)}")
        else:
            self.highlight_tips = set(entry.clade_tips)
            self.highlight_intruders = set(entry.intruders)
            status = (f"{entry.name} ({entry.rank}): "
                      f"{_plural(entry.member_count, 'member')}")
            if entry.is_monophyletic:
                status += " - monophyletic"
            else:
                status += (f", {_plural(len(entry.intruders), 'intruder')}, "
                           f"{_plural(len(entry.outliers), 'outlier')}"
                           f" - move {_plural(entry.conflict_cost, 'tip')} "
                           f"({entry.conflict_type})")

        self._draw_tree()
        self._scroll_to_highlight()
        self._update_status(status)

    def _scroll_to_highlight(self):
        """Centre the canvas on whatever is currently highlighted."""
        positions, _width, height = self._layout_tree(self.zoom_var.get())
        ys = [positions[id(n)][1] for n in self.current_tree.get_tip_nodes()
              if n.label in self.highlight_tips]
        if ys:
            self._scroll_canvas_to((min(ys) + max(ys)) / 2, height)

    def _clear_results(self):
        self.summary_text.delete("1.0", tk.END)
        for table in (self.nonmono_tree, self.mono_tree, self.conflict_tree):
            for item in table.get_children():
                table.delete(item)
        self._row_results.clear()
        self.conflict_banner.configure(text="")
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
            # Always write branch lengths and support: the drawing checkboxes
            # control the view, not what lands in the exported file.
            newick = self.current_tree.to_newick(
                branch_lengths=True, label_format="taxon_support") + ";"
            with open(filepath, "w", encoding="utf-8") as f:
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
            report = checker.get_summary(
                self.results, self.taxon_results, self.tip_matches)

            if self.unresolved_tips:
                report += "\n\n--- Unresolved Tips ---\n"
                for tip in sorted(self.unresolved_tips):
                    report += f"  {tip}\n"

            with open(filepath, "w", encoding="utf-8") as f:
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

    def _load_project_guide(self):
        """Load the onboarding guide from disk, with a safe fallback."""
        guide_path = Path(__file__).resolve().parent / "docs" / "ONBOARDING.md"
        try:
            return guide_path.read_text(encoding="utf-8")
        except OSError as exc:
            return (
                "Project guide could not be loaded.\n\n"
                f"Expected path: {guide_path}\n"
                f"Error: {exc}\n"
            )

    def _style_text_widget(self, widget):
        """Apply consistent colors to tk text widgets."""
        widget.configure(
            bg=COLORS["field_bg"],
            fg=COLORS["field_fg"],
            insertbackground=COLORS["field_fg"],
            selectbackground=COLORS["highlight"],
            selectforeground=COLORS["button_fg"],
            highlightbackground=COLORS["panel_edge"],
            highlightcolor=COLORS["panel_edge"],
            relief=tk.FLAT,
        )

    # ---- Helpers ----


    def _bind_mousewheel(self, widget, yview=True, xview=False):
        """Enable mouse wheel scrolling for a widget."""
        def on_scroll(event):
            step = mousewheel_units(event)
            if step == 0:
                return
            use_x = xview and (event.state & 0x0001)
            if use_x:
                widget.xview_scroll(step, "units")
            elif yview:
                widget.yview_scroll(step, "units")
            return "break"

        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>",
                    "<Shift-MouseWheel>", "<Shift-Button-4>",
                    "<Shift-Button-5>"):
            widget.bind(seq, on_scroll, add="+")

    def _bind_scroll_wheels(self):
        """Attach mouse wheel scrolling to all scrollable panels."""
        self._bind_mousewheel(self.tree_listbox, yview=True)
        self._bind_mousewheel(self.tree_canvas, yview=True, xview=True)
        self._bind_mousewheel(self.summary_text, yview=True)
        self._bind_mousewheel(self.nonmono_tree, yview=True)
        self._bind_mousewheel(self.mono_tree, yview=True)
        self._bind_mousewheel(self.conflict_tree, yview=True)
        self._bind_mousewheel(self.unresolved_text, yview=True)

    def _search_tips(self):
        """Highlight tips matching the search box and scroll to the first."""
        query = self.search_var.get().strip().lower()
        if not self.current_tree:
            return
        if not query:
            self.highlight_tips = set()
            self._redraw_tree()
            return

        matches = [n for n in self.current_tree.get_tip_nodes()
                   if query in n.label.lower()]
        self.highlight_tips = {n.label for n in matches}
        self._draw_tree()

        if matches:
            positions, _width, height = self._layout_tree(self.zoom_var.get())
            self._scroll_canvas_to(positions[id(matches[0])][1], height)
            self._update_status(f"{len(matches)} tip(s) match '{query}'")
        else:
            self._update_status(f"No tips match '{query}'")

    def _set_busy(self, busy):
        """Disable the analysis buttons while work is in flight.

        Without this a second run can start mid-analysis and the two threads
        overwrite each other's results.
        """
        self._busy = busy
        state = tk.DISABLED if busy else tk.NORMAL
        for widget in (self.run_button, self.batch_button):
            widget.configure(state=state)

    def _run_threaded(self, target, status_msg="Working..."):
        """Run a function in a background thread."""
        self._update_status(status_msg)
        self.root.after(0, lambda: self._set_busy(True))

        def wrapper():
            try:
                target()
            except Exception as e:
                message = str(e)
                self.root.after(0, lambda: messagebox.showerror(
                    "Error", message))
            finally:
                self.root.after(0, lambda: self._set_busy(False))
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
