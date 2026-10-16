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

import json
import math
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
from monophyly import (DEFAULT_MIN_BRANCH_LENGTH, MAIN_RANKS, CollapsedGroup,
                       ConflictingTip, MonophylyChecker, plural as _plural)
from help_content import HELP_SECTIONS, PANEL_HELP

# File-dialog filters, kept next to each other so CLI and GUI stay in sync.
# macOS Tk only compares a file's last extension against these patterns, so
# "*.tre.gz" alone never matches tree.tre.gz there. The bare "*.gz *.bz2 *.xz"
# entries make compressed files selectable everywhere; the content is checked
# by magic bytes after picking anyway.
COMPRESSED_PATTERNS = "*.gz *.bz2 *.xz"
ARCHIVE_PATTERNS = "*.tar.gz *.tgz *.tar *.zip *.gz"

TREE_FILETYPES = [
    ("Tree files and archives",
     "*.nwk *.tre *.tree *.newick *.tgz *.tar *.zip " + COMPRESSED_PATTERNS),
    ("Newick files", "*.nwk *.tre *.tree *.newick"),
    ("Compressed Newick",
     "*.nwk.gz *.tre.gz *.tree.gz *.newick.gz " + COMPRESSED_PATTERNS),
    ("Tree archives", ARCHIVE_PATTERNS),
    ("All files", "*.*"),
]

TAXONOMY_FILETYPES = [
    ("Taxonomy dump or archive",
     "*.dmp *.tgz *.tar *.zip " + COMPRESSED_PATTERNS),
    ("DMP files", "*.dmp *.dmp.gz *.gz"),
    ("NCBI archives", ARCHIVE_PATTERNS),
    ("All files", "*.*"),
]


# ---- Tree layout ----
TIP_SPACING = 16          # px between adjacent tips at zoom 1.0
LEVEL_STEP = 26           # px per depth level in cladogram mode
LEVEL_STEP_MIN = 10       # never squeeze levels below this
PHYLOGRAM_WIDTH = 620     # px the longest root-to-tip path spans
MARGIN_LEFT = 46          # room for support values left of the root
MARGIN_TOP = 20
LABEL_GUTTER_PAD = 24     # breathing room right of the longest tip name

# The file list grows with its contents between these bounds. Below the max
# the panel is only as deep as it needs to be, so a session with one tree
# stops reserving a third of the left column for nothing.
TREE_LIST_MIN_ROWS = 3
TREE_LIST_MAX_ROWS = 12

TICK_ON = "☑"        # ballot box with check
TICK_OFF = "☐"       # empty ballot box
TICK_HIT_WIDTH = 26       # px from a row's left edge that toggles the tick

SASH_THICKNESS = 9        # thick enough to aim at, thin enough to ignore
SASH_HIT_SLOP = 14        # px from the sash that still counts as hitting it
SASH_EDGE_ZONE = 26       # within this of an end, a pane counts as collapsed

# Circular layout. The sweep stops short of a full turn so the root stays
# visible and the reading order is obvious.
CIRCLE_RADIUS = 340       # outer radius at zoom 1.0
CIRCLE_INNER = 18         # the root is a disc, not a point
CIRCLE_START = 95         # degrees; just past 12 o'clock
CIRCLE_SWEEP = 350        # degrees of arc the tips occupy
CIRCLE_LABEL_OFFSET = 14  # push clade names out along their own spoke

# Discrete zoom stops, so each click is a visible change and the ends are
# hard limits.
ZOOM_LEVELS = (0.1, 0.15, 0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 2.8, 4.0)
ZOOM_DEFAULT = 1.0


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


# ---- Colour themes ----
# COLORS is mutated in place by apply_theme(), so every widget and drawing
# routine keeps reading the same dict and picks up the change without being
# handed a new one. Verified: nothing reads COLORS at import time.
#
# Every theme must define the same keys. The semantic ones - mono, non_mono,
# unresolved, inexact - carry meaning, not decoration, so a theme may restyle
# them but must keep them mutually distinguishable.

THEMES = {}

THEMES["Dark (default)"] = {
    "conflict_dash": None,
    "mono": "#2d8a4e", "non_mono": "#d63031", "unresolved": "#6c757d",
    "inexact": "#b8860b", "support": "#8a8a8a", "clade_band": "#fff3bf",
    "chip_edge": "#dddddd", "tip": "#2c3e50", "branch": "#555555",
    "bg": "#ffffff", "highlight": "#ffeaa7",
    "header_bg": "#34495e", "header_fg": "#ecf0f1",
    "app_bg": "#2f2f2f", "panel_bg": "#3a3a3a", "panel_edge": "#565656",
    "text": "#f1f1f1", "muted_text": "#d0d0d0",
    "button_bg": "#f2f2f2", "button_fg": "#1f1f1f",
    "tab_bg": "#4a4a4a", "tab_active": "#f2f2f2", "tab_active_fg": "#202020",
    "field_bg": "#1f1f1f", "field_fg": "#f5f5f5",
    "info_icon": "#7fb3ff", "tooltip_bg": "#f7f7f2", "tooltip_fg": "#1f1f1f",
}

THEMES["Light (cream)"] = {
    "conflict_dash": None,
    "mono": "#2f7d4f", "non_mono": "#c0392b", "unresolved": "#8d8577",
    "inexact": "#a8761b", "support": "#9a9284", "clade_band": "#f6e7bd",
    "chip_edge": "#d8cdb4", "tip": "#3b342a", "branch": "#7a7266",
    "bg": "#fffdf6", "highlight": "#f2dfa4",
    "header_bg": "#e8ddc4", "header_fg": "#3b342a",
    "app_bg": "#f4ecd8", "panel_bg": "#efe5cd", "panel_edge": "#d6c9ab",
    "text": "#3b342a", "muted_text": "#6f6657",
    "button_bg": "#fffdf6", "button_fg": "#3b342a",
    "tab_bg": "#e4d8ba", "tab_active": "#fffdf6", "tab_active_fg": "#3b342a",
    "field_bg": "#fffdf6", "field_fg": "#3b342a",
    "info_icon": "#1f6f8b", "tooltip_bg": "#fffdf6", "tooltip_fg": "#3b342a",
}

THEMES["Coffee (warm brown)"] = {
    "conflict_dash": None,
    "mono": "#6a9a4a", "non_mono": "#d1603d", "unresolved": "#8a7c6d",
    "inexact": "#c99a3c", "support": "#9c8c7a", "clade_band": "#5a4632",
    "chip_edge": "#6b563f", "tip": "#f0e4d4", "branch": "#8a7562",
    "bg": "#2b211a", "highlight": "#7a5c34",
    "header_bg": "#3d2f24", "header_fg": "#f0e4d4",
    "app_bg": "#2a211b", "panel_bg": "#352a22", "panel_edge": "#4d3d31",
    "text": "#f0e4d4", "muted_text": "#c3b2a0",
    "button_bg": "#e8d9c5", "button_fg": "#2a211b",
    "tab_bg": "#41332a", "tab_active": "#e8d9c5", "tab_active_fg": "#2a211b",
    "field_bg": "#1f1813", "field_fg": "#f0e4d4",
    "info_icon": "#e0a860", "tooltip_bg": "#f4e9d8", "tooltip_fg": "#2a211b",
}

THEMES["Black and white"] = {
    # No hue at all, so conflict cannot ride on red-vs-green. It rides on a
    # dashed line plus weight instead: agreeing branches are mid-grey and
    # solid, conflicting ones black and dashed. Both survive greyscale
    # printing and photocopying, which colour never does.
    "conflict_dash": (5, 3),
    # The four greys are spread on lightness alone, since that is the only
    # channel this theme has. #606060 for inexact sat 5 dE from the #6e6e6e
    # of mono - indistinguishable to anyone, colour vision or not - which the
    # legend puts side by side and makes obvious.
    "mono": "#6e6e6e", "non_mono": "#000000", "unresolved": "#909090",
    "inexact": "#3a3a3a", "support": "#707070", "clade_band": "#e6e6e6",
    "chip_edge": "#000000", "tip": "#000000", "branch": "#000000",
    "bg": "#ffffff", "highlight": "#d0d0d0",
    "header_bg": "#000000", "header_fg": "#ffffff",
    "app_bg": "#ffffff", "panel_bg": "#f2f2f2", "panel_edge": "#000000",
    "text": "#000000", "muted_text": "#4a4a4a",
    "button_bg": "#ffffff", "button_fg": "#000000",
    "tab_bg": "#e0e0e0", "tab_active": "#000000", "tab_active_fg": "#ffffff",
    "field_bg": "#ffffff", "field_fg": "#000000",
    "info_icon": "#000000", "tooltip_bg": "#ffffff", "tooltip_fg": "#000000",
}

THEMES["Terminal green"] = {
    "conflict_dash": None,
    "mono": "#39ff14", "non_mono": "#ff4d4d", "unresolved": "#4f6f4f",
    "inexact": "#d7ff5c", "support": "#4e8a4e", "clade_band": "#0d3f16",
    "chip_edge": "#1f7a2e", "tip": "#9dff9d", "branch": "#2f8f3f",
    "bg": "#04120a", "highlight": "#14572a",
    "header_bg": "#0a2413", "header_fg": "#9dff9d",
    "app_bg": "#05170c", "panel_bg": "#0a2413", "panel_edge": "#1f7a2e",
    "text": "#9dff9d", "muted_text": "#5fbf6f",
    "button_bg": "#0f3d1e", "button_fg": "#c6ffc6",
    "tab_bg": "#0a2413", "tab_active": "#1f7a2e", "tab_active_fg": "#c6ffc6",
    "field_bg": "#04120a", "field_fg": "#9dff9d",
    "info_icon": "#39ff14", "tooltip_bg": "#0a2413", "tooltip_fg": "#9dff9d",
}

THEMES["Colour-blind safe"] = {
    "conflict_dash": None,
    # Conflict status normally rides on red vs green, which is exactly the
    # pair deuteranopia and protanopia collapse - roughly 8% of men. Orange
    # and blue stay distinct under every common form of colour blindness.
    "mono": "#0072b2", "non_mono": "#d55e00", "unresolved": "#7a7a7a",
    "inexact": "#cc79a7", "support": "#8a8a8a", "clade_band": "#ffe9c9",
    "chip_edge": "#c8c8c8", "tip": "#1b2a38", "branch": "#5a5a5a",
    "bg": "#ffffff", "highlight": "#f0e442",
    "header_bg": "#1b2a38", "header_fg": "#ffffff",
    "app_bg": "#eef1f4", "panel_bg": "#e2e7ec", "panel_edge": "#b4bec8",
    "text": "#1b2a38", "muted_text": "#4f5b66",
    "button_bg": "#ffffff", "button_fg": "#1b2a38",
    "tab_bg": "#d5dce3", "tab_active": "#ffffff", "tab_active_fg": "#1b2a38",
    "field_bg": "#ffffff", "field_fg": "#1b2a38",
    "info_icon": "#0072b2", "tooltip_bg": "#ffffff", "tooltip_fg": "#1b2a38",
}

THEMES["Colour-blind safe (dark)"] = {
    "conflict_dash": None,
    # The same Okabe-Ito reasoning as the light version, lifted onto a dark
    # canvas. #0072b2 reads as near-black against a dark background, so
    # agreement moves up to sky blue and conflict to orange; that pair stays
    # separable under deuteranopia, protanopia and tritanopia alike. Pink for
    # inexact and neutral grey for unresolved keep all four statuses apart on
    # lightness as well as hue.
    "mono": "#56b4e9", "non_mono": "#e69f00", "unresolved": "#707c87",
    "inexact": "#d4a5c9", "support": "#7f8b98", "clade_band": "#24374a",
    "chip_edge": "#4a5a6a", "tip": "#e9eef3", "branch": "#8592a0",
    "bg": "#10161c", "highlight": "#f0e442",
    "header_bg": "#1b2733", "header_fg": "#eef3f7",
    "app_bg": "#141b22", "panel_bg": "#1d252e", "panel_edge": "#37424e",
    "text": "#e9eef3", "muted_text": "#b0bcc7",
    "button_bg": "#dfe6ec", "button_fg": "#141b22",
    "tab_bg": "#222c37", "tab_active": "#dfe6ec", "tab_active_fg": "#141b22",
    "field_bg": "#0d1319", "field_fg": "#eef3f7",
    "info_icon": "#56b4e9", "tooltip_bg": "#eef3f7", "tooltip_fg": "#141b22",
}

DEFAULT_THEME = "Dark (default)"

# Themes that have been renamed. Without this a settings file naming the old
# theme would fall through to the default and quietly discard a choice the
# user had already made.
THEME_ALIASES = {"Dark (warm brown)": "Coffee (warm brown)"}

# The live palette. Mutated, never rebound.
COLORS = dict(THEMES[DEFAULT_THEME])


def apply_theme(name):
    """Swap the live palette in place. Returns the name actually applied."""
    name = THEME_ALIASES.get(name, name)
    palette = THEMES.get(name)
    if palette is None:
        name, palette = DEFAULT_THEME, THEMES[DEFAULT_THEME]
    COLORS.clear()
    COLORS.update(palette)
    return name



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
                                   style="Stepper.TButton",
                                   command=lambda: self.zoom(-1))
        self._smaller.pack(side=tk.LEFT, padx=(6, 0))
        Tooltip(self._smaller, "Smaller text  (Ctrl/Cmd −)")
        self._bigger = ttk.Button(bar, text="A+", width=3,
                                  style="Stepper.TButton",
                                  command=lambda: self.zoom(1))
        self._bigger.pack(side=tk.LEFT, padx=(1, 0))
        Tooltip(self._bigger, "Larger text  (Ctrl/Cmd +)")
        # Readout first, then Reset: the number belongs to the stepper that
        # changes it, not stranded on the far side of another button.
        self._label = ttk.Label(bar, text="", width=6, anchor="center",
                                style="Readout.TLabel")
        self._label.pack(side=tk.LEFT, padx=(8, 6))
        ttk.Button(bar, text="Reset", style="Quiet.TButton",
                   command=self.reset).pack(side=tk.LEFT)
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


class LabelSpace:
    """Tracks the rectangles labels already occupy, so they can avoid each other.

    Bucketed into horizontal bands: a plain list means every new label scans
    every placed one, which gets slow on a tree with hundreds of clades.
    """

    BAND = 24

    def __init__(self):
        self._bands = {}

    @staticmethod
    def pad(bbox, x=2, y=1):
        x0, y0, x1, y1 = bbox
        return (x0 - x, y0 - y, x1 + x, y1 + y)

    def _bands_for(self, y0, y1):
        return range(int(y0 // self.BAND), int(y1 // self.BAND) + 1)

    def collides(self, rect):
        x0, y0, x1, y1 = rect
        for band in self._bands_for(y0, y1):
            for a0, b0, a1, b1 in self._bands.get(band, ()):
                if x0 < a1 and a0 < x1 and y0 < b1 and b0 < y1:
                    return True
        return False

    def add(self, rect):
        for band in self._bands_for(rect[1], rect[3]):
            self._bands.setdefault(band, []).append(rect)


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


class LegendBar(tk.Canvas):
    """Key to the marks the tree draws, sitting directly under the tree.

    A canvas rather than a row of ttk.Labels because not every distinction
    the tree makes is a colour. The black-and-white theme carries conflict on
    a dash pattern, and an intruding tip is called out by weight as well as
    hue - a coloured square would quietly misreport both. Each swatch is
    drawn the same way the tree draws the thing it stands for.

    Colours are read from COLORS at draw time and never cached, so a theme
    switch is one redraw away and the colour-blind palettes are described in
    their own terms rather than in the default theme's.
    """

    # (swatch kind, palette key, caption). Order runs branches, then tips,
    # then selection - the same order the eye meets them on the canvas.
    ENTRIES = (
        ("line", "mono", "Monophyletic"),
        ("conflict", "non_mono", "Not monophyletic"),
        ("line", "branch", "Not assessed"),
        ("bold", "non_mono", "Intruding tip"),
        ("text", "inexact", "Inexact match"),
        ("text", "unresolved", "Unmatched tip"),
        ("band", "clade_band", "Selected clade"),
    )

    SWATCH_W = 26       # long enough to show two dashes of a dash pattern
    GAP = 5             # swatch to its caption
    SPACING = 15        # between entries
    PAD = 5
    ROW_H = 18

    def __init__(self, parent, **pack_kwargs):
        super().__init__(parent, highlightthickness=0,
                         height=self.ROW_H + 2 * self.PAD)
        self._pack_kwargs = pack_kwargs
        self._font = tkfont.Font(family="Helvetica", size=9)
        self._heading_font = tkfont.Font(family="Helvetica", size=9,
                                         weight="bold")
        self._height = None
        self.bind("<Configure>", lambda _e: self.redraw())

    # ---- visibility ----

    def show(self):
        """Pack back into the slot it was built for, not onto the end."""
        self.pack(**self._pack_kwargs)
        self.redraw()

    def hide(self):
        self.pack_forget()

    def visible(self):
        return bool(self.winfo_manager())

    # ---- drawing ----

    def redraw(self):
        """Lay the entries out across the available width, wrapping as needed."""
        self.delete("all")
        self.configure(bg=COLORS["bg"])
        width = self.winfo_width()
        if width <= 1:
            return      # not laid out yet; the <Configure> binding calls back

        heading = "Key"
        self.create_text(self.PAD, self.PAD + self.ROW_H / 2, text=heading,
                         anchor="w", fill=COLORS["support"],
                         font=self._heading_font)
        left = self.PAD + self._heading_font.measure(heading) + self.SPACING

        x, y, rows = left, self.PAD, 1
        for kind, key, caption in self.ENTRIES:
            entry_w = self.SWATCH_W + self.GAP + self._font.measure(caption)
            if x > left and x + entry_w > width - self.PAD:
                x, y = self.PAD, y + self.ROW_H
                rows += 1
            self._draw_entry(kind, key, caption, x, y)
            x += entry_w + self.SPACING

        # Only resize when the row count actually changed: configure() fires
        # another <Configure>, and reasserting the same height every time
        # would leave the two bouncing off each other forever.
        needed = rows * self.ROW_H + 2 * self.PAD
        if needed != self._height:
            self._height = needed
            self.configure(height=needed)

    def _draw_entry(self, kind, key, caption, x, y):
        colour = COLORS[key]
        mid = y + self.ROW_H / 2
        if kind in ("line", "conflict"):
            dash = COLORS.get("conflict_dash") if kind == "conflict" else None
            self.create_line(x, mid, x + self.SWATCH_W, mid, fill=colour,
                             width=2.5, **({"dash": dash} if dash else {}))
        elif kind == "band":
            self.create_rectangle(x, mid - 5, x + self.SWATCH_W, mid + 5,
                                  fill=colour, outline=COLORS["chip_edge"])
        else:
            self.create_text(
                x + self.SWATCH_W / 2, mid, text="Aa", fill=colour,
                font=("Helvetica", 9, "bold" if kind == "bold" else "normal"))
        # Captions take the canvas text colour, not the panel's: this widget
        # sits on the tree's background, which in several themes is the
        # opposite lightness to the surrounding chrome.
        self.create_text(x + self.SWATCH_W + self.GAP, mid, text=caption,
                         anchor="w", fill=COLORS["tip"], font=self._font)


class PhyLabelerApp:
    """Main application window."""

    def __init__(self, root):
        self.root = root
        self.root.title("PhyLabeler - Phylogenetic Tree Labeler")
        self.root.geometry("1200x800")
        self.root.minsize(900, 600)

        # State
        self.db = TaxonomyDB()
        settings = self._load_settings()
        self.theme_name = apply_theme(settings.get("theme", DEFAULT_THEME))
        self.theme_var = tk.StringVar(value=self.theme_name)
        self.show_legend_var = tk.BooleanVar(
            value=bool(settings.get("legend", True)))
        # What counts as a weak branch; see MonophylyChecker.
        self.min_branch_length = settings.get(
            "min_branch_length", DEFAULT_MIN_BRANCH_LENGTH)
        self.min_support = settings.get("min_support")
        self._configure_styles()
        self.current_tree = None
        self.results = None
        self.unresolved_tips = []
        self.tree_files = []        # For batch mode
        self.file_enabled = {}      # path -> included in batch runs
        self.batch_results = {}     # filename -> analysis dict
        self.taxon_results = []     # One entry per named taxonomic group
        self.groups = []            # Taxa merged when they cover the same tips
        self.conflicting_tips = []  # Tips ranked by groups they would repair
        self.tip_matches = {}       # tip label -> TipMatch
        self.highlight_tips = set()      # Tips of the selected clade
        self.highlight_intruders = set()  # Tips that break the selected group
        self.search_hits = set()         # id(node) for every search match
        self._temp_dirs = []        # Scratch dirs for archive extraction
        self._redraw_job = None     # Pending debounced canvas redraw
        self._hidden_labels = 0     # Clade chips omitted for lack of room
        self.tree_expanded = False  # Tree pane given the whole column
        self._results_sash = None   # Sash position to restore
        self._sash_memory = {}      # Pane key -> width to restore on reopen
        self._fullscreen = False
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
        self.root.after_idle(self._fit_vertical_split)
        self._check_cache()

    def _fit_vertical_split(self):
        """Open the tree/results split so the controls column fits unscrolled.

        The column can scroll, but the primary action lives at its bottom and
        a Run button you have to go looking for is a worse default than a
        slightly shorter results panel. Capped, so a tall controls column
        cannot squeeze the results away entirely.
        """
        try:
            self.root.update_idletasks()
            total = self.right_pane.winfo_height()
            if total <= 1:
                return
            needed = self.controls_frame.winfo_reqheight() + 24
            self.right_pane.sashpos(
                0, max(int(total * 0.5), min(needed, int(total * 0.72))))
        except tk.TclError:
            pass    # window torn down before it was ever mapped

    def _sync_left_scroll(self):
        """Keep the controls column's scrollregion and scrollbar honest."""
        canvas, inner = self.left_canvas, self.controls_frame
        canvas.itemconfigure(self._controls_window,
                             width=max(canvas.winfo_width(), 1))
        canvas.configure(scrollregion=canvas.bbox("all"))
        overflowing = inner.winfo_reqheight() > canvas.winfo_height()
        if overflowing and not self.left_scroll.winfo_manager():
            self.left_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        elif not overflowing and self.left_scroll.winfo_manager():
            self.left_scroll.pack_forget()
            canvas.yview_moveto(0.0)   # nothing above to scroll back to

    def _fit_left_pane(self):
        """Widen the split so no control in the left panel is cut off."""
        try:
            self.left_frame.update_idletasks()
            # Measure the panels, not their container: the column is a canvas
            # now, and a canvas reports the width it was given rather than the
            # width its contents need.
            needed = (self.controls_frame.winfo_reqwidth()
                      + self.left_scroll.winfo_reqwidth() + 16)
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

    def _set_theme(self, name):
        """Switch palette at runtime and repaint everything that caches one.

        ttk styles are re-declared wholesale; the plain-Tk widgets (canvas,
        listbox, text) bake their colours in at creation, so they are walked
        and re-applied rather than rebuilt.
        """
        self.theme_name = apply_theme(name)
        self.theme_var.set(self.theme_name)
        self._configure_styles()
        self._restyle_widgets(self.root)
        self.results_scaler.apply()      # re-asserts the table style
        self.legend.redraw()             # swatches are the new palette's
        self._draw_now()
        self._save_settings()

    def _make_collapsible(self, pane, horizontal, key, default):
        """Make a sash collapse and restore the panes it divides.

        Dragging always worked, but there was no way to push a panel fully out
        of the way and no way to get it back, so every pane that needed
        closing grew a button. Double-clicking the divider closes the panel it
        is pointing at; double-clicking again restores it to the size it had.
        The edge is the control, which leaves Expand as the only button that
        moves a panel - and that one earns its place by being a shortcut for
        the common case rather than the only way in.
        """
        pane.bind(
            "<Double-Button-1>",
            lambda event: self._toggle_pane(event, pane, horizontal, key,
                                            default))
        Tooltip(pane, "Drag to resize  ·  double-click to collapse or restore")

    def _toggle_pane(self, event, pane, horizontal, key, default):
        try:
            current = pane.sashpos(0)
        except tk.TclError:
            return None
        # A double-click anywhere inside a child lands on the child, but be
        # explicit: only a hit near the divider itself is this gesture.
        offset = event.x if horizontal else event.y
        if abs(offset - current) > SASH_HIT_SLOP:
            return None

        extent = (pane.winfo_width() if horizontal
                  else pane.winfo_height()) or 1
        if current <= SASH_EDGE_ZONE or current >= extent - SASH_EDGE_ZONE:
            pane.sashpos(0, self._sash_memory.get(key) or default)
        else:
            self._sash_memory[key] = current
            # Collapse toward the nearer end, so the gesture always closes the
            # smaller panel rather than whichever one the code happened to
            # list first.
            pane.sashpos(0, 0 if current < extent / 2 else extent)
        self._redraw_tree()
        return "break"

    def _toggle_legend(self):
        """Show or hide the colour key under the tree, and remember which."""
        if self.show_legend_var.get():
            self.legend.show()
        else:
            self.legend.hide()
        self._save_settings()

    def _restyle_widgets(self, widget):
        """Re-apply the palette to every plain-Tk widget beneath `widget`."""
        for child in widget.winfo_children():
            if isinstance(child, scrolledtext.ScrolledText) or \
                    isinstance(child, tk.Text):
                self._style_text_widget(child)
            elif isinstance(child, tk.Listbox):
                child.configure(
                    bg=COLORS["field_bg"], fg=COLORS["field_fg"],
                    selectbackground=COLORS["highlight"],
                    selectforeground=COLORS["button_fg"],
                    highlightbackground=COLORS["panel_edge"],
                    highlightcolor=COLORS["panel_edge"])
            elif isinstance(child, tk.Canvas):
                # The controls column is a canvas too, but it is chrome, not
                # a drawing surface - given the tree's background it turned
                # white behind the panels in every dark theme.
                child.configure(
                    bg=COLORS["app_bg"]
                    if child is getattr(self, "left_canvas", None)
                    else COLORS["bg"])
            elif isinstance(child, tk.Toplevel):
                child.configure(bg=COLORS["app_bg"])
            self._restyle_widgets(child)

    # ---- Settings ----

    def _make_checker(self):
        return MonophylyChecker(self.db,
                                min_branch_length=self.min_branch_length,
                                min_support=self.min_support)

    def _edit_weak_branch_settings(self):
        """Let the user say which branches are too weak to count."""
        dialog = tk.Toplevel(self.root)
        dialog.title("Weak Branch Settings")
        dialog.transient(self.root)
        dialog.resizable(False, False)
        frame = ttk.Frame(dialog, padding=12)
        frame.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frame, wraplength=380, justify=tk.LEFT, text=(
            "A group split only by weak branches is reported as "
            "'weak branches only' instead of broken: collapse those "
            "branches into polytomies and the group fits the tree. Leave a "
            "box empty to switch that test off.")).grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 10))

        length_var = tk.StringVar(value="" if self.min_branch_length is None
                                  else f"{self.min_branch_length:g}")
        support_var = tk.StringVar(value="" if self.min_support is None
                                   else f"{self.min_support:g}")
        ttk.Label(frame, text="Branch length at or below:").grid(
            row=1, column=0, sticky="w")
        ttk.Entry(frame, textvariable=length_var, width=12).grid(
            row=1, column=1, sticky="w", padx=6)
        ttk.Label(frame, text="Support below (tree's own scale,\n"
                              "e.g. 70 bootstrap or 0.95 posterior):").grid(
            row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(frame, textvariable=support_var, width=12).grid(
            row=2, column=1, sticky="w", padx=6, pady=(6, 0))

        def parse(text, label):
            text = text.strip()
            if not text:
                return None
            try:
                return float(text)
            except ValueError:
                raise ValueError(f"{label} must be a number or empty.")

        def apply():
            try:
                length = parse(length_var.get(), "Branch length")
                support = parse(support_var.get(), "Support")
            except ValueError as e:
                messagebox.showerror("Weak Branch Settings", str(e),
                                     parent=dialog)
                return
            self.min_branch_length, self.min_support = length, support
            self._save_settings()
            dialog.destroy()
            if self.results:
                self._update_status(
                    "Weak-branch settings changed - run the analysis again "
                    "to apply them.")

        buttons = ttk.Frame(frame)
        buttons.grid(row=3, column=0, columnspan=2, sticky="e", pady=(12, 0))
        ttk.Button(buttons, text="Cancel", command=dialog.destroy).pack(
            side=tk.RIGHT)
        ttk.Button(buttons, text="Apply", command=apply).pack(
            side=tk.RIGHT, padx=6)
        dialog.grab_set()

    def _settings_path(self):
        return os.path.join(self.db.cache_dir, "settings.json")

    def _load_settings(self):
        """Remember the theme between runs; a colour choice is not a per-run
        decision and re-picking it every launch would be tedious."""
        try:
            with open(self._settings_path(), encoding="utf-8") as handle:
                return json.load(handle)
        except (OSError, ValueError):
            return {}

    def _save_settings(self):
        try:
            with open(self._settings_path(), "w", encoding="utf-8") as handle:
                json.dump({"theme": self.theme_name,
                           "legend": self.show_legend_var.get(),
                           "min_branch_length": self.min_branch_length,
                           "min_support": self.min_support}, handle)
        except OSError:
            pass          # a read-only cache dir is not worth an error dialog

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
            sashthickness=SASH_THICKNESS,
        )
        # The divider has to look like something you can grab. It was a 6px
        # strip of the window background - indistinguishable from padding -
        # which is why panels needed buttons to resize them at all. Given a
        # contrasting fill and a row of grip dots it reads as a handle, and
        # the double-click-to-collapse gesture has something to sit on.
        style.configure(
            "Sash",
            background=COLORS["panel_edge"],
            lightcolor=COLORS["panel_edge"],
            bordercolor=COLORS["app_bg"],
            sashthickness=SASH_THICKNESS,
            gripcount=14,
            handlesize=0,
        )
        # One button surface for the whole app. button_bg is a near-white chip
        # in every dark theme, so painting every control with it turned the
        # buttons into the brightest objects on screen and left the content
        # they act on looking like the background. Buttons take the panel
        # surface and the panel's own text colour, and read as buttons through
        # a border and a hover state instead of through sheer brightness.
        # Hardcoded greys are gone too: they were invisible in half the themes.
        chrome = {
            "background": COLORS["panel_bg"],
            "foreground": COLORS["text"],
            "bordercolor": COLORS["panel_edge"],
            "lightcolor": COLORS["panel_edge"],
            "darkcolor": COLORS["panel_edge"],
            "focuscolor": "none",
        }
        chrome_map = {
            "background": [
                ("pressed", COLORS["highlight"]),
                ("active", COLORS["tab_bg"]),
                ("disabled", COLORS["app_bg"]),
            ],
            "foreground": [
                ("pressed", COLORS["button_fg"]),
                ("disabled", COLORS["panel_edge"]),
            ],
        }
        style.configure("App.TButton", padding=(10, 6), **chrome)
        style.map("App.TButton", **chrome_map)
        # A stepper is not a general-purpose button. Zoom −/+ and text A−/A+
        # are one control with two ends, clicked over and over, so they get a
        # tight square footprint and sit flush against each other. Given the
        # App.TButton treatment they read as two more peers of Reset and Fit,
        # which is exactly the wrong grouping.
        style.configure(
            "Stepper.TButton",
            # Vertical padding tuned so a stepper stands exactly as tall as
            # the Quiet buttons beside it - a 2px mismatch on one row reads
            # as a misaligned control, not as a deliberate hierarchy.
            padding=(0, 5),
            anchor="center",
            font=("Helvetica", 12, "bold"),
            **chrome,
        )
        style.map("Stepper.TButton", **chrome_map)
        # Secondary actions sitting beside a stepper - Reset, Fit, Expand.
        # Flat until hovered, because a row where five controls shout equally
        # loudly gives the eye nowhere to land; the stepper is the control
        # being used, these are the ones being reached for occasionally.
        style.configure(
            "Quiet.TButton",
            background=COLORS["app_bg"],
            # Quiet is about the *fill*, not the label. Dimming the text to
            # muted_text as well left Reset and Fit almost unreadable on the
            # dark bars - a control you cannot find is not restrained, it is
            # broken. Full text colour, no chrome until you hover.
            foreground=COLORS["text"],
            bordercolor=COLORS["app_bg"],
            lightcolor=COLORS["app_bg"],
            darkcolor=COLORS["app_bg"],
            relief="flat",
            focuscolor="none",
            padding=(9, 4),
            # ttk's stock TButton carries width=-11, a minimum of eleven
            # characters. That is what stretched "Fit" to the width of
            # "Reset" and left the row a wall of identical slabs; 0 lets each
            # label size its own button.
            width=0,
        )
        style.map(
            "Quiet.TButton",
            background=[
                ("active", COLORS["button_bg"]),
                ("pressed", COLORS["button_bg"]),
            ],
            foreground=[
                ("active", COLORS["button_fg"]),
                ("pressed", COLORS["button_fg"]),
                ("disabled", COLORS["muted_text"]),
            ],
            bordercolor=[("active", COLORS["panel_edge"])],
            lightcolor=[("active", COLORS["button_bg"])],
            darkcolor=[("active", COLORS["button_bg"])],
        )
        # The value a stepper is reporting - "100%", "10 pt". State, not a
        # control, so it is muted and centred in a fixed width; left to size
        # itself it shoves the row sideways every time it crosses a digit.
        style.configure(
            "Readout.TLabel",
            background=COLORS["app_bg"],
            foreground=COLORS["muted_text"],
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

        # Scrollbars: clam draws the thumb with its own bevel colours, and
        # leaving lightcolor/darkcolor at the theme default produced a pale
        # slab that was the brightest object on a dark window - worse, right
        # next to the tree it was framing.
        for scrollbar_style in ("TScrollbar", "Vertical.TScrollbar",
                                "Horizontal.TScrollbar"):
            style.configure(
                scrollbar_style,
                background=COLORS["panel_edge"],
                troughcolor=COLORS["app_bg"],
                bordercolor=COLORS["app_bg"],
                lightcolor=COLORS["panel_edge"],
                darkcolor=COLORS["panel_edge"],
                arrowcolor=COLORS["muted_text"],
                gripcount=0,
                relief="flat",
            )
            style.map(
                scrollbar_style,
                background=[("active", COLORS["muted_text"]),
                            ("pressed", COLORS["muted_text"])],
                lightcolor=[("active", COLORS["muted_text"])],
                darkcolor=[("active", COLORS["muted_text"])],
                arrowcolor=[("active", COLORS["text"])],
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

        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Fit Tree to Window",
                              command=self._zoom_to_fit)
        view_menu.add_command(label="Expand Tree Panel",
                              command=self._toggle_tree_expanded)
        view_menu.add_separator()
        theme_menu = tk.Menu(view_menu, tearoff=0)
        for name in THEMES:
            theme_menu.add_radiobutton(
                label=name, value=name, variable=self.theme_var,
                command=lambda n=name: self._set_theme(n))
        view_menu.add_cascade(label="Colour Theme", menu=theme_menu)
        view_menu.add_separator()
        view_menu.add_command(label="Full Screen", accelerator="F11",
                              command=self._toggle_fullscreen)
        menubar.add_cascade(label="View", menu=view_menu)

        analysis_menu = tk.Menu(menubar, tearoff=0)
        analysis_menu.add_command(label="Weak Branch Settings...",
                                  command=self._edit_weak_branch_settings)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

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
        self.root.bind("<F11>", self._toggle_fullscreen)
        self.root.bind("<Escape>", lambda e: (self._fullscreen
                                              and self._toggle_fullscreen()))

        # Main layout: a controls-and-tree row above, results across the full
        # width below.
        #
        # Results used to be nested inside the right-hand column, which cost
        # it the width of the controls panel for its whole height - so the
        # tables were clipped mid-column while the bottom of the controls
        # column sat empty. The tables are the widest thing in the app and the
        # controls are the shortest, so the split that suits both is
        # horizontal at the top and full-bleed underneath.
        outer_pane = ttk.PanedWindow(self.root, orient=tk.VERTICAL)
        outer_pane.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        # Named right_pane still: it is the pane whose sash Expand drives, and
        # everything that reaches for it wants the tree/results divider.
        self.right_pane = outer_pane
        self._make_collapsible(outer_pane, horizontal=False, key="results",
                               default=420)

        top_row = ttk.Frame(outer_pane)
        outer_pane.add(top_row, weight=3)

        main_pane = ttk.PanedWindow(top_row, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True)
        self.main_pane = main_pane
        self._make_collapsible(main_pane, horizontal=True, key="controls",
                               default=320)

        # ---- Left panel: controls + file list ----
        # Scrollable, because this column no longer owns the window's full
        # height - it shares the top row with the tree. Its three panels plus
        # a long file list can exceed that, and the alternative to scrolling
        # is silently clipping whichever control is packed last. The scrollbar
        # only appears when the content actually overflows.
        left_frame = ttk.Frame(main_pane, width=320)
        main_pane.add(left_frame, weight=1)
        self.left_frame = left_frame

        self.left_canvas = tk.Canvas(left_frame, highlightthickness=0,
                                     bg=COLORS["app_bg"], width=300)
        self.left_scroll = ttk.Scrollbar(left_frame, orient=tk.VERTICAL,
                                         command=self.left_canvas.yview)
        self.left_canvas.configure(yscrollcommand=self.left_scroll.set)
        self.left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        controls = ttk.Frame(self.left_canvas)
        self._controls_window = self.left_canvas.create_window(
            (0, 0), window=controls, anchor="nw")
        controls.bind("<Configure>", lambda _e: self._sync_left_scroll())
        self.left_canvas.bind("<Configure>", lambda _e: self._sync_left_scroll())
        self.controls_frame = controls

        # Taxonomy DB status
        db_frame = self._panel(controls, "Taxonomy Database", "taxonomy")
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
        tree_frame = self._panel(controls, "Tree Files", "trees")
        tree_frame.pack(fill=tk.X, padx=5, pady=5)
        self.tree_panel = tree_frame

        tree_btn_frame = ttk.Frame(tree_frame)
        tree_btn_frame.pack(fill=tk.X, padx=5, pady=2)
        for column, (label, command) in enumerate((
                ("Add Tree(s)", self._open_tree),
                ("Clear", self._clear_trees))):
            ttk.Button(tree_btn_frame, text=label, command=command,
                       style="App.TButton").grid(
                           row=0, column=column, sticky="ew", padx=2)
            tree_btn_frame.grid_columnconfigure(column, weight=1)

        # The list is sized to its contents rather than given the whole
        # column. One tree used to reserve a panel deep enough for thirty and
        # left most of the left-hand side as an empty black box; the space it
        # was hoarding now goes to the panels below it and, past
        # TREE_LIST_MAX_ROWS, to a scrollbar instead of more emptiness.
        # A Treeview rather than a Listbox, because each row now carries a
        # tick as well as a name. Excluding a tree from a batch used to mean
        # clearing the list and re-adding everything except that one; now it
        # is one click, and Run Batch Analysis re-runs whatever is still
        # ticked.
        list_frame = ttk.Frame(tree_frame)
        list_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=2)
        self.tree_listbox = ttk.Treeview(
            list_frame, show="tree", selectmode="browse",
            height=TREE_LIST_MIN_ROWS)
        self.tree_list_scroll = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self.tree_listbox.yview)
        self.tree_listbox.configure(yscrollcommand=self.tree_list_scroll.set)
        self.tree_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.tree_listbox.column("#0", width=200, stretch=True)
        self.tree_listbox.bind("<<TreeviewSelect>>", self._on_tree_select)
        # Button-1 before the selection handler: a click on the tick toggles
        # it and stops there, so ticking a tree does not also load it.
        self.tree_listbox.bind("<Button-1>", self._on_tree_click, add="")
        self.tree_listbox.bind("<space>", self._toggle_selected_file)

        # Analysis controls.
        #
        analysis_frame = self._panel(controls, "Analysis", "analysis")
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
            values=["all conflicts", "supported only", "weak branches only",
                    "intruded", "scattered", "both"])
        conflict_filter.pack(side=tk.LEFT, padx=4)
        conflict_filter.bind(
            "<<ComboboxSelected>>", lambda e: self._populate_group_tables())

        self.collapse_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            analysis_frame, text="Merge groups with identical tips",
            variable=self.collapse_var,
            command=self._populate_group_tables).pack(
                anchor="w", padx=5, pady=(0, 4))

        # Search lives in the Results panel, not here - see _build_ui below.
        self.search_var = tk.StringVar()

        # ---- Tree view: fills the rest of the top row ----
        tree_view_frame = self._panel(main_pane, "Tree Visualization",
                                     "tree_view")
        main_pane.add(tree_view_frame, weight=3)

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
        # Discrete steps rather than a slider. A slider fires once per pixel
        # of drag, and debouncing that meant the tree did not move until the
        # drag stopped - which reads as lag. One click, one redraw.
        zoom_frame = ttk.Frame(tree_view_frame)
        zoom_frame.pack(fill=tk.X, padx=5, pady=2)
        self.zoom_var = tk.DoubleVar(value=ZOOM_DEFAULT)
        ttk.Label(zoom_frame, text="Zoom:").pack(side=tk.LEFT)
        self.zoom_out_button = ttk.Button(
            zoom_frame, text="−", width=2, style="Stepper.TButton",
            command=lambda: self._step_zoom(-1))
        self.zoom_out_button.pack(side=tk.LEFT, padx=(6, 0))
        Tooltip(self.zoom_out_button, "Zoom out  (Ctrl/Cmd −)")
        self.zoom_in_button = ttk.Button(
            zoom_frame, text="+", width=2, style="Stepper.TButton",
            command=lambda: self._step_zoom(1))
        self.zoom_in_button.pack(side=tk.LEFT, padx=(1, 0))
        Tooltip(self.zoom_in_button, "Zoom in  (Ctrl/Cmd +)")
        self.zoom_label = ttk.Label(zoom_frame, text="100%", width=6,
                                    anchor="center", style="Readout.TLabel")
        self.zoom_label.pack(side=tk.LEFT, padx=(8, 6))
        ttk.Button(zoom_frame, text="Reset", style="Quiet.TButton",
                   command=self._reset_zoom).pack(side=tk.LEFT)
        ttk.Button(zoom_frame, text="Fit", style="Quiet.TButton",
                   command=self._zoom_to_fit).pack(side=tk.LEFT, padx=(2, 0))
        # Expand resizes the pane, not the tree, so it sits apart from the
        # zoom cluster rather than reading as one more zoom control.
        self.expand_button = ttk.Button(
            zoom_frame, text="Expand ⤢", style="Quiet.TButton",
            command=self._toggle_tree_expanded)
        self.expand_button.pack(side=tk.RIGHT, padx=(10, 2))

        # Display toggles get their own row. Sharing one with the zoom
        # controls meant the row overflowed once Fit and Expand were added,
        # and "Labels" was silently clipped off the right edge - which looks
        # exactly like the labels feature being broken.
        display_frame = ttk.Frame(tree_view_frame)
        display_frame.pack(fill=tk.X, padx=5, pady=(0, 3))
        ttk.Label(display_frame, text="Show:").pack(side=tk.LEFT)

        self.show_labels_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Labels",
                        variable=self.show_labels_var,
                        command=self._redraw_tree).pack(side=tk.LEFT, padx=6)

        self.show_support_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(display_frame, text="Support",
                        variable=self.show_support_var,
                        command=self._redraw_tree).pack(side=tk.LEFT, padx=6)

        self.show_branch_lengths_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(display_frame, text="Branch lengths",
                        variable=self.show_branch_lengths_var,
                        command=self._redraw_tree).pack(side=tk.LEFT, padx=6)

        ttk.Checkbutton(display_frame, text="Legend",
                        variable=self.show_legend_var,
                        command=self._toggle_legend).pack(side=tk.LEFT, padx=6)

        ttk.Label(display_frame, text="  View:").pack(side=tk.LEFT,
                                                     padx=(12, 0))
        self.view_mode_var = tk.StringVar(value="Rectangular")
        view_box = ttk.Combobox(
            display_frame, textvariable=self.view_mode_var, state="readonly",
            width=12, values=["Rectangular", "Circular"])
        view_box.pack(side=tk.LEFT, padx=4)
        view_box.bind("<<ComboboxSelected>>", lambda e: self._draw_now())

        self.label_note = ttk.Label(display_frame, text="",
                                    foreground=COLORS["muted_text"])
        self.label_note.pack(side=tk.RIGHT)

        # The key goes beside the tree, not only in Help. Which hue means
        # conflict is a property of the current theme, and a reader looking at
        # an orange branch should not have to leave the tree to find out what
        # orange means. Built last so it can be packed into its slot above the
        # zoom row, and re-packed there when toggled back on.
        self.legend = LegendBar(tree_view_frame, fill=tk.X, padx=5,
                                pady=(0, 2), before=zoom_frame)
        if self.show_legend_var.get():
            self.legend.show()

        # Results panel
        results_frame = self._panel(outer_pane, "Results", "results")
        outer_pane.add(results_frame, weight=2)

        # Text size applies to the whole results panel - the two text tabs
        # and the row/heading font of the three tables - so switching tabs
        # never lands on a control that does nothing.
        results_bar = self.results_scaler.build_controls(
            results_frame, hint="(Ctrl/Cmd +, −, 0)")
        results_bar.pack(fill=tk.X, padx=5, pady=(4, 0))

        # Find rides along on the Results bar rather than owning a panel in
        # the left column. It searches what the results describe and marks it
        # in the tree, so it belongs next to them - and the row it joins had
        # a wide empty middle, so the whole panel's worth of height came back
        # for free.
        ttk.Separator(results_bar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=(14, 10), pady=2)
        ttk.Label(results_bar, text="Find:").pack(side=tk.LEFT)
        search_entry = ttk.Entry(results_bar, textvariable=self.search_var,
                                 width=22)
        search_entry.pack(side=tk.LEFT, padx=(6, 0))
        search_entry.bind("<Return>", lambda e: self._search_tips())
        self.search_var.trace_add("write", lambda *a: self._search_tips())
        Tooltip(search_entry, PANEL_HELP["search"][0])
        self.search_note = ttk.Label(results_bar, text="",
                                     style="Readout.TLabel")
        self.search_note.pack(side=tk.LEFT, padx=(8, 0))

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
                     ("support", "Support", 110),
                     ("culprit", "One tip explains it", 230),
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
        self._update_zoom_controls()

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
                ("NCBI archives", ARCHIVE_PATTERNS),
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
                ("Tree archives", ARCHIVE_PATTERNS),
                ("All files", "*.*"),
            ]
        )
        if not archive:
            return
        if not is_archive(archive):
            # A single compressed tree such as tree.tre.gz.
            self._add_tree_files([archive], select_last=True)
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

    def _resize_tree_list(self):
        """Fit the file list to how many files there actually are.

        Under the cap the panel takes only the depth it needs and the spare
        height goes to the panels below it. At the cap it stops growing and a
        scrollbar appears, which is a better use of the space than a taller
        box that is still mostly empty.
        """
        count = len(self.tree_files)
        rows = max(TREE_LIST_MIN_ROWS, min(count, TREE_LIST_MAX_ROWS))
        self.tree_listbox.configure(height=rows)

        if count > rows and not self.tree_list_scroll.winfo_manager():
            self.tree_list_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        elif count <= rows and self.tree_list_scroll.winfo_manager():
            self.tree_list_scroll.pack_forget()
        # The panel always takes exactly its content's height: the column
        # around it scrolls now, so there is no spare height here to claim and
        # nothing below it to squeeze out.
        self._sync_left_scroll()

    def _row_text(self, filepath):
        tick = TICK_ON if self.file_enabled.get(filepath, True) else TICK_OFF
        return f" {tick}  {os.path.basename(filepath)}"

    def _enabled_files(self):
        """The trees a batch run should cover, in list order."""
        return [fp for fp in self.tree_files if self.file_enabled.get(fp, True)]

    def _on_tree_click(self, event):
        """Toggle a row's tick when the click lands on it."""
        row = self.tree_listbox.identify_row(event.y)
        if not row or event.x > TICK_HIT_WIDTH:
            return None         # a normal click: let it select and load
        self._toggle_file(row)
        return "break"

    def _toggle_selected_file(self, _event=None):
        for row in self.tree_listbox.selection():
            self._toggle_file(row)
        return "break"

    def _toggle_file(self, filepath):
        """Include or exclude one tree without disturbing the rest of the list."""
        if filepath not in self.file_enabled:
            return
        self.file_enabled[filepath] = not self.file_enabled[filepath]
        self.tree_listbox.item(filepath, text=self._row_text(filepath))
        enabled = len(self._enabled_files())
        self._update_status(
            f"{os.path.basename(filepath)} "
            f"{'included' if self.file_enabled[filepath] else 'excluded'}"
            f" - {enabled} of {_plural(len(self.tree_files), 'tree')} "
            f"in the next batch run")

    def _add_tree_files(self, filepaths, select_last=False):
        """Add tree paths to the file list, skipping duplicates."""
        for fp in filepaths:
            if fp not in self.tree_files:
                self.tree_files.append(fp)
                self.file_enabled[fp] = True
                # The path is the row id, so every later lookup is a dict hit
                # rather than an index that shifts when the list changes.
                self.tree_listbox.insert("", tk.END, iid=fp,
                                         text=self._row_text(fp))
        self._resize_tree_list()

        if select_last and self.tree_files:
            last = self.tree_files[-1]
            self.tree_listbox.selection_set(last)
            self.tree_listbox.see(last)
            self._load_tree(last)

    def _clear_trees(self):
        self.tree_files.clear()
        self.file_enabled.clear()
        self.tree_listbox.delete(*self.tree_listbox.get_children())
        self._resize_tree_list()
        self.current_tree = None
        self.batch_results.clear()
        self.tree_canvas.delete("all")
        self._clear_analysis_state()
        self._clear_results()

    def _on_tree_select(self, event):
        selection = self.tree_listbox.selection()
        if selection:
            self._load_tree(selection[0])

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
                checker = self._make_checker()
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
        self.search_hits = set()

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
        # Size the right-hand gutter to the longest tip name actually drawn.
        # A flat 320*zoom reserved ~900px of empty canvas at high zoom, which
        # is pure scrolling distance for no content.
        font_size = max(6, int(9 * zoom))
        metrics = tkfont.Font(family="Helvetica", size=font_size)
        longest = max((metrics.measure(tip.label) for tip in tips), default=0)
        width = (max(x for x, _ in positions.values())
                 + longest + LABEL_GUTTER_PAD)
        height = MARGIN_TOP + max(1, len(tips)) * tip_spacing + MARGIN_TOP
        return positions, width, height

    def _layout_circular(self, zoom):
        """Polar layout: tips fanned around a circle, depth as radius.

        Worth having for its own sake, but it also answers the large-tree
        problem directly - a circle puts n tips around a circumference rather
        than down a column, so a tree that needs 10,000px of vertical
        scrolling as a phylogram fits on one screen as a fan.

        A wedge is left open so the root is visible and the reading order is
        obvious, which is why the sweep is not a full turn.
        """
        tips = self.current_tree.get_tip_nodes()
        count = max(1, len(tips))

        depth_of, dist_of = {}, {}
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
        use_lengths = self.show_branch_lengths_var.get()
        outer = CIRCLE_RADIUS * zoom

        radius_of = {}
        for key in depth_of:
            fraction = (dist_of[key] / max_dist if use_lengths
                        else depth_of[key] / max_depth)
            radius_of[key] = CIRCLE_INNER + fraction * (outer - CIRCLE_INNER)

        angle_of = {}
        step = CIRCLE_SWEEP / count
        for index, tip in enumerate(tips):
            angle_of[id(tip)] = CIRCLE_START + index * step
        for node in self.current_tree.post_order():
            if node.children:
                angles = [angle_of[id(child)] for child in node.children]
                angle_of[id(node)] = (min(angles) + max(angles)) / 2.0

        # Room for the tip names, which radiate outward.
        font_size = max(6, int(9 * zoom))
        metrics = tkfont.Font(family="Helvetica", size=font_size)
        longest = max((metrics.measure(tip.label) for tip in tips), default=0)
        margin = outer + longest + 30
        centre = (margin, margin)

        positions = {}
        for key, radius in radius_of.items():
            radians = math.radians(angle_of[key])
            positions[key] = (centre[0] + radius * math.cos(radians),
                              centre[1] - radius * math.sin(radians))
        return positions, angle_of, radius_of, centre, margin * 2

    def _draw_circular(self, canvas, zoom, space):
        """Render the polar layout: an arc per node, a spoke per child."""
        positions, angle_of, radius_of, centre, size = self._layout_circular(zoom)
        cx, cy = centre
        highlight = self.highlight_tips
        line_width = max(0.7, 1.2 * zoom)

        # Band behind the selected clade, drawn first so branches sit on top.
        # The rectangular view has always had one; without its equivalent here
        # a click in the results table moved the selection to somewhere on a
        # circle with nothing to say where, which is the whole reason picking
        # a clade felt ambiguous in this view. A clade occupies a contiguous
        # run of angles, so its band is the wedge those angles sweep - carried
        # out past the tips so the names in the group are inside it too.
        if highlight and self._highlight_is_contiguous():
            tips = self.current_tree.get_tip_nodes()
            angles = [angle_of[id(t)] for t in tips if t.label in highlight]
            if angles:
                # Half a tip-slot of padding at each end, so the outermost
                # members sit inside the wedge rather than on its edge.
                pad = CIRCLE_SWEEP / max(1, len(tips)) / 2
                reach = cx * 0.98        # centre to the edge of the label ring
                canvas.create_arc(
                    cx - reach, cy - reach, cx + reach, cy + reach,
                    start=min(angles) - pad,
                    extent=max(max(angles) - min(angles) + 2 * pad, 0.75),
                    style=tk.PIESLICE, fill=COLORS["clade_band"], outline="")
                # Punch the hub back out. A pie slice runs to the centre, and
                # the disc inside the first ring is the root's, not the
                # selected clade's - leaving it filled says the selection
                # reaches all the way down, which it does not.
                canvas.create_oval(
                    cx - CIRCLE_INNER, cy - CIRCLE_INNER,
                    cx + CIRCLE_INNER, cy + CIRCLE_INNER,
                    fill=COLORS["bg"], outline="")

        for node in self.current_tree.pre_order():
            if not node.children:
                continue
            radius = radius_of[id(node)]
            angles = [angle_of[id(child)] for child in node.children]
            low, high = min(angles), max(angles)
            colour = self._branch_color(node)
            if high > low:
                canvas.create_arc(
                    cx - radius, cy - radius, cx + radius, cy + radius,
                    start=low, extent=high - low, style=tk.ARC,
                    outline=colour, width=line_width,
                    **self._branch_dash(node))
            for child in node.children:
                theta = math.radians(angle_of[id(child)])
                child_radius = radius_of[id(child)]
                emphasis = highlight and child.label in highlight
                canvas.create_line(
                    cx + radius * math.cos(theta),
                    cy - radius * math.sin(theta),
                    cx + child_radius * math.cos(theta),
                    cy - child_radius * math.sin(theta),
                    fill=self._branch_color(child),
                    width=(2.4 if emphasis else 1.0) * max(0.7, zoom),
                    **self._branch_dash(child))

        if not self.show_labels_var.get():
            self.label_note.configure(text="labels off")
            self._hidden_labels = 0
            # Marks still go on with labels off - that is the case where a
            # ring is the only thing telling you where a hit is.
            self._draw_search_marks(canvas, positions, zoom)
            canvas.configure(scrollregion=(0, 0, size, size))
            return

        font_size = max(6, int(9 * zoom))
        tips = self.current_tree.get_tip_nodes()
        spacing = (2 * math.pi * CIRCLE_RADIUS * zoom
                   * (CIRCLE_SWEEP / 360.0) / max(1, len(tips)))
        room = spacing >= 7
        if room:
            for tip in tips:
                angle = angle_of[id(tip)]
                theta = math.radians(angle)
                radius = radius_of[id(tip)] + 6
                match = self.tip_matches.get(tip.label)
                if tip.label in self.highlight_intruders:
                    colour, weight = COLORS["non_mono"], "bold"
                elif match is None:
                    colour, weight = COLORS["unresolved"], "normal"
                elif not match.is_exact:
                    colour, weight = COLORS["inexact"], "normal"
                else:
                    colour, weight = COLORS["tip"], "normal"
                # Flip the text on the left half so no name reads upside down.
                flipped = 90 < (angle % 360) < 270
                item = canvas.create_text(
                    cx + radius * math.cos(theta),
                    cy - radius * math.sin(theta),
                    text=tip.label, anchor="e" if flipped else "w",
                    angle=(angle + 180) % 360 if flipped else angle % 360,
                    fill=colour, font=("Helvetica", font_size, weight))
                # Deliberately not reserved. Tk reports the axis-aligned box
                # of rotated text, which is far larger than the glyphs and,
                # around a circle, overlaps its neighbours even when the
                # names visibly do not. Reserving those boxes starved the
                # clade names of anywhere to sit.

        # Precomputed in one pass: calling count_tips() per node inside a
        # sort key is O(tips) each time and turns this into a quadratic.
        tip_counts = {}
        for node in self.current_tree.post_order():
            tip_counts[id(node)] = (1 if node.is_tip else
                                    sum(tip_counts[id(c)]
                                        for c in node.children))

        hidden = 0
        labelled = sorted((n for n in self.current_tree.pre_order()
                           if not n.is_tip and n.label),
                          key=lambda n: -tip_counts[id(n)])
        for node in labelled:
            # Push the name a little outward along its own spoke. Internal
            # nodes bunch up near the centre, so placing labels exactly on
            # the node loses most of them to collisions.
            theta = math.radians(angle_of[id(node)])
            radius = radius_of[id(node)] + CIRCLE_LABEL_OFFSET
            x = cx + radius * math.cos(theta)
            y = cy - radius * math.sin(theta)
            colour = self._branch_color(node)
            if colour == COLORS["branch"]:
                colour = COLORS["unresolved"]
            item = canvas.create_text(
                x, y, text=node.label, anchor="center", fill=colour,
                font=("Helvetica", max(9, font_size - 1), "bold"))
            rect = space.pad(canvas.bbox(item))
            if space.collides(rect):
                canvas.delete(item)
                hidden += 1
                continue
            chip = canvas.create_rectangle(
                *rect, fill=COLORS["bg"], outline=COLORS["chip_edge"])
            canvas.tag_lower(chip, item)
            space.add(rect)

        self._hidden_labels = hidden
        self._describe_hidden_labels(hidden, room, spacing)
        self._draw_search_marks(canvas, positions, zoom)
        canvas.configure(scrollregion=(0, 0, size, size))

    def _highlight_is_contiguous(self):
        """True when the highlighted tips form one unbroken run.

        The band behind a selection only means anything for a clade, which is
        always contiguous in tip order. A search can match tips scattered all
        over the tree, and banding min-to-max there paints over most of the
        drawing and claims everything between the hits is selected.
        """
        highlight = self.highlight_tips
        if not highlight:
            return False
        indices = [i for i, tip in enumerate(self.current_tree.get_tip_nodes())
                   if tip.label in highlight]
        return bool(indices) and \
            indices[-1] - indices[0] + 1 == len(indices)

    def _draw_search_marks(self, canvas, positions, zoom):
        """Ring every node the search matched, so a hit is findable by eye.

        Colour alone cannot do this: a matched tip may already be coloured by
        match quality or by conflict, and a matched clade is a name on a chip
        somewhere in the drawing. A ring is a mark nothing else uses.
        """
        radius = max(4.0, 6.0 * zoom)
        for node in self.current_tree.pre_order():
            if id(node) not in self.search_hits:
                continue
            spot = positions.get(id(node))
            if spot is None:
                continue
            x, y = spot
            canvas.create_oval(x - radius, y - radius, x + radius, y + radius,
                               outline=COLORS["highlight"], width=2)

    def _branch_dash(self, node):
        """Dash pattern for a branch, when the theme encodes status without hue."""
        pattern = COLORS.get("conflict_dash")
        if pattern and node.monophyly_status == "not_monophyletic":
            return {"dash": pattern}
        return {}

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
        if self.view_mode_var.get() == "Circular":
            self._draw_circular(canvas, zoom, LabelSpace())
            return

        positions, width, height = self._layout_tree(zoom)
        tip_spacing = max(4.0, TIP_SPACING * zoom)
        show_labels = self.show_labels_var.get()
        font_size = max(6, int(9 * zoom))
        # Clade names keep a legible size however far out you zoom - they are
        # the summary you zoom out in order to read. Only their number is
        # thinned, by the collision pass below.
        clade_font_size = max(9, font_size - 1)
        # Tip names are different: below this they physically cannot fit
        # between adjacent tips, and drawing them just paints the canvas solid.
        room_for_tip_labels = show_labels and tip_spacing >= 7

        highlight = self.highlight_tips
        intruders = self.highlight_intruders

        # Band behind the selected clade, drawn first so branches sit on top.
        if highlight and self._highlight_is_contiguous():
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
                style = self._branch_dash(child)
                canvas.create_line(x, y, x, cy, fill=color,
                                   width=line_width, **style)
                canvas.create_line(x, cy, cx, cy, fill=color,
                                   width=line_width, **style)

        # Everything textual is placed through here so labels stop landing on
        # top of one another. Order is by priority: tip names are the data and
        # always win, support values are pinned to a branch and cannot move,
        # clade chips are free to shift and so give way.
        space = LabelSpace()

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
                item = canvas.create_text(
                    x + 6, y, text=tip.label, anchor="w", fill=color,
                    font=("Helvetica", font_size, weight))
                space.add(space.pad(canvas.bbox(item)))

        if not show_labels:
            self._hidden_labels = 0
            self.label_note.configure(text="labels off")
            self._draw_search_marks(canvas, positions, zoom)
            canvas.configure(scrollregion=(0, 0, width, height))
            return

        # Support values, small and grey, above and left of their node.
        if self.show_support_var.get() and tip_spacing >= 9:
            for node in self.current_tree.pre_order():
                if node.is_tip or not node.support:
                    continue
                x, y = positions[id(node)]
                item = canvas.create_text(
                    x - 3, y - 6 * max(0.7, zoom), text=str(node.support),
                    anchor="e", fill=COLORS["support"],
                    font=("Helvetica", max(6, font_size - 2)))
                rect = space.pad(canvas.bbox(item), x=1, y=0)
                if space.collides(rect):
                    canvas.delete(item)   # crowded branch; the value is noise
                else:
                    space.add(rect)

        # Clade names, drawn to the RIGHT of the node on an opaque chip.
        # Drawing them to the left clipped every shallow clade off the canvas.
        #
        # Biggest clades are placed first, so when two nested groups compete
        # for the same spot the more inclusive one keeps it. A chip that
        # cannot sit on its node is nudged vertically and gets a leader line;
        # if nothing fits it is dropped rather than drawn over its neighbour.
        tip_counts = {}
        for node in self.current_tree.post_order():
            tip_counts[id(node)] = (1 if node.is_tip else
                                    sum(tip_counts[id(c)]
                                        for c in node.children))

        labelled = [n for n in self.current_tree.pre_order()
                    if not n.is_tip and n.label]
        labelled.sort(key=lambda n: -tip_counts[id(n)])

        step = max(5.0, tip_spacing * 0.55)
        offsets = [0.0]
        for multiple in (1, 2, 3, 4):
            offsets.extend((-step * multiple, step * multiple))

        hidden = 0
        for node in labelled:
            x, y = positions[id(node)]
            color = self._branch_color(node)
            if color == COLORS["branch"]:
                color = COLORS["unresolved"]

            item = canvas.create_text(
                x + 5, y, text=node.label, anchor="w", fill=color,
                font=("Helvetica", clade_font_size, "bold"))
            rect = space.pad(canvas.bbox(item))

            shift = next(
                (dy for dy in offsets
                 if not space.collides((rect[0], rect[1] + dy,
                                        rect[2], rect[3] + dy))),
                None)
            if shift is None:
                canvas.delete(item)
                hidden += 1
                continue

            if shift:
                canvas.move(item, 0, shift)
                rect = (rect[0], rect[1] + shift, rect[2], rect[3] + shift)
                canvas.create_line(
                    x, y, rect[0], (rect[1] + rect[3]) / 2,
                    fill=COLORS["chip_edge"], width=1)

            chip = canvas.create_rectangle(
                *rect, fill=COLORS["bg"], outline=COLORS["chip_edge"])
            canvas.tag_lower(chip, item)
            space.add(rect)

        self._hidden_labels = hidden
        self._describe_hidden_labels(hidden, room_for_tip_labels, tip_spacing)
        self._draw_search_marks(canvas, positions, zoom)
        canvas.configure(scrollregion=(0, 0, width, height))

    # ---- Zoom ----

    def _step_zoom(self, direction):
        """Move one stop along the zoom ladder, keeping the view centred."""
        current = self.zoom_var.get()
        index = min(range(len(ZOOM_LEVELS)),
                    key=lambda i: abs(ZOOM_LEVELS[i] - current))
        target = ZOOM_LEVELS[max(0, min(len(ZOOM_LEVELS) - 1,
                                        index + direction))]
        self._apply_zoom(target)

    def _reset_zoom(self):
        self._apply_zoom(ZOOM_DEFAULT)

    def _fit_zoom_for(self, tip_count):
        """Largest zoom stop at which `tip_count` tips fit the viewport.

        Tip spacing scales linearly with zoom, so this is arithmetic rather
        than a trial layout.
        """
        if tip_count <= 0:
            return ZOOM_DEFAULT
        viewport = self.tree_canvas.winfo_height() or 420
        usable = max(120, viewport - 2 * MARGIN_TOP)
        ideal = usable / (tip_count * TIP_SPACING)
        fitting = [level for level in ZOOM_LEVELS if level <= ideal]
        return fitting[-1] if fitting else ZOOM_LEVELS[0]

    def _zoom_to_fit(self):
        """Zoom so the whole tree is visible, or as close as the range allows."""
        if not self.current_tree:
            return
        tips = self.current_tree.count_tips()
        self._apply_zoom(self._fit_zoom_for(tips))
        self.tree_canvas.yview_moveto(0.0)
        self.tree_canvas.xview_moveto(0.0)
        self.tree_canvas.update_idletasks()

        needed = tips * TIP_SPACING * self.zoom_var.get()
        viewport = self.tree_canvas.winfo_height() or 1
        if needed > viewport:
            # Honest about it rather than silently leaving them scrolled.
            self._update_status(
                f"{tips} tips will not fit even at "
                f"{self.zoom_var.get() * 100:.0f}% - use Expand for a taller "
                f"panel, or scroll")
        else:
            self._update_status(f"Fitted {_plural(tips, 'tip')} in view")

    def _toggle_tree_expanded(self):
        """Give the tree the whole right-hand column, or share it again.

        The panes have always been draggable, but a sash is not a control
        anyone finds. This makes "let me actually look at the tree" a button.
        """
        self.tree_expanded = not self.tree_expanded
        try:
            if self.tree_expanded:
                self._results_sash = self.right_pane.sashpos(0)
                self.right_pane.sashpos(0, self.right_pane.winfo_height() - 60)
                self.expand_button.configure(text="Restore ⤡")
            else:
                self.right_pane.sashpos(0, self._results_sash or 420)
                self.expand_button.configure(text="Expand ⤢")
        except tk.TclError:
            pass
        self._redraw_tree()

    def _toggle_fullscreen(self, _event=None):
        self._fullscreen = not getattr(self, "_fullscreen", False)
        try:
            self.root.attributes("-fullscreen", self._fullscreen)
        except tk.TclError:
            # Not every window manager supports it; fall back to zoomed.
            try:
                self.root.state("zoomed" if self._fullscreen else "normal")
            except tk.TclError:
                pass
        return "break"

    def _apply_zoom(self, level):
        if level == self.zoom_var.get():
            return
        # Hold the middle of the viewport still, so zooming does not throw
        # the user to a different part of a 3000px-tall tree.
        try:
            top, bottom = self.tree_canvas.yview()
            anchor = (top + bottom) / 2
        except (tk.TclError, ValueError):
            anchor = None

        self.zoom_var.set(level)
        self._draw_now()

        if anchor is not None:
            view = self.tree_canvas.yview()
            self.tree_canvas.yview_moveto(
                max(0.0, anchor - (view[1] - view[0]) / 2))

    def _zoom_by_wheel(self, event):
        """Ctrl/Cmd + wheel over the canvas zooms; the bare wheel scrolls."""
        step = mousewheel_units(event)
        if step:
            self._step_zoom(-step)
        return "break"

    def _describe_hidden_labels(self, hidden, tips_drawn, tip_spacing):
        """Say what is not being shown, and why, next to the toggles."""
        notes = []
        if not tips_drawn:
            notes.append("tip names need more zoom")
        if hidden:
            notes.append(f"{_plural(hidden, 'clade name')} hidden")
        self.label_note.configure(text="  ·  ".join(notes))

    def _redraw_tree(self):
        """Coalesce rapid redraws from the display checkboxes."""
        if self._redraw_job is not None:
            self.root.after_cancel(self._redraw_job)
        self._redraw_job = self.root.after(40, self._do_redraw)

    def _draw_now(self):
        """Redraw immediately, cancelling anything already scheduled."""
        if self._redraw_job is not None:
            self.root.after_cancel(self._redraw_job)
        self._do_redraw()

    def _do_redraw(self):
        self._redraw_job = None
        self._draw_tree()
        self._update_zoom_controls()

    def _update_zoom_controls(self):
        level = self.zoom_var.get()
        self.zoom_label.configure(text=f"{level * 100:.0f}%")
        self.zoom_out_button.configure(
            state=tk.DISABLED if level <= ZOOM_LEVELS[0] else tk.NORMAL)
        self.zoom_in_button.configure(
            state=tk.DISABLED if level >= ZOOM_LEVELS[-1] else tk.NORMAL)
        # Say so when labels had to be dropped, rather than leaving the user
        # wondering why a clade they can see has no name on it.
        if self._hidden_labels:
            self._update_status(
                f"{_plural(self._hidden_labels, 'clade label')} hidden for "
                f"lack of room - zoom in to see them")

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

        # Read the selection HERE, on the main thread. Tk widgets are not
        # thread-safe and querying one from the worker is undefined behaviour
        # that happens to work most of the time.
        selection = self.tree_listbox.selection()
        selected_name = (os.path.basename(selection[0]) if selection else None)

        def analyze():
            checker = self._make_checker()
            analysis = checker.analyze(
                self.current_tree, progress_callback=self._update_status)
            checker.label_tree(self.current_tree, analysis["node_results"])
            self._apply_analysis(analysis)

            if selected_name:
                self.batch_results[selected_name] = analysis

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

        # Only the ticked trees. Untick one and press this again to re-run
        # without it - the point of the ticks is that excluding a tree costs
        # a click rather than rebuilding the whole list.
        targets = self._enabled_files()
        if not targets:
            messagebox.showwarning(
                "No Trees Included",
                "Every tree in the list is unticked. Tick at least one to "
                "include it in the batch run.")
            return

        selection = self.tree_listbox.selection()
        selected_path = selection[0] if selection else None

        def analyze_batch():
            checker = self._make_checker()
            total = len(targets)
            for idx, fp in enumerate(targets):
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

            skipped = len(self.tree_files) - total
            self._update_status(
                f"Batch analysis complete: {_plural(total, 'tree')}"
                + (f", {skipped} unticked and skipped" if skipped else ""))

            # Refresh current view, using the path captured before the
            # worker started rather than querying the listbox from here.
            if selected_path:
                self.root.after(0, lambda: self._load_tree(selected_path))

        self._run_threaded(analyze_batch, "Running batch analysis...")

    def _show_results(self):
        """Update the results panels."""
        self._clear_results()
        if not self.results:
            return

        checker = self._make_checker()
        self.summary_text.insert(tk.END, checker.get_summary(
            self.results, self.taxon_results, self.tip_matches))

        self._populate_group_tables()
        self._populate_tip_matching()

        if any(g.is_broken for g in self.groups):
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
            if wanted_conflict == "supported only":
                if not group.is_broken:
                    continue
            elif wanted_conflict == "weak branches only":
                if not group.is_unsupported_conflict:
                    continue
            elif (wanted_conflict != "all conflicts"
                    and group.conflict_type != wanted_conflict):
                continue
            listed = ", ".join(group.intruders[:12])
            if len(group.intruders) > 12:
                listed += f"  ... (+{len(group.intruders) - 12} more)"
            row = self.nonmono_tree.insert("", tk.END, values=(
                group.name, group.rank, group.member_count,
                group.conflict_cost, group.conflict_type,
                "weak branches only" if group.is_unsupported_conflict
                else "supported",
                group.culprit_summary(), len(group.intruders),
                len(group.outliers), listed))
            self._remember_row(self.nonmono_tree, row, group)
            self._add_constituent_rows(self.nonmono_tree, row, group, 10)
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

        broken = sum(1 for g in self.groups if g.is_broken)
        fixable = sum(len(t.repairs) for t in repairing)
        if repairing:
            self.conflict_banner.configure(
                text=f"Removing any one of these "
                     f"{_plural(len(repairing), 'tip')} repairs groups; "
                     f"together they account for {fixable} of {broken} "
                     f"supported conflicts. A conflict can also mean "
                     f"outdated taxonomy, hybridization or lineage sorting, "
                     f"not only a bad sample.")
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

        retired = sorted((m for m in self.tip_matches.values()
                          if m.method == "merged"), key=lambda m: m.label)
        if retired:
            text.insert(tk.END, f"--- Outdated taxon IDs ({len(retired)}) ---\n")
            text.insert(tk.END, "NCBI has merged these IDs into newer ones; "
                                "the current taxon was used.\n\n")
            for m in retired:
                text.insert(tk.END, f"  {m.label} -> {m.tax_id} "
                                    f"({m.matched_name})\n")
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
                if entry.is_unsupported_conflict:
                    status += ("  |  only weak branches contradict it, so "
                               "the data do not reject this group")
                # Say plainly when the named tip is a choice, not the answer.
                if entry.is_ambiguous:
                    status += ("  |  any ONE of these repairs it: "
                               + ", ".join(entry.repair_options))

        self._draw_tree()
        if isinstance(entry, ConflictingTip):
            self._scroll_to_highlight()
        else:
            self._frame_highlight()
        self._update_status(status)

    def _frame_highlight(self):
        """Bring the highlighted clade into view, zooming out if it will not fit.

        Selecting a 26-member group at 280% used to drop the user inside the
        clade: the highlight band filled the viewport and there was nothing to
        see. Zoom is only ever reduced here, never raised, so an explicit zoom
        choice is respected - use the Fit button to go the other way.
        """
        if not self.highlight_tips or not self.current_tree:
            return
        fitting = self._fit_zoom_for(len(self.highlight_tips))
        if fitting < self.zoom_var.get():
            self._apply_zoom(fitting)

        positions, _width, height = self._layout_tree(self.zoom_var.get())
        ys = [positions[id(n)][1] for n in self.current_tree.get_tip_nodes()
              if n.label in self.highlight_tips]
        if ys:
            self._scroll_canvas_to((min(ys) + max(ys)) / 2, height)

    # Kept for the tip search, which frames a single tip and must not rezoom.
    def _scroll_to_nodes(self, nodes):
        """Centre the viewport on specific nodes, in whichever layout shows.

        Both axes, deliberately. A search hit can be off-screen sideways as
        easily as vertically - the internal node carrying a clade name sits
        far to the left of the tips it spans - and scrolling only vertically
        was why a match could be reported and still not be anywhere on
        screen.
        """
        if not nodes:
            return
        zoom = self.zoom_var.get()
        if self.view_mode_var.get() == "Circular":
            positions, _angles, _radii, _centre, size = \
                self._layout_circular(zoom)
            width = height = size
        else:
            positions, width, height = self._layout_tree(zoom)

        points = [positions[id(node)] for node in nodes
                  if id(node) in positions]
        if not points:
            return
        canvas = self.tree_canvas
        for moveto, values, first, view, extent in (
                (canvas.xview_moveto, [x for x, _ in points], points[0][0],
                 canvas.winfo_width() or 1, width),
                (canvas.yview_moveto, [y for _, y in points], points[0][1],
                 canvas.winfo_height() or 1, height)):
            # Centring the bounding box only helps while the box fits. Nested
            # matches - Rosales, Rosaceae, Rosoideae - span most of the tree,
            # and centring that box framed the empty middle with every hit
            # off both edges. Wider than the viewport, fall back to the first
            # hit, which at least guarantees one of them is on screen.
            if max(values) - min(values) > view:
                target = first
            else:
                target = (min(values) + max(values)) / 2
            moveto(max(0.0, min(1.0, (target - view / 2) / max(extent, 1))))

    def _scroll_to_highlight(self):
        """Bring the selected clade into view, in whichever layout is showing.

        This used to lay the tree out rectangularly no matter what was on
        screen, so in circular view it scrolled to a y computed for a drawing
        that was not there - landing somewhere arbitrary and making a clade
        click look broken. A fan needs both axes moved, not just the vertical.
        """
        zoom = self.zoom_var.get()
        selected = self.highlight_tips
        if not selected:
            return

        if self.view_mode_var.get() == "Circular":
            positions, _angles, _radii, _centre, size = \
                self._layout_circular(zoom)
            points = [positions[id(n)]
                      for n in self.current_tree.get_tip_nodes()
                      if n.label in selected]
            if not points:
                return
            x_mid = (min(x for x, _ in points) + max(x for x, _ in points)) / 2
            y_mid = (min(y for _, y in points) + max(y for _, y in points)) / 2
            for axis, centre_of, extent in (
                    (self.tree_canvas.xview_moveto, x_mid,
                     self.tree_canvas.winfo_width() or 1),
                    (self.tree_canvas.yview_moveto, y_mid,
                     self.tree_canvas.winfo_height() or 1)):
                axis(max(0.0, min(1.0, (centre_of - extent / 2) / max(size, 1))))
            return

        positions, _width, height = self._layout_tree(zoom)
        ys = [positions[id(n)][1] for n in self.current_tree.get_tip_nodes()
              if n.label in selected]
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
            checker = self._make_checker()
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
            "Accepts plain or compressed Newick (.nwk, .tre, .gz),\n"
            "tips as species names or NCBI taxon IDs, and treats\n"
            "near-zero or low-support branches as unresolved.\n\n"
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

    def _bind_mousewheel_to_tree(self, widget, target):
        """Scroll `target` from anywhere inside `widget`.

        Labels and frames do not scroll themselves, so without this the wheel
        did nothing over most of the controls column - only over the gaps
        between panels, which reads as the wheel being broken.
        """
        if widget is self.tree_listbox or widget is self.tree_list_scroll:
            return          # the file list scrolls itself
        def on_scroll(event):
            step = mousewheel_units(event)
            if step:
                target.yview_scroll(step, "units")
            return "break"
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            widget.bind(seq, on_scroll, add="+")
        for child in widget.winfo_children():
            self._bind_mousewheel_to_tree(child, target)

    def _bind_scroll_wheels(self):
        """Attach mouse wheel scrolling to all scrollable panels."""
        self._bind_mousewheel(self.tree_listbox, yview=True)
        self._bind_mousewheel(self.tree_canvas, yview=True, xview=True)
        # The controls column and every static panel inside it, so the wheel
        # works wherever the pointer happens to be resting in that column.
        self._bind_mousewheel(self.left_canvas, yview=True)
        for panel in self.controls_frame.winfo_children():
            self._bind_mousewheel_to_tree(panel, self.left_canvas)
        for seq in ("<Control-MouseWheel>", "<Command-MouseWheel>"):
            self.tree_canvas.bind(seq, self._zoom_by_wheel, add="+")
        self._bind_mousewheel(self.summary_text, yview=True)
        self._bind_mousewheel(self.nonmono_tree, yview=True)
        self._bind_mousewheel(self.mono_tree, yview=True)
        self._bind_mousewheel(self.conflict_tree, yview=True)
        self._bind_mousewheel(self.unresolved_text, yview=True)

    def _search_tips(self):
        """Highlight everything matching the search box and scroll to it.

        Two things were wrong with this. It looked only at tips, so typing a
        group that is drawn right there on the canvas - "hypn" for Hypnanae -
        reported no match at all. And it scrolled using the rectangular
        layout whatever view was showing, so in the circular view it jumped
        to a coordinate belonging to a drawing that was not on screen. Both
        added up to a Find that did not visibly find anything.
        """
        query = self.search_var.get().strip().lower()
        self.search_hits = set()
        if not self.current_tree:
            return
        if not query:
            self.highlight_tips = set()
            self.search_note.configure(text="")
            self._redraw_tree()
            return

        tips, clades = [], []
        for node in self.current_tree.pre_order():
            if node.label and query in node.label.lower():
                (tips if node.is_tip else clades).append(node)

        self.search_hits = {id(n) for n in tips + clades}
        selected = {node.label for node in tips}
        for clade in clades:
            selected.update(tip.label for tip in clade.get_tip_nodes())
        self.highlight_tips = selected
        self._draw_tree()

        if not self.search_hits:
            self.search_note.configure(text="no match")
            self._update_status(f"Nothing in the tree matches '{query}'")
            return

        # Scroll to the matches themselves, not to the selection they imply.
        # Nested clades like Rosales span most of the tree, so centring on the
        # tips they cover lands between the hits rather than on one.
        self._scroll_to_nodes(tips + clades)
        found = [part for part in
                 (_plural(len(tips), "tip") if tips else "",
                  _plural(len(clades), "clade") if clades else "") if part]
        summary = " + ".join(found)
        self.search_note.configure(text=summary)
        self._update_status(f"{summary} matching '{query}'")

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
