"""
webapp.py - Web-based UI for PhyLabeler

A local web application that runs in the browser. Uses only Python stdlib
(http.server + json). No Flask, no npm, no external dependencies.

Launch with: python3 webapp.py
Then open http://localhost:8080 in your browser.
"""

import http.server
import json
import os
import sys
import threading
import urllib.parse
import webbrowser
import socketserver
import io

from taxonomy_db import TaxonomyDB
from tree_parser import parse_newick, parse_newick_file
from monophyly import MonophylyChecker

# ---- Global state ----
db = TaxonomyDB()
current_tree = None
current_results = None
current_unresolved = []
status_message = "Ready"

PORT = 8080


def get_html():
    """Return the complete single-page application HTML."""
    return HTML_PAGE


class PhyLabelerHandler(http.server.BaseHTTPRequestHandler):
    """HTTP request handler for the PhyLabeler web app."""

    def log_message(self, format, *args):
        """Suppress default logging noise."""
        pass

    def do_GET(self):
        path = urllib.parse.urlparse(self.path).path

        if path == "/" or path == "/index.html":
            self._send_html(get_html())
        elif path == "/api/status":
            self._send_json({"status": status_message, "db_loaded": db.loaded,
                             "db_taxa": len(db.code_to_parent) if db.loaded else 0})
        elif path == "/api/cache-info":
            self._send_json(db.get_cache_info())
        else:
            self._send_error(404, "Not found")

    def do_POST(self):
        global current_tree, current_results, current_unresolved, status_message

        path = urllib.parse.urlparse(self.path).path
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length).decode("utf-8") if content_length else ""

        try:
            if path == "/api/load-cache":
                status_message = "Loading taxonomy cache..."
                success = db.load()
                if success:
                    status_message = f"Loaded {len(db.code_to_parent):,} taxa from cache"
                    self._send_json({"success": True, "taxa": len(db.code_to_parent)})
                else:
                    status_message = "No cache found"
                    self._send_json({"success": False, "error": "No cache found. Load files or download from NCBI."})

            elif path == "/api/load-files":
                data = json.loads(body)
                names_path = data.get("names_file", "")
                nodes_path = data.get("nodes_file", "")
                if not names_path or not nodes_path:
                    self._send_json({"success": False, "error": "Both files required"})
                    return
                status_message = "Parsing taxonomy files..."
                db.load(names_file=names_path, nodes_file=nodes_path)
                status_message = f"Loaded {len(db.code_to_parent):,} taxa"
                self._send_json({"success": True, "taxa": len(db.code_to_parent)})

            elif path == "/api/download-taxonomy":
                status_message = "Downloading NCBI taxonomy (this may take a few minutes)..."
                def _dl():
                    global status_message
                    try:
                        db.download_and_load(progress_callback=lambda m: None)
                        status_message = f"Downloaded and cached {len(db.code_to_parent):,} taxa"
                    except Exception as e:
                        status_message = f"Download failed: {e}"
                t = threading.Thread(target=_dl, daemon=True)
                t.start()
                self._send_json({"success": True, "message": "Download started"})

            elif path == "/api/load-tree":
                data = json.loads(body)
                newick_str = data.get("newick", "")
                filename = data.get("filename", "uploaded tree")
                if not newick_str.strip():
                    self._send_json({"success": False, "error": "Empty tree data"})
                    return
                current_tree = parse_newick(newick_str)
                current_results = None
                current_unresolved = []
                tree_data = _tree_to_dict(current_tree)
                tips = current_tree.count_tips()
                status_message = f"Loaded: {filename} ({tips} tips)"
                self._send_json({"success": True, "tree": tree_data, "tips": tips,
                                 "internal": current_tree.count_internal()})

            elif path == "/api/load-tree-file":
                data = json.loads(body)
                filepath = data.get("filepath", "")
                if not os.path.exists(filepath):
                    self._send_json({"success": False, "error": f"File not found: {filepath}"})
                    return
                current_tree = parse_newick_file(filepath)
                current_results = None
                current_unresolved = []
                tree_data = _tree_to_dict(current_tree)
                tips = current_tree.count_tips()
                status_message = f"Loaded: {os.path.basename(filepath)} ({tips} tips)"
                self._send_json({"success": True, "tree": tree_data, "tips": tips,
                                 "internal": current_tree.count_internal()})

            elif path == "/api/analyze":
                if not db.loaded:
                    self._send_json({"success": False, "error": "Load taxonomy database first"})
                    return
                if not current_tree:
                    self._send_json({"success": False, "error": "Load a tree first"})
                    return

                status_message = "Running monophyly analysis..."
                checker = MonophylyChecker(db)
                current_results, current_unresolved = checker.check_tree(current_tree)
                checker.label_tree(current_tree, current_results)

                results_data = _results_to_dict(current_results, current_unresolved)
                tree_data = _tree_to_dict(current_tree)
                summary = checker.get_summary(current_results)

                mono = sum(1 for r in current_results if r.status == "monophyletic")
                nonmono = sum(1 for r in current_results if r.status == "not_monophyletic")
                status_message = f"Done: {mono} monophyletic, {nonmono} non-monophyletic"

                self._send_json({
                    "success": True, "tree": tree_data, "results": results_data,
                    "summary": summary
                })

            elif path == "/api/export-newick":
                if not current_tree:
                    self._send_json({"success": False, "error": "No tree loaded"})
                    return
                newick = current_tree.to_newick(branch_lengths=True) + ";"
                self._send_json({"success": True, "newick": newick})

            elif path == "/api/export-report":
                if not current_results:
                    self._send_json({"success": False, "error": "Run analysis first"})
                    return
                checker = MonophylyChecker(db)
                report = checker.get_summary(current_results)
                if current_unresolved:
                    report += "\n\n--- Unresolved Tips ---\n"
                    for tip in sorted(current_unresolved):
                        report += f"  {tip}\n"
                self._send_json({"success": True, "report": report})

            else:
                self._send_error(404, "Not found")

        except Exception as e:
            status_message = f"Error: {e}"
            self._send_json({"success": False, "error": str(e)})

    def _send_html(self, html):
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(html.encode("utf-8"))

    def _send_json(self, data):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps(data).encode("utf-8"))

    def _send_error(self, code, message):
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps({"error": message}).encode("utf-8"))


def _tree_to_dict(node):
    """Convert tree to JSON-serializable dict for the frontend."""
    d = {
        "label": node.label or "",
        "branch_length": node.branch_length,
        "is_tip": node.is_tip,
        "monophyly": node.monophyly_status,
        "tax_name": node.taxonomy_name or "",
    }
    if node.children:
        d["children"] = [_tree_to_dict(c) for c in node.children]
    return d


def _results_to_dict(results, unresolved):
    """Convert analysis results to JSON-serializable dict."""
    mono = []
    non_mono = []

    for r in results:
        item = {
            "name": r.mrca_name,
            "rank": r.mrca_rank,
            "tips": len(r.tip_labels),
            "status": r.status,
        }
        if r.status == "monophyletic":
            mono.append(item)
        elif r.status == "not_monophyletic":
            item["intruders"] = r.intruders[:20]
            item["total_intruders"] = len(r.intruders)
            non_mono.append(item)

    mono.sort(key=lambda x: x["tips"], reverse=True)
    non_mono.sort(key=lambda x: x["tips"], reverse=True)

    return {
        "monophyletic": mono,
        "non_monophyletic": non_mono,
        "unresolved_tips": sorted(unresolved),
        "total_mono": len(mono),
        "total_non_mono": len(non_mono),
        "total_unresolved": len(unresolved),
    }


# ---- HTML Single Page App ----

HTML_PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PhyLabeler</title>
<style>
:root {
  --mono: #27ae60;
  --non-mono: #e74c3c;
  --unresolved: #95a5a6;
  --primary: #2c3e50;
  --primary-light: #34495e;
  --accent: #3498db;
  --bg: #f5f6fa;
  --card-bg: #ffffff;
  --border: #dcdde1;
  --text: #2c3e50;
  --text-light: #636e72;
  --shadow: 0 2px 8px rgba(0,0,0,0.08);
}

* { margin: 0; padding: 0; box-sizing: border-box; }

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, sans-serif;
  background: var(--bg);
  color: var(--text);
  height: 100vh;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}

/* Header */
header {
  background: var(--primary);
  color: white;
  padding: 12px 24px;
  display: flex;
  align-items: center;
  gap: 20px;
  box-shadow: 0 2px 10px rgba(0,0,0,0.15);
  z-index: 100;
}
header h1 { font-size: 20px; font-weight: 600; letter-spacing: -0.5px; }
header h1 span { color: var(--accent); }
.header-status {
  margin-left: auto;
  font-size: 13px;
  opacity: 0.85;
  background: rgba(255,255,255,0.1);
  padding: 4px 12px;
  border-radius: 12px;
}

/* Main layout */
.main-container {
  display: flex;
  flex: 1;
  overflow: hidden;
}

/* Sidebar */
.sidebar {
  width: 300px;
  min-width: 300px;
  background: var(--card-bg);
  border-right: 1px solid var(--border);
  overflow-y: auto;
  padding: 16px;
  display: flex;
  flex-direction: column;
  gap: 16px;
}

.panel {
  border: 1px solid var(--border);
  border-radius: 8px;
  overflow: hidden;
}
.panel-header {
  background: var(--bg);
  padding: 10px 14px;
  font-size: 13px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  color: var(--text-light);
  border-bottom: 1px solid var(--border);
}
.panel-body { padding: 12px 14px; }

.status-dot {
  display: inline-block;
  width: 8px;
  height: 8px;
  border-radius: 50%;
  margin-right: 6px;
}
.status-dot.loaded { background: var(--mono); }
.status-dot.empty { background: var(--non-mono); }
.status-dot.loading { background: #f39c12; animation: pulse 1s infinite; }
@keyframes pulse { 50% { opacity: 0.4; } }

.db-info { font-size: 13px; margin-bottom: 10px; color: var(--text-light); }

.btn-group { display: flex; gap: 6px; flex-wrap: wrap; }

button {
  padding: 8px 14px;
  border: 1px solid var(--border);
  border-radius: 6px;
  background: var(--card-bg);
  color: var(--text);
  font-size: 13px;
  cursor: pointer;
  transition: all 0.15s;
  font-family: inherit;
}
button:hover { background: var(--bg); border-color: var(--accent); }
button:active { transform: scale(0.97); }

.btn-primary {
  background: var(--accent);
  color: white;
  border-color: var(--accent);
}
.btn-primary:hover { background: #2980b9; }

.btn-success {
  background: var(--mono);
  color: white;
  border-color: var(--mono);
}
.btn-success:hover { background: #219a52; }

.btn-full { width: 100%; }
.btn-sm { padding: 6px 10px; font-size: 12px; }

input[type="file"] { display: none; }

.file-path-input {
  width: 100%;
  padding: 8px 10px;
  border: 1px solid var(--border);
  border-radius: 6px;
  font-size: 13px;
  font-family: 'SF Mono', Monaco, monospace;
  margin-bottom: 8px;
}
.file-path-input:focus { outline: none; border-color: var(--accent); }

/* Content area */
.content {
  flex: 1;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}

/* Tabs */
.tab-bar {
  display: flex;
  background: var(--card-bg);
  border-bottom: 1px solid var(--border);
  padding: 0 16px;
}
.tab {
  padding: 12px 20px;
  font-size: 13px;
  font-weight: 500;
  cursor: pointer;
  border-bottom: 2px solid transparent;
  color: var(--text-light);
  transition: all 0.15s;
}
.tab:hover { color: var(--text); }
.tab.active {
  color: var(--accent);
  border-bottom-color: var(--accent);
}
.tab .badge {
  display: inline-block;
  padding: 1px 7px;
  border-radius: 10px;
  font-size: 11px;
  margin-left: 6px;
  font-weight: 600;
}
.badge-mono { background: #d5f5e3; color: var(--mono); }
.badge-nonmono { background: #fadbd8; color: var(--non-mono); }
.badge-unresolved { background: #eaecee; color: var(--unresolved); }

/* Tab content */
.tab-content { flex: 1; overflow: hidden; display: none; }
.tab-content.active { display: flex; flex-direction: column; }

/* Tree view */
#tree-container {
  flex: 1;
  overflow: auto;
  padding: 20px;
  background: white;
}
#tree-svg { min-width: 100%; }

.tree-tip { font: 12px 'SF Mono', Monaco, monospace; fill: var(--text); cursor: pointer; }
.tree-tip:hover { fill: var(--accent); font-weight: bold; }
.tree-label {
  font: bold 11px -apple-system, sans-serif;
  cursor: pointer;
}
.tree-label.mono { fill: var(--mono); }
.tree-label.non-mono { fill: var(--non-mono); }
.tree-label.unresolved { fill: var(--unresolved); }
.tree-branch { stroke: #999; stroke-width: 1.5; fill: none; }
.tree-branch.mono { stroke: var(--mono); stroke-width: 2; }
.tree-branch.non-mono { stroke: var(--non-mono); stroke-width: 2; }

/* Zoom controls */
.zoom-bar {
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 8px 16px;
  background: var(--card-bg);
  border-bottom: 1px solid var(--border);
  font-size: 13px;
}
.zoom-bar input[type="range"] { width: 150px; }
.zoom-bar label { display: flex; align-items: center; gap: 4px; cursor: pointer; }

/* Results table */
.results-table-container { flex: 1; overflow: auto; padding: 16px; }

table {
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}
th {
  background: var(--bg);
  padding: 10px 14px;
  text-align: left;
  font-weight: 600;
  border-bottom: 2px solid var(--border);
  position: sticky;
  top: 0;
}
td {
  padding: 8px 14px;
  border-bottom: 1px solid var(--border);
}
tr:hover td { background: #f8f9fa; }
.intruder-list {
  font-size: 12px;
  color: var(--text-light);
  max-width: 400px;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}

/* Summary */
.summary-container {
  flex: 1;
  overflow: auto;
  padding: 20px;
}
.summary-box {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 20px;
  white-space: pre-wrap;
  font-family: 'SF Mono', Monaco, monospace;
  font-size: 13px;
  line-height: 1.6;
}

/* Stats cards */
.stats-row {
  display: flex;
  gap: 12px;
  padding: 16px;
  border-bottom: 1px solid var(--border);
}
.stat-card {
  flex: 1;
  padding: 16px;
  border-radius: 8px;
  text-align: center;
}
.stat-card.mono-card { background: #d5f5e3; }
.stat-card.nonmono-card { background: #fadbd8; }
.stat-card.unresolved-card { background: #eaecee; }
.stat-number {
  font-size: 28px;
  font-weight: 700;
  line-height: 1;
}
.stat-label {
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  margin-top: 4px;
  opacity: 0.7;
}

/* Empty state */
.empty-state {
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100%;
  color: var(--text-light);
  text-align: center;
  padding: 40px;
}
.empty-state .icon { font-size: 48px; margin-bottom: 16px; opacity: 0.3; }
.empty-state h3 { margin-bottom: 8px; font-weight: 500; }
.empty-state p { font-size: 14px; max-width: 400px; line-height: 1.5; }

/* Tooltip */
.tooltip {
  position: absolute;
  background: var(--primary);
  color: white;
  padding: 8px 12px;
  border-radius: 6px;
  font-size: 12px;
  pointer-events: none;
  z-index: 1000;
  max-width: 300px;
  box-shadow: 0 4px 12px rgba(0,0,0,0.2);
  display: none;
}

/* Loading overlay */
.loading-overlay {
  position: fixed;
  inset: 0;
  background: rgba(255,255,255,0.85);
  display: none;
  align-items: center;
  justify-content: center;
  z-index: 999;
  flex-direction: column;
  gap: 16px;
}
.loading-overlay.active { display: flex; }
.spinner {
  width: 36px;
  height: 36px;
  border: 3px solid var(--border);
  border-top-color: var(--accent);
  border-radius: 50%;
  animation: spin 0.7s linear infinite;
}
@keyframes spin { to { transform: rotate(360deg); } }
.loading-text { font-size: 14px; color: var(--text-light); }
</style>
</head>
<body>

<header>
  <h1>Phy<span>Labeler</span></h1>
  <div class="header-status" id="header-status">Ready</div>
</header>

<div class="main-container">
  <!-- Sidebar -->
  <div class="sidebar">
    <!-- Taxonomy DB -->
    <div class="panel">
      <div class="panel-header">Taxonomy Database</div>
      <div class="panel-body">
        <div class="db-info" id="db-info">
          <span class="status-dot empty" id="db-dot"></span>
          <span id="db-status-text">Not loaded</span>
        </div>
        <div class="btn-group">
          <button onclick="loadCache()" class="btn-sm">Load Cache</button>
          <button onclick="downloadNCBI()" class="btn-sm">Download NCBI</button>
        </div>
        <div style="margin-top:10px">
          <input type="text" class="file-path-input" id="names-path" placeholder="Path to names.dmp">
          <input type="text" class="file-path-input" id="nodes-path" placeholder="Path to nodes.dmp">
          <button onclick="loadFiles()" class="btn-sm btn-full">Load From Files</button>
        </div>
      </div>
    </div>

    <!-- Tree Input -->
    <div class="panel">
      <div class="panel-header">Tree Input</div>
      <div class="panel-body">
        <input type="text" class="file-path-input" id="tree-path"
               placeholder="Path to .tre / .nwk file"
               onkeydown="if(event.key==='Enter') loadTreeFile()">
        <button onclick="loadTreeFile()" class="btn-sm btn-full" style="margin-bottom:8px">Load Tree File</button>

        <div style="text-align:center; color:var(--text-light); font-size:12px; margin:4px 0;">or paste Newick below</div>
        <textarea id="newick-input" rows="4"
          style="width:100%; padding:8px; border:1px solid var(--border); border-radius:6px; font-family:monospace; font-size:11px; resize:vertical;"
          placeholder="((A:0.1,B:0.2):0.3,C:0.4);"></textarea>
        <button onclick="loadNewick()" class="btn-sm btn-full" style="margin-top:6px">Load Pasted Tree</button>
      </div>
    </div>

    <!-- Analysis -->
    <div class="panel">
      <div class="panel-header">Analysis</div>
      <div class="panel-body">
        <button onclick="runAnalysis()" class="btn-primary btn-full" style="margin-bottom:8px; font-size:14px; padding:12px;">
          Run Monophyly Check
        </button>
        <div class="btn-group">
          <button onclick="exportNewick()" class="btn-sm">Export Tree</button>
          <button onclick="exportReport()" class="btn-sm">Export Report</button>
        </div>
      </div>
    </div>

    <!-- Tree Info -->
    <div class="panel" id="tree-info-panel" style="display:none">
      <div class="panel-header">Tree Info</div>
      <div class="panel-body">
        <div id="tree-info" style="font-size:13px; line-height:1.8;"></div>
      </div>
    </div>
  </div>

  <!-- Main content -->
  <div class="content">
    <!-- Tab bar -->
    <div class="tab-bar">
      <div class="tab active" data-tab="tree" onclick="switchTab('tree')">Tree View</div>
      <div class="tab" data-tab="nonmono" onclick="switchTab('nonmono')">
        Non-Monophyletic <span class="badge badge-nonmono" id="badge-nonmono" style="display:none">0</span>
      </div>
      <div class="tab" data-tab="mono" onclick="switchTab('mono')">
        Monophyletic <span class="badge badge-mono" id="badge-mono" style="display:none">0</span>
      </div>
      <div class="tab" data-tab="unresolved" onclick="switchTab('unresolved')">
        Unresolved <span class="badge badge-unresolved" id="badge-unresolved" style="display:none">0</span>
      </div>
      <div class="tab" data-tab="summary" onclick="switchTab('summary')">Summary</div>
    </div>

    <!-- Tree View -->
    <div class="tab-content active" id="tab-tree">
      <div class="zoom-bar">
        <span>Zoom:</span>
        <input type="range" id="zoom-slider" min="0.3" max="4" step="0.1" value="1" oninput="redrawTree()">
        <span id="zoom-label">1.0x</span>
        <label><input type="checkbox" id="show-labels" checked onchange="redrawTree()"> Labels</label>
        <label><input type="checkbox" id="show-bl" checked onchange="redrawTree()"> Branch Lengths</label>
      </div>
      <div id="tree-container">
        <div class="empty-state" id="tree-empty">
          <div class="icon">&#127795;</div>
          <h3>No tree loaded</h3>
          <p>Load a Newick tree file using the sidebar, or paste a Newick string directly.</p>
        </div>
        <svg id="tree-svg" style="display:none"></svg>
      </div>
    </div>

    <!-- Non-mono tab -->
    <div class="tab-content" id="tab-nonmono">
      <div class="stats-row" id="stats-row" style="display:none">
        <div class="stat-card mono-card">
          <div class="stat-number" id="stat-mono">0</div>
          <div class="stat-label">Monophyletic</div>
        </div>
        <div class="stat-card nonmono-card">
          <div class="stat-number" id="stat-nonmono">0</div>
          <div class="stat-label">Non-Monophyletic</div>
        </div>
        <div class="stat-card unresolved-card">
          <div class="stat-number" id="stat-unresolved">0</div>
          <div class="stat-label">Unresolved Tips</div>
        </div>
      </div>
      <div class="results-table-container">
        <table id="nonmono-table">
          <thead>
            <tr><th>Clade Name</th><th>Rank</th><th># Tips</th><th>Intruders</th></tr>
          </thead>
          <tbody id="nonmono-tbody"></tbody>
        </table>
      </div>
    </div>

    <!-- Mono tab -->
    <div class="tab-content" id="tab-mono">
      <div class="results-table-container">
        <table id="mono-table">
          <thead>
            <tr><th>Clade Name</th><th>Rank</th><th># Tips</th></tr>
          </thead>
          <tbody id="mono-tbody"></tbody>
        </table>
      </div>
    </div>

    <!-- Unresolved tab -->
    <div class="tab-content" id="tab-unresolved">
      <div class="results-table-container">
        <div id="unresolved-content" class="summary-box" style="white-space:pre-wrap;"></div>
      </div>
    </div>

    <!-- Summary tab -->
    <div class="tab-content" id="tab-summary">
      <div class="summary-container">
        <div class="summary-box" id="summary-content">Run an analysis to see results here.</div>
      </div>
    </div>
  </div>
</div>

<!-- Tooltip -->
<div class="tooltip" id="tooltip"></div>

<!-- Loading overlay -->
<div class="loading-overlay" id="loading">
  <div class="spinner"></div>
  <div class="loading-text" id="loading-text">Loading...</div>
</div>

<script>
// ---- State ----
let treeData = null;
let resultsData = null;

// ---- API helpers ----
async function api(path, method='GET', body=null) {
  const opts = { method };
  if (body) {
    opts.headers = { 'Content-Type': 'application/json' };
    opts.body = JSON.stringify(body);
  }
  const res = await fetch(path, opts);
  return res.json();
}

function showLoading(msg) {
  document.getElementById('loading-text').textContent = msg;
  document.getElementById('loading').classList.add('active');
}
function hideLoading() {
  document.getElementById('loading').classList.remove('active');
}
function setStatus(msg) {
  document.getElementById('header-status').textContent = msg;
}

// ---- Taxonomy ----
async function loadCache() {
  showLoading('Loading taxonomy cache...');
  const r = await api('/api/load-cache', 'POST');
  hideLoading();
  if (r.success) {
    updateDbStatus(true, r.taxa);
  } else {
    updateDbStatus(false);
    alert(r.error);
  }
}

async function downloadNCBI() {
  if (!confirm('Download NCBI taxonomy (~60 MB)? This may take a few minutes.')) return;
  showLoading('Downloading NCBI taxonomy...');
  setStatus('Downloading NCBI taxonomy...');
  await api('/api/download-taxonomy', 'POST');
  // Poll for completion
  const poll = setInterval(async () => {
    const s = await api('/api/status');
    setStatus(s.status);
    if (s.db_loaded) {
      clearInterval(poll);
      hideLoading();
      updateDbStatus(true, s.db_taxa);
    }
  }, 2000);
}

async function loadFiles() {
  const names = document.getElementById('names-path').value.trim();
  const nodes = document.getElementById('nodes-path').value.trim();
  if (!names || !nodes) { alert('Enter paths to both names.dmp and nodes.dmp'); return; }
  showLoading('Loading taxonomy files...');
  const r = await api('/api/load-files', 'POST', { names_file: names, nodes_file: nodes });
  hideLoading();
  if (r.success) {
    updateDbStatus(true, r.taxa);
  } else {
    alert(r.error);
  }
}

function updateDbStatus(loaded, taxa) {
  const dot = document.getElementById('db-dot');
  const text = document.getElementById('db-status-text');
  dot.className = 'status-dot ' + (loaded ? 'loaded' : 'empty');
  text.textContent = loaded ? `${taxa.toLocaleString()} taxa loaded` : 'Not loaded';
  if (loaded) setStatus(`Taxonomy: ${taxa.toLocaleString()} taxa`);
}

// ---- Tree loading ----
async function loadTreeFile() {
  const path = document.getElementById('tree-path').value.trim();
  if (!path) { alert('Enter a tree file path'); return; }
  showLoading('Loading tree...');
  const r = await api('/api/load-tree-file', 'POST', { filepath: path });
  hideLoading();
  if (r.success) {
    treeData = r.tree;
    showTreeInfo(r.tips, r.internal);
    drawTree();
    setStatus(`Loaded: ${r.tips} tips, ${r.internal} internal nodes`);
  } else {
    alert(r.error);
  }
}

async function loadNewick() {
  const nwk = document.getElementById('newick-input').value.trim();
  if (!nwk) { alert('Paste a Newick string first'); return; }
  showLoading('Parsing tree...');
  const r = await api('/api/load-tree', 'POST', { newick: nwk, filename: 'pasted tree' });
  hideLoading();
  if (r.success) {
    treeData = r.tree;
    showTreeInfo(r.tips, r.internal);
    drawTree();
    setStatus(`Loaded: ${r.tips} tips`);
  } else {
    alert(r.error);
  }
}

function showTreeInfo(tips, internal) {
  document.getElementById('tree-info-panel').style.display = '';
  document.getElementById('tree-info').innerHTML =
    `<b>Tips:</b> ${tips}<br><b>Internal nodes:</b> ${internal}<br><b>Total:</b> ${tips + internal}`;
}

// ---- Analysis ----
async function runAnalysis() {
  showLoading('Running monophyly analysis...');
  const r = await api('/api/analyze', 'POST');
  hideLoading();
  if (!r.success) { alert(r.error); return; }

  treeData = r.tree;
  resultsData = r.results;
  drawTree();
  showResults(r.results, r.summary);

  const nm = r.results.total_non_mono;
  setStatus(`Done: ${r.results.total_mono} monophyletic, ${nm} non-monophyletic`);

  // Switch to non-mono tab (PI preference: mismatches first)
  if (nm > 0) switchTab('nonmono');
}

function showResults(results, summary) {
  // Badges
  const bm = document.getElementById('badge-mono');
  const bnm = document.getElementById('badge-nonmono');
  const bu = document.getElementById('badge-unresolved');
  bm.textContent = results.total_mono; bm.style.display = '';
  bnm.textContent = results.total_non_mono; bnm.style.display = '';
  bu.textContent = results.total_unresolved; bu.style.display = '';

  // Stats
  document.getElementById('stats-row').style.display = '';
  document.getElementById('stat-mono').textContent = results.total_mono;
  document.getElementById('stat-nonmono').textContent = results.total_non_mono;
  document.getElementById('stat-unresolved').textContent = results.total_unresolved;

  // Non-mono table
  const nmtb = document.getElementById('nonmono-tbody');
  nmtb.innerHTML = '';
  results.non_monophyletic.forEach(r => {
    const intr = r.intruders.join(', ') + (r.total_intruders > 20 ? ` ... (+${r.total_intruders - 20} more)` : '');
    nmtb.innerHTML += `<tr><td><b>${esc(r.name)}</b></td><td>${esc(r.rank)}</td><td>${r.tips}</td><td class="intruder-list" title="${esc(intr)}">${esc(intr)}</td></tr>`;
  });

  // Mono table
  const mtb = document.getElementById('mono-tbody');
  mtb.innerHTML = '';
  results.monophyletic.forEach(r => {
    mtb.innerHTML += `<tr><td><b>${esc(r.name)}</b></td><td>${esc(r.rank)}</td><td>${r.tips}</td></tr>`;
  });

  // Unresolved
  const uc = document.getElementById('unresolved-content');
  if (results.unresolved_tips.length) {
    uc.textContent = `${results.total_unresolved} tip(s) could not be matched to NCBI taxonomy:\n\n` +
      results.unresolved_tips.map(t => '  ' + t).join('\n');
  } else {
    uc.textContent = 'All tips resolved successfully.';
  }

  // Summary
  document.getElementById('summary-content').textContent = summary;
}

function esc(s) { const d = document.createElement('div'); d.textContent = s; return d.innerHTML; }

// ---- Export ----
async function exportNewick() {
  const r = await api('/api/export-newick', 'POST');
  if (!r.success) { alert(r.error); return; }
  downloadFile('labeled_tree.nwk', r.newick);
}

async function exportReport() {
  const r = await api('/api/export-report', 'POST');
  if (!r.success) { alert(r.error); return; }
  downloadFile('analysis_report.txt', r.report);
}

function downloadFile(name, content) {
  const a = document.createElement('a');
  a.href = 'data:text/plain;charset=utf-8,' + encodeURIComponent(content);
  a.download = name;
  a.click();
}

// ---- Tabs ----
function switchTab(name) {
  document.querySelectorAll('.tab').forEach(t => t.classList.toggle('active', t.dataset.tab === name));
  document.querySelectorAll('.tab-content').forEach(t => t.classList.toggle('active', t.id === 'tab-' + name));
}

// ---- Tree Drawing (SVG) ----
function drawTree() {
  if (!treeData) return;

  const svg = document.getElementById('tree-svg');
  const container = document.getElementById('tree-container');
  const empty = document.getElementById('tree-empty');
  empty.style.display = 'none';
  svg.style.display = '';

  const zoom = parseFloat(document.getElementById('zoom-slider').value);
  document.getElementById('zoom-label').textContent = zoom.toFixed(1) + 'x';
  const showLabels = document.getElementById('show-labels').checked;

  const tipSpacing = 20 * zoom;
  const marginLeft = 30;
  const marginTop = 20;

  // Count tips
  function countTips(node) {
    if (node.is_tip) return 1;
    return (node.children || []).reduce((s, c) => s + countTips(c), 0);
  }
  const nTips = countTips(treeData);

  // Get max depth
  function getDepth(node) {
    if (node.is_tip) return 0;
    return 1 + Math.max(...(node.children || []).map(getDepth));
  }
  const maxDepth = Math.max(1, getDepth(treeData));
  const xScale = Math.max(300, 160 * zoom * Math.min(maxDepth, 30) / 10);

  // Assign y positions
  let tipIdx = 0;
  function assignY(node) {
    if (node.is_tip) {
      node._y = marginTop + tipIdx * tipSpacing;
      tipIdx++;
      return;
    }
    (node.children || []).forEach(assignY);
    const ys = node.children.map(c => c._y);
    node._y = (Math.min(...ys) + Math.max(...ys)) / 2;
  }

  // Assign x positions (depth-based)
  function assignX(node, depth) {
    node._x = marginLeft + (depth / maxDepth) * xScale;
    node._depth = depth;
    (node.children || []).forEach(c => assignX(c, depth + 1));
  }

  assignY(treeData);
  assignX(treeData, 0);

  // Build SVG
  const totalW = marginLeft + xScale + 250 * zoom;
  const totalH = marginTop + nTips * tipSpacing + 30;
  svg.setAttribute('width', totalW);
  svg.setAttribute('height', totalH);
  svg.setAttribute('viewBox', `0 0 ${totalW} ${totalH}`);

  let html = '';

  function drawNode(node) {
    (node.children || []).forEach(child => {
      let cls = 'tree-branch';
      if (!child.is_tip && child.monophyly === 'monophyletic') cls += ' mono';
      else if (!child.is_tip && child.monophyly === 'not_monophyletic') cls += ' non-mono';

      // Rectangular phylogram: vertical then horizontal
      html += `<line class="${cls}" x1="${node._x}" y1="${node._y}" x2="${node._x}" y2="${child._y}"/>`;
      html += `<line class="${cls}" x1="${node._x}" y1="${child._y}" x2="${child._x}" y2="${child._y}"/>`;

      drawNode(child);
    });

    if (showLabels) {
      if (node.is_tip) {
        html += `<text class="tree-tip" x="${node._x + 5}" y="${node._y + 4}" font-size="${Math.max(9, 12 * zoom)}">${esc(node.label)}</text>`;
      } else if (node.label || node.tax_name) {
        const name = node.tax_name || node.label;
        let cls = 'tree-label';
        if (node.monophyly === 'monophyletic') cls += ' mono';
        else if (node.monophyly === 'not_monophyletic') cls += ' non-mono';
        else cls += ' unresolved';
        html += `<text class="${cls}" x="${node._x - 4}" y="${node._y - 6 * zoom}" text-anchor="end" font-size="${Math.max(8, 10 * zoom)}">${esc(name)}</text>`;
      }
    }
  }

  drawNode(treeData);
  svg.innerHTML = html;
}

function redrawTree() { drawTree(); }

// ---- Init: check if cache exists ----
(async function init() {
  const s = await api('/api/status');
  if (s.db_loaded) {
    updateDbStatus(true, s.db_taxa);
  } else {
    const info = await api('/api/cache-info');
    if (info.exists) loadCache();
  }
})();
</script>
</body>
</html>
"""


def main():
    """Start the PhyLabeler web server and open the browser."""
    # Use SO_REUSEADDR to avoid "address already in use"
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("", PORT), PhyLabelerHandler) as httpd:
        url = f"http://localhost:{PORT}"
        print(f"")
        print(f"  PhyLabeler is running at: {url}")
        print(f"  Open this URL in your browser.")
        print(f"  Press Ctrl+C to stop.")
        print(f"")

        # Open browser automatically
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()

        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nShutting down.")
            httpd.shutdown()


if __name__ == "__main__":
    main()
