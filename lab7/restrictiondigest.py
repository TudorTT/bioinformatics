"""
Restriction enzyme digestion simulator + virtual agarose gel

Run the GUI:      python restriction_digest.py
Run in console:   python restriction_digest.py --cli genome.fasta [--circular]

Enzyme notation: the recognition site with '^' marking where the TOP strand
(5'->3') is cut, e.g. EcoRI = G^AATTC  ->  5'...G | AATTC...3'
"""

import re
import sys
import threading
import urllib.request
import urllib.parse
from itertools import combinations

# --------------------------------------------------------------------------
# 1. ENZYME DATA
# --------------------------------------------------------------------------
# name -> site written with '^' at the top-strand cleavage point
DEFAULT_ENZYMES = {
    "EcoRI":   "G^AATTC",   # sticky ends, 5' overhang AATT
    "BamHI":   "G^GATCC",   # sticky ends, 5' overhang GATC
    "HindIII": "A^AGCTT",   # sticky ends, 5' overhang AGCT
    "TaqI":    "T^CGA",     # sticky ends, 5' overhang CG
    "HaeIII":  "GG^CC",     # blunt ends
}

# IUPAC codes so custom enzymes like "GANTC" or "CCWGG" also work
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def parse_site(site_with_cut):
    """'G^AATTC' -> ('GAATTC', 1)   (1 = cut after the 1st base of the site)"""
    s = site_with_cut.strip().upper().replace(" ", "")
    if s.count("^") != 1:
        raise ValueError(f"Site '{site_with_cut}' must contain exactly one '^' (e.g. G^AATTC)")
    offset = s.index("^")
    site = s.replace("^", "")
    bad = [c for c in site if c not in IUPAC]
    if not site or bad:
        raise ValueError(f"Site '{site_with_cut}' has invalid letters: {''.join(bad)}")
    return site, offset


def site_to_regex(site):
    """Builds a regex with a lookahead so OVERLAPPING sites are found too."""
    return re.compile("(?=" + "".join(IUPAC[c] for c in site) + ")")


def show_cut(site_with_cut):
    """Pretty two-strand drawing of the cut, like in the lab sheet."""
    site, off = parse_site(site_with_cut)
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    bottom = "".join(comp.get(c, "N") for c in site)
    bottom_off = len(site) - off      # palindromic sites: bottom cut is mirrored
    top = f"5'-{site[:off]} {site[off:]}-3'"
    bot = f"3'-{bottom[:bottom_off]} {bottom[bottom_off:]}-5'"
    return top, bot


# --------------------------------------------------------------------------
# 2. SEQUENCE INPUT
# --------------------------------------------------------------------------
def clean_sequence(text):
    """Accepts raw text or FASTA (header lines start with '>'), returns ACGT... string."""
    lines = [ln for ln in text.splitlines() if not ln.startswith(">")]
    return re.sub(r"[^A-Za-z]", "", "".join(lines)).upper()


def read_fasta_header(text):
    for ln in text.splitlines():
        if ln.startswith(">"):
            return ln[1:].strip()
    return ""


def fetch_from_ncbi(accession):
    """Downloads a FASTA record from NCBI (needs internet). E.g. NC_001416.1 = phage lambda."""
    params = urllib.parse.urlencode({
        "db": "nuccore", "id": accession.strip(),
        "rettype": "fasta", "retmode": "text",
    })
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + params
    with urllib.request.urlopen(url, timeout=60) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    if not text.startswith(">"):
        raise ValueError(f"NCBI did not return a FASTA record for '{accession}':\n{text[:200]}")
    return text


# --------------------------------------------------------------------------
# 3. DIGESTION LOGIC
# --------------------------------------------------------------------------
def find_cut_positions(seq, site_with_cut, circular=False):
    """
    Returns the sorted list of cut positions.
    A cut position p means: the top strand is cut between base p and base p+1
    (1-based), i.e. the first fragment ends at base p.
    """
    site, offset = parse_site(site_with_cut)
    pattern = site_to_regex(site)
    n = len(seq)
    # circular DNA: also look for sites that cross the end -> start junction
    search_in = seq + seq[:len(site) - 1] if circular else seq
    cuts = set()
    for m in pattern.finditer(search_in):
        start = m.start()
        if start >= n:                      # only matches that START inside the sequence
            continue
        p = start + offset                  # 0-based index of first base AFTER the cut
        if circular:
            p %= n
        elif p == 0 or p == n:              # cut at the very end of a linear molecule = no cut
            continue
        cuts.add(p)
    return sorted(cuts)


def make_fragments(cuts, n, circular=False):
    """
    cuts: sorted 0-based 'first base after cut' positions.
    Returns list of (start, end, length) with 1-based inclusive coordinates.
    """
    if not cuts:
        return [(1, n, n)]                  # uncut: one piece (circular stays a circle)
    frags = []
    if circular:
        # pieces between consecutive cuts + the piece that wraps around the origin
        for a, b in zip(cuts, cuts[1:]):
            frags.append((a + 1, b, b - a))
        last, first = cuts[-1], cuts[0]
        frags.append((last + 1, first if first > 0 else n, n - last + first))
    else:
        bounds = [0] + cuts + [n]
        for a, b in zip(bounds, bounds[1:]):
            frags.append((a + 1, b, b - a))
    return frags


def digest(seq, enzyme_names, enzymes, circular=False):
    """Digest with one OR several enzymes at the same time (a combination)."""
    cut_owner = {}                          # position -> list of enzymes that cut there
    per_enzyme = {}
    for name in enzyme_names:
        pos = find_cut_positions(seq, enzymes[name], circular)
        per_enzyme[name] = len(pos)
        for p in pos:
            cut_owner.setdefault(p, []).append(name)
    cuts = sorted(cut_owner)
    return {
        "label": " + ".join(enzyme_names),
        "enzymes": list(enzyme_names),
        "per_enzyme": per_enzyme,
        "cuts": cuts,
        "cut_owner": cut_owner,
        "fragments": make_fragments(cuts, len(seq), circular),
        "circular": circular,
    }


def all_digests(seq, enzyme_names, enzymes, circular=False, combos=True, max_combo=None):
    """Single digests first, then every combination of 2, 3, ... enzymes."""
    results = []
    top = len(enzyme_names) if combos else 1
    if max_combo:
        top = min(top, max_combo)
    for k in range(1, top + 1):
        for combo in combinations(enzyme_names, k):
            results.append(digest(seq, combo, enzymes, circular))
    return results


def format_result(res, seq_len, max_rows=None):
    """Text report for one digestion."""
    out = []
    out.append("=" * 70)
    out.append(f"  {res['label']}")
    out.append("=" * 70)
    if len(res["enzymes"]) > 1:
        detail = ", ".join(f"{e}: {c}" for e, c in res["per_enzyme"].items())
        out.append(f"Cleavages per enzyme : {detail}")
    shared = sum(1 for v in res["cut_owner"].values() if len(v) > 1)
    out.append(f"Number of cleavages  : {len(res['cuts'])}"
               + (f"  ({shared} position(s) cut by more than one enzyme)" if shared else ""))
    out.append(f"Number of fragments  : {len(res['fragments'])}"
               + ("  (circular DNA)" if res["circular"] else "  (linear DNA)"))

    if res["cuts"]:
        out.append("\nCleavage positions (cut after base p, top strand):")
        rows = res["cuts"] if max_rows is None else res["cuts"][:max_rows]
        for i, p in enumerate(rows, 1):
            shown = p if p > 0 else seq_len  # circular cut at origin -> after last base
            out.append(f"  {i:>5}. after base {shown:>10,}   [{', '.join(res['cut_owner'][p])}]")
        if max_rows is not None and len(res["cuts"]) > max_rows:
            out.append(f"  ... {len(res['cuts']) - max_rows:,} more (use Export for the full list)")

    out.append("\nFragments (in sequence order):")
    out.append(f"  {'#':>5}  {'start':>10}  {'end':>10}  {'length (bp)':>12}")
    frags = res["fragments"] if max_rows is None else res["fragments"][:max_rows]
    for i, (s, e, L) in enumerate(frags, 1):
        out.append(f"  {i:>5}  {s:>10,}  {e:>10,}  {L:>12,}")
    if max_rows is not None and len(res["fragments"]) > max_rows:
        out.append(f"  ... {len(res['fragments']) - max_rows:,} more (use Export for the full list)")

    lengths = sorted((f[2] for f in res["fragments"]), reverse=True)
    preview = ", ".join(f"{L:,}" for L in lengths[:15])
    out.append(f"\nSizes, largest -> smallest (gel order): {preview}"
               + (" ..." if len(lengths) > 15 else ""))
    out.append("")
    return "\n".join(out)


def summary_table(results):
    out = ["SUMMARY", f"{'Digestion':<45}{'cuts':>8}{'fragments':>11}{'largest':>12}{'smallest':>10}"]
    for r in results:
        L = [f[2] for f in r["fragments"]]
        out.append(f"{r['label']:<45}{len(r['cuts']):>8}{len(L):>11}{max(L):>12,}{min(L):>10,}")
    return "\n".join(out) + "\n\n"


# --------------------------------------------------------------------------
# 4. GEL MODEL  (migration distance ~ log10(fragment length))
# --------------------------------------------------------------------------
LADDER = [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000,
          3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000]


def migration(length, min_bp, max_bp):
    """0.0 = at the well (big fragments), 1.0 = far end (small fragments)."""
    import math
    L = min(max(length, min_bp), max_bp)
    return (math.log10(max_bp) - math.log10(L)) / (math.log10(max_bp) - math.log10(min_bp))


# --------------------------------------------------------------------------
# 5. CONSOLE MODE
# --------------------------------------------------------------------------
def run_cli(argv):
    if len(argv) < 1:
        print("usage: python restriction_digest.py --cli file.fasta [--circular] [--no-combos]")
        return
    with open(argv[0], encoding="utf-8") as f:
        text = f.read()
    seq = clean_sequence(text)
    circular = "--circular" in argv
    combos = "--no-combos" not in argv
    print(f"Sequence: {read_fasta_header(text) or argv[0]}  |  {len(seq):,} bp  |  "
          f"{'circular' if circular else 'linear'}\n")
    print("Enzymes:")
    for name, s in DEFAULT_ENZYMES.items():
        t, b = show_cut(s)
        print(f"  {name:<8} {s:<9} {t}\n  {'':<8} {'':<9} {b}")
    print()
    results = all_digests(seq, list(DEFAULT_ENZYMES), DEFAULT_ENZYMES, circular, combos)
    print(summary_table(results))
    for r in results:
        print(format_result(r, len(seq), max_rows=50))


# --------------------------------------------------------------------------
# 6. TKINTER GUI
# --------------------------------------------------------------------------
def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("Restriction Digest Simulator")
            self.geometry("1250x780")
            self.enzymes = dict(DEFAULT_ENZYMES)
            self.enzyme_vars = {}
            self.sequence = ""          # source of truth when loaded from file / NCBI
            self.seq_name = ""
            self.results = []
            self._build()

        # ---------------- layout ----------------
        def _build(self):
            left = ttk.Frame(self, padding=6)
            left.pack(side="left", fill="y")
            right = ttk.Frame(self, padding=6)
            right.pack(side="left", fill="both", expand=True)

            # --- sequence input
            box = ttk.LabelFrame(left, text="1. DNA sequence", padding=6)
            box.pack(fill="x")
            row = ttk.Frame(box)
            row.pack(fill="x")
            ttk.Button(row, text="Load FASTA...", command=self.load_fasta).pack(side="left")
            ttk.Label(row, text="  NCBI accession:").pack(side="left")
            self.acc = ttk.Entry(row, width=13)
            self.acc.insert(0, "NC_001416.1")
            self.acc.pack(side="left")
            ttk.Button(row, text="Fetch", command=self.fetch).pack(side="left", padx=3)
            ttk.Label(box, text="...or paste a sequence / FASTA here:").pack(anchor="w", pady=(6, 0))
            self.seq_text = scrolledtext.ScrolledText(box, width=52, height=9, font=("Courier", 9))
            self.seq_text.pack(fill="x")
            self.seq_text.edit_modified(False)
            self.seq_info = ttk.Label(box, text="No sequence loaded", foreground="gray")
            self.seq_info.pack(anchor="w")
            self.circular = tk.BooleanVar(value=False)
            ttk.Checkbutton(box, text="Circular DNA (plasmid / circular genome)",
                            variable=self.circular).pack(anchor="w")

            # --- enzymes
            ebox = ttk.LabelFrame(left, text="2. Restriction enzymes", padding=6)
            ebox.pack(fill="x", pady=6)
            self.enz_frame = ttk.Frame(ebox)
            self.enz_frame.pack(fill="x")
            self._refresh_enzyme_list()
            add = ttk.Frame(ebox)
            add.pack(fill="x", pady=(6, 0))
            ttk.Label(add, text="Name").grid(row=0, column=0)
            ttk.Label(add, text="Site with ^ (e.g. GT^MKAC)").grid(row=0, column=1)
            self.new_name = ttk.Entry(add, width=10)
            self.new_name.grid(row=1, column=0, padx=2)
            self.new_site = ttk.Entry(add, width=18)
            self.new_site.grid(row=1, column=1, padx=2)
            ttk.Button(add, text="Add", command=self.add_enzyme).grid(row=1, column=2)

            # --- run
            rbox = ttk.LabelFrame(left, text="3. Digest", padding=6)
            rbox.pack(fill="x")
            self.combos = tk.BooleanVar(value=True)
            ttk.Checkbutton(rbox, text="Also every combination of the selected enzymes",
                            variable=self.combos).pack(anchor="w")
            row = ttk.Frame(rbox)
            row.pack(fill="x", pady=4)
            ttk.Button(row, text="DIGEST", command=self.run_digest).pack(side="left")
            ttk.Button(row, text="Export full results...", command=self.export).pack(side="left", padx=6)

            # --- right side: results + gel
            nb = ttk.Notebook(right)
            nb.pack(fill="both", expand=True)
            self.nb = nb
            res_tab = ttk.Frame(nb)
            nb.add(res_tab, text="Results")
            self.out = scrolledtext.ScrolledText(res_tab, font=("Courier", 9), wrap="none")
            self.out.pack(fill="both", expand=True)

            gel_tab = ttk.Frame(nb)
            nb.add(gel_tab, text="Gel")
            side = ttk.Frame(gel_tab, padding=4)
            side.pack(side="left", fill="y")
            ttk.Label(side, text="Lanes to load\n(Ctrl/Shift-click):").pack(anchor="w")
            self.lane_list = tk.Listbox(side, selectmode="extended", width=34, height=30,
                                        exportselection=False)
            self.lane_list.pack(fill="y", expand=True)
            ttk.Button(side, text="Run gel", command=self.draw_gel).pack(fill="x", pady=4)
            self.gel = tk.Canvas(gel_tab, bg="#14141c", highlightthickness=0)
            self.gel.pack(side="left", fill="both", expand=True)
            self.gel.bind("<Configure>", lambda e: self.draw_gel())

        def _refresh_enzyme_list(self):
            for w in self.enz_frame.winfo_children():
                w.destroy()
            for name, site in self.enzymes.items():
                var = self.enzyme_vars.setdefault(name, tk.BooleanVar(value=True))
                top, bot = show_cut(site)
                ttk.Checkbutton(self.enz_frame, text=f"{name:<8} {site:<10} {top}",
                                variable=var).pack(anchor="w")

        # ---------------- input ----------------
        def _set_sequence(self, text, name):
            self.sequence = clean_sequence(text)
            self.seq_name = read_fasta_header(text) or name
            self.seq_text.delete("1.0", "end")
            if len(self.sequence) <= 200_000:
                self.seq_text.insert("1.0", self.sequence)
            else:   # don't put megabases into a Text widget - show a preview instead
                self.seq_text.insert("1.0", f">{self.seq_name}\n(large sequence - preview only)\n"
                                     + self.sequence[:3000] + "...")
            self.seq_text.edit_modified(False)
            gc = 100 * (self.sequence.count("G") + self.sequence.count("C")) / max(1, len(self.sequence))
            self.seq_info.config(text=f"{self.seq_name[:60]}\n{len(self.sequence):,} bp, GC = {gc:.1f}%",
                                 foreground="black")

        def load_fasta(self):
            path = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta *.fa *.fna *.txt"), ("All", "*.*")])
            if path:
                with open(path, encoding="utf-8", errors="replace") as f:
                    self._set_sequence(f.read(), path.split("/")[-1])

        def fetch(self):
            acc = self.acc.get().strip()
            if not acc:
                return
            self.seq_info.config(text=f"Downloading {acc} from NCBI...", foreground="blue")
            self.config(cursor="watch")

            def worker():
                try:
                    text = fetch_from_ncbi(acc)
                    self.after(0, lambda: (self._set_sequence(text, acc), self.config(cursor="")))
                except Exception as exc:
                    msg = str(exc)
                    self.after(0, lambda: (self.config(cursor=""),
                                           self.seq_info.config(text="Download failed", foreground="red"),
                                           messagebox.showerror("NCBI", msg)))
            threading.Thread(target=worker, daemon=True).start()

        def _current_sequence(self):
            # if the user typed/pasted in the box, that wins over the loaded file
            if self.seq_text.edit_modified() or not self.sequence:
                text = self.seq_text.get("1.0", "end")
                self.sequence = clean_sequence(text)
                self.seq_name = read_fasta_header(text) or "pasted sequence"
                self.seq_text.edit_modified(False)
            return self.sequence

        def add_enzyme(self):
            name, site = self.new_name.get().strip(), self.new_site.get().strip()
            try:
                if not name:
                    raise ValueError("Give the enzyme a name.")
                parse_site(site)
            except ValueError as exc:
                messagebox.showerror("Enzyme", str(exc))
                return
            self.enzymes[name] = site.upper()
            self._refresh_enzyme_list()
            self.new_name.delete(0, "end")
            self.new_site.delete(0, "end")

        # ---------------- digestion ----------------
        def run_digest(self):
            seq = self._current_sequence()
            if not seq:
                messagebox.showwarning("Digest", "Load, fetch or paste a DNA sequence first.")
                return
            chosen = [n for n in self.enzymes if self.enzyme_vars[n].get()]
            if not chosen:
                messagebox.showwarning("Digest", "Select at least one enzyme.")
                return
            self.config(cursor="watch")
            self.update_idletasks()
            self.results = all_digests(seq, chosen, self.enzymes, self.circular.get(), self.combos.get())
            self.config(cursor="")

            self.out.delete("1.0", "end")
            self.out.insert("end", f"Sequence: {self.seq_name}  |  {len(seq):,} bp  |  "
                                   f"{'circular' if self.circular.get() else 'linear'}\n\n")
            self.out.insert("end", summary_table(self.results))
            for r in self.results:
                self.out.insert("end", format_result(r, len(seq), max_rows=200))

            # gel lanes: preselect single digests + the all-enzymes digest
            self.lane_list.delete(0, "end")
            for i, r in enumerate(self.results):
                self.lane_list.insert("end", r["label"])
                if len(r["enzymes"]) == 1 or len(r["enzymes"]) == len(chosen):
                    self.lane_list.selection_set(i)
            self.draw_gel()

        def export(self):
            if not self.results:
                messagebox.showinfo("Export", "Run a digestion first.")
                return
            path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text report", "*.txt"), ("CSV fragments", "*.csv")])
            if not path:
                return
            with open(path, "w", encoding="utf-8") as f:
                if path.lower().endswith(".csv"):
                    f.write("digestion,fragment_no,start,end,length_bp\n")
                    for r in self.results:
                        for i, (s, e, L) in enumerate(r["fragments"], 1):
                            f.write(f"\"{r['label']}\",{i},{s},{e},{L}\n")
                else:
                    f.write(f"Sequence: {self.seq_name} | {len(self.sequence):,} bp\n\n")
                    f.write(summary_table(self.results))
                    for r in self.results:
                        f.write(format_result(r, len(self.sequence)))
            messagebox.showinfo("Export", f"Saved to {path}")

        # ---------------- gel ----------------
        def draw_gel(self):
            c = self.gel
            c.delete("all")
            W, H = c.winfo_width(), c.winfo_height()
            if W < 50 or not self.results:
                c.create_text(W / 2, H / 2, text="Run a digestion to load the gel", fill="#888")
                return
            lanes = [self.results[i] for i in self.lane_list.curselection()]
            all_lengths = [f[2] for r in lanes for f in r["fragments"]] + LADDER
            min_bp = max(50, min(all_lengths) * 0.8)
            max_bp = max(all_lengths) * 1.2

            n_lanes = len(lanes) + 1
            top, bottom, left = 70, H - 30, 60
            lane_w = min(90, (W - left - 20) / n_lanes)
            band_w = lane_w * 0.7

            def y_of(L):
                return top + migration(L, min_bp, max_bp) * (bottom - top)

            c.create_text(W / 2, 14, text=f"Virtual 1% agarose gel - {self.seq_name[:70]}",
                          fill="#ddd", font=("Arial", 10, "bold"))
            c.create_text(W / 2, 30, text="(-) wells                 migration -> (+)",
                          fill="#777", font=("Arial", 8))

            all_lanes = [("Ladder", [(0, 0, L) for L in LADDER], True)] + \
                        [(r["label"], r["fragments"], False) for r in lanes]
            for i, (label, frags, is_ladder) in enumerate(all_lanes):
                x = left + i * lane_w + lane_w / 2
                c.create_rectangle(x - band_w / 2, top - 12, x + band_w / 2, top - 4,
                                   outline="#555")                       # the well
                short = label if len(label) < 14 else label.replace(" + ", "+")[:22]
                c.create_text(x, top - 26, text=short, fill="#ccc", font=("Arial", 7),
                              width=lane_w - 2)
                # merge co-migrating fragments (same pixel row) -> brighter band
                rows = {}
                for (_, _, L) in frags:
                    y = round(y_of(L))
                    rows[y] = rows.get(y, 0) + L          # DNA mass ~ length (EtBr staining)
                peak = max(rows.values())
                for y, mass in rows.items():
                    k = 0.35 + 0.65 * (mass / peak) ** 0.5
                    col = "#%02x%02x%02x" % (int(255 * k), int(170 * k), int(60 * k)) if not is_ladder \
                        else "#%02x%02x%02x" % (int(150 * k), int(200 * k), int(255 * k))
                    c.create_line(x - band_w / 2, y, x + band_w / 2, y, fill=col, width=3)

            for L in LADDER:                              # ladder size labels
                c.create_text(left - 6, y_of(L), text=f"{L:,}", fill="#9ac", anchor="e",
                              font=("Arial", 7))

    App().mainloop()


if __name__ == "__main__":
    if "--cli" in sys.argv:
        run_cli(sys.argv[sys.argv.index("--cli") + 1:])
    else:
        run_gui()