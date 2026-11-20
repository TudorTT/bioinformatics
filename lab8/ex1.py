"""
Needleman-Wunsch global alignment of two DNA sequences - Tkinter app.
Plain Python (no Biopython).

Run:  python needleman_wunsch_app.py
"""
import tkinter as tk
from tkinter import ttk, messagebox

# ---------------------------------------------------------------------------
# 1. The algorithm
# ---------------------------------------------------------------------------

def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    """Return (aligned_s1, aligned_s2, score, score_matrix)."""
    n, m = len(s1), len(s2)

    # Step 1: build the (n+1) x (m+1) score matrix
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # Step 2: initialise first row and column with cumulative gap penalties
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap

    # Step 3: fill the matrix - each cell takes the best of 3 moves
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if s1[i - 1] == s2[j - 1] else mismatch
            diag = F[i - 1][j - 1] + s      # match / mismatch
            up = F[i - 1][j] + gap          # gap in s2
            left = F[i][j - 1] + gap        # gap in s1
            F[i][j] = max(diag, up, left)

    # Step 4: traceback from bottom-right corner to top-left
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if s1[i - 1] == s2[j - 1] else mismatch
            if F[i][j] == F[i - 1][j - 1] + s:
                a1.append(s1[i - 1]); a2.append(s2[j - 1])
                i -= 1; j -= 1
                continue
        if i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1.append(s1[i - 1]); a2.append("-")
            i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1])
            j -= 1

    return "".join(reversed(a1)), "".join(reversed(a2)), F[n][m], F


def match_line(a1, a2):
    """'|' for a match, '.' for a mismatch, ' ' for a gap."""
    return "".join("|" if x == y else (" " if "-" in (x, y) else ".")
                   for x, y in zip(a1, a2))


def stats(a1, a2):
    matches = sum(x == y for x, y in zip(a1, a2))
    gaps = sum("-" in (x, y) for x, y in zip(a1, a2))
    mismatches = len(a1) - matches - gaps
    identity = 100 * matches / len(a1) if a1 else 0
    return matches, mismatches, gaps, identity


# ---------------------------------------------------------------------------
# 2. The GUI
# ---------------------------------------------------------------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Needleman-Wunsch DNA Aligner")
        self.geometry("900x650")

        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)

        # Sequence inputs
        ttk.Label(frm, text="Sequence 1:").grid(row=0, column=0, sticky="w")
        self.e1 = ttk.Entry(frm, width=70, font=("Courier", 11))
        self.e1.grid(row=0, column=1, columnspan=5, sticky="we", pady=2)
        self.e1.insert(0, "ACCGTGAAGCCAATAC")

        ttk.Label(frm, text="Sequence 2:").grid(row=1, column=0, sticky="w")
        self.e2 = ttk.Entry(frm, width=70, font=("Courier", 11))
        self.e2.grid(row=1, column=1, columnspan=5, sticky="we", pady=2)
        self.e2.insert(0, "AGCGTGCAGCCAATAC")

        # Scoring parameters
        self.vm = tk.IntVar(value=1)
        self.vx = tk.IntVar(value=-1)
        self.vg = tk.IntVar(value=-2)
        for col, (lbl, var) in enumerate(
                [("Match", self.vm), ("Mismatch", self.vx), ("Gap", self.vg)]):
            ttk.Label(frm, text=lbl + ":").grid(row=2, column=col * 2, sticky="e", pady=6)
            ttk.Spinbox(frm, from_=-10, to=10, textvariable=var, width=5)\
                .grid(row=2, column=col * 2 + 1, sticky="w")

        ttk.Button(frm, text="Align", command=self.run)\
            .grid(row=3, column=0, columnspan=6, pady=6)

        # Output
        self.out = tk.Text(frm, font=("Courier", 11), wrap="none")
        self.out.grid(row=4, column=0, columnspan=6, sticky="nsew")
        sy = ttk.Scrollbar(frm, orient="vertical", command=self.out.yview)
        sx = ttk.Scrollbar(frm, orient="horizontal", command=self.out.xview)
        self.out.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.grid(row=4, column=6, sticky="ns")
        sx.grid(row=5, column=0, columnspan=6, sticky="we")

        frm.columnconfigure(1, weight=1)
        frm.rowconfigure(4, weight=1)

    def run(self):
        s1 = self.e1.get().strip().upper()
        s2 = self.e2.get().strip().upper()
        if not s1 or not s2:
            messagebox.showerror("Error", "Both sequences are required.")
            return
        if set(s1 + s2) - set("ACGT"):
            messagebox.showerror("Error", "Sequences may only contain A, C, G, T.")
            return
        try:
            m, x, g = self.vm.get(), self.vx.get(), self.vg.get()
        except tk.TclError:
            messagebox.showerror("Error", "Scores must be integers.")
            return

        a1, a2, score, F = needleman_wunsch(s1, s2, m, x, g)
        mt, mm, gp, ident = stats(a1, a2)

        lines = [
            "=== ALIGNMENT ===",
            f"Seq1: {a1}",
            f"      {match_line(a1, a2)}",
            f"Seq2: {a2}",
            "",
            f"Score: {score}",
            f"Length: {len(a1)} | Matches: {mt} | Mismatches: {mm} | Gaps: {gp}",
            f"Identity: {ident:.2f}%",
            "",
            "=== SCORE MATRIX ===",
            "      " + "".join(f"{c:>4}" for c in "-" + s2),
        ]
        for i, row in enumerate(F):
            label = "-" if i == 0 else s1[i - 1]
            lines.append(f"{label:>4}  " + "".join(f"{v:>4}" for v in row))

        self.out.delete("1.0", "end")
        self.out.insert("1.0", "\n".join(lines))


if __name__ == "__main__":
    App().mainloop()