"""
Transposable element (TE) simulator + detector
Plain Python, standard library only (random, itertools, array, tkinter) - no Biopython.

PART 1  - build an artificial 200-400 bp DNA sequence with 3-4 TEs (some of them
          nested / overlapping) and detect them using their KNOWN inverted repeats.
PART 2  - search a whole genome (FASTA, up to ~6 Mb) for NEW candidate TEs whose
          inverted repeats are UNKNOWN, with the s1 / s2 sliding-window idea.

A transposon (DNA transposon / bacterial IS element) looks like this:

   TSD        left IR                                   right IR        TSD
  [ATGC] 5'-[GTAACGGCTA]--------- transposase gene ----[TAGCCGTTAC]-3' [ATGC]
                  s1                                         s2 = reverse complement of s1

  - the two ends are INVERTED REPEATS: right IR = reverse complement of left IR
  - insertion duplicates a few bases of the target -> Target Site Duplication (TSD)

Run:  python transposon_finder.py                     (GUI)
      python transposon_finder.py --cli-part1 [scenario] [seed]
      python transposon_finder.py --cli-genome genome.fasta [k] [min_gap] [max_gap] [min_ir]
"""

import random
import sys
import time
from array import array

COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def revcomp(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def clean_sequence(text):
    lines = [ln for ln in text.splitlines() if not ln.startswith(">")]
    return "".join("".join(lines).split()).upper()


def first_record(text):
    """Multi-record FASTA (chromosome + plasmids): keep only the first record, otherwise
    the records would be glued together and fake IR pairs could span the junction."""
    if text.count(">") <= 1:
        return text, 1
    recs = text.split(">")[1:]
    return ">" + recs[0], len(recs)


def fasta_header(text):
    for ln in text.splitlines():
        if ln.startswith(">"):
            return ln[1:].strip()
    return ""


# ==========================================================================
# PART 1 - ARTIFICIAL SEQUENCE WITH TRANSPOSONS
# ==========================================================================
SCENARIOS = {
    "separate":  "4 TEs, no intersections (baseline)",
    "nested2":   "TE2 inserted inside TE1 (2 intersect) + 2 separate TEs",
    "nested3":   "TE3 inside TE2 inside TE1 (3 intersect) + 1 separate TE",
    "overlap2":  "TE1 and TE2 partially overlap (interleaved IRs) + 1 separate TE",
    "mixed3":    "TE1/TE2 overlap and TE3 is inserted inside both (3 intersect)",
}


def random_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_te(rng, name, length, ir_len=10):
    """TE = left IR + body + reverse complement of left IR."""
    while True:
        ir = random_dna(rng, ir_len)
        if ir != revcomp(ir):            # avoid palindromic IRs (left == right, ambiguous)
            break
    body = random_dna(rng, length - 2 * ir_len)
    return {"name": name, "ir": ir, "seq": ir + body + revcomp(ir)}


def insert_te(seq, ann, pos, te, tsd_len=4):
    """
    Real transposition: the target site seq[pos:pos+tsd] gets duplicated
        before:  ...xxxx[TTAG]yyyy...
        after:   ...xxxx[TTAG]==== TE ====[TTAG]yyyy...
    All existing annotations to the right of the insertion are shifted.
    """
    cut = pos + tsd_len
    add = len(te["seq"]) + tsd_len
    for a in ann:                              # shift coordinates (0-based, inclusive)
        if a["start"] >= cut:
            a["start"] += add
        if a["end"] >= cut:
            a["end"] += add
    new = seq[:cut] + te["seq"] + seq[pos:cut] + seq[cut:]
    ann.append({"name": te["name"], "ir": te["ir"], "orig": te["seq"],
                "start": cut, "end": cut + len(te["seq"]) - 1,
                "tsd": seq[pos:cut] if tsd_len else ""})
    return new


def insert_overlapping(seq, ann, pos, rng, nameA, nameB, lenA, lenB, overlap, ir_len=10):
    """
    Two TEs whose ends are interleaved:   L_A ..... L_B ..... R_A ..... R_B
                                          |------- TE A ------|
                                                    |------- TE B -------|
    (cannot happen by a single insertion - it's artificial, to stress the algorithm)
    """
    irA = make_te(rng, nameA, 2 * ir_len + 1, ir_len)["ir"]
    irB = make_te(rng, nameB, 2 * ir_len + 1, ir_len)["ir"]
    a1 = random_dna(rng, lenA - overlap - ir_len)          # part of A before L_B
    shared = random_dna(rng, overlap - 2 * ir_len)         # between L_B and R_A
    b2 = random_dna(rng, lenB - overlap - ir_len)          # part of B after R_A
    block = irA + a1 + irB + shared + revcomp(irA) + b2 + revcomp(irB)
    for a in ann:
        if a["start"] >= pos:
            a["start"] += len(block)
        if a["end"] >= pos:
            a["end"] += len(block)
    sA = pos
    eA = pos + len(irA + a1 + irB + shared + revcomp(irA)) - 1
    sB = pos + len(irA + a1)
    eB = pos + len(block) - 1
    ann.append({"name": nameA, "ir": irA, "orig": block[:eA - sA + 1], "start": sA, "end": eA, "tsd": ""})
    ann.append({"name": nameB, "ir": irB, "orig": block[sB - pos:], "start": sB, "end": eB, "tsd": ""})
    return seq[:pos] + block + seq[pos:]


def inside_body(a, rng, ir_len=10, margin=3):
    """Random insertion point inside a TE body (not touching its IRs)."""
    return rng.randint(a["start"] + ir_len + margin, a["end"] - ir_len - margin - 4)


def generate_sequence(scenario="nested3", seed=None, ir_len=10):
    """Returns (sequence, annotations). Retries until IRs are unique and length is 200-400."""
    for attempt in range(1000):
        rng = random.Random(None if seed is None else seed + attempt * 7919)
        ann = []
        if scenario == "separate":
            seq = random_dna(rng, 110)
            for i, L in enumerate([55, 60, 50, 65], 1):
                # insert left-to-right at free spots of the background
                pos = rng.randint(5, 20) + (ann[-1]["end"] + 5 if ann else 0)
                seq = insert_te(seq, ann, min(pos, len(seq) - 5), make_te(rng, f"TE{i}", L, ir_len))
        elif scenario == "nested2":
            seq = random_dna(rng, 120)
            seq = insert_te(seq, ann, 20, make_te(rng, "TE1", 90, ir_len))
            seq = insert_te(seq, ann, inside_body(ann[0], rng), make_te(rng, "TE2", 45, ir_len))
            seq = insert_te(seq, ann, ann[0]["end"] + 15, make_te(rng, "TE3", 45, ir_len))
            seq = insert_te(seq, ann, len(seq) - 20, make_te(rng, "TE4", 40, ir_len))
        elif scenario == "nested3":
            seq = random_dna(rng, 110)
            seq = insert_te(seq, ann, 15, make_te(rng, "TE1", 95, ir_len))
            seq = insert_te(seq, ann, inside_body(ann[0], rng), make_te(rng, "TE2", 70, ir_len))
            seq = insert_te(seq, ann, inside_body(ann[1], rng), make_te(rng, "TE3", 40, ir_len))
            seq = insert_te(seq, ann, len(seq) - 20, make_te(rng, "TE4", 50, ir_len))
        elif scenario == "overlap2":
            seq = random_dna(rng, 140)
            seq = insert_overlapping(seq, ann, 20, rng, "TE1", "TE2", 80, 80, 35, ir_len)
            seq = insert_te(seq, ann, len(seq) - 25, make_te(rng, "TE3", 55, ir_len))
        elif scenario == "mixed3":
            seq = random_dna(rng, 140)
            seq = insert_overlapping(seq, ann, 20, rng, "TE1", "TE2", 90, 90, 45, ir_len)
            # TE3 lands in the shared region -> inside TE1 AND TE2
            lo = ann[1]["start"] + ir_len + 3
            hi = ann[0]["end"] - ir_len - 7
            seq = insert_te(seq, ann, rng.randint(lo, hi), make_te(rng, "TE3", 40, ir_len))
        else:
            raise ValueError(f"unknown scenario {scenario}")

        if not 200 <= len(seq) <= 400:
            continue
        # every IR (and its reverse complement) must occur exactly once - no accidental copies
        ok = all(seq.count(a["ir"]) == 1 and seq.count(revcomp(a["ir"])) == 1 for a in ann)
        if ok:
            ann.sort(key=lambda a: a["start"])
            return seq, ann
    raise RuntimeError("could not generate a valid sequence")


def find_all(seq, pattern):
    """All (overlapping) occurrences of pattern."""
    out, i = [], seq.find(pattern)
    while i != -1:
        out.append(i)
        i = seq.find(pattern, i + 1)
    return out


def detect_known(seq, families, min_len=1, max_len=None):
    """
    PART 1 detection - inverted repeats KNOWN.
    families: {name: left_IR}. For each family find the left IR and its reverse complement
    (right IR), then pair them like brackets: each right IR closes the most recent open
    left IR of the same family (a stack) -> works for separate AND nested copies.
    Returns (elements, orphans).
    """
    found, orphans = [], []
    for name, ir in families.items():
        rc = revcomp(ir)
        events = [(p, "L") for p in find_all(seq, ir)] + [(p, "R") for p in find_all(seq, rc)]
        events.sort()
        stack = []
        for p, kind in events:
            if kind == "L":
                stack.append(p)
            elif stack:
                s = stack.pop()
                e = p + len(ir) - 1
                L = e - s + 1
                if L >= min_len and (max_len is None or L <= max_len):
                    found.append({"name": name, "start": s, "end": e, "length": L})
            else:
                orphans.append((name, "right IR without left IR", p))
        for s in stack:
            orphans.append((name, "left IR without right IR", s))
    found.sort(key=lambda d: d["start"])
    return found, orphans


def relation(a, b):
    """How two elements intersect."""
    if a["end"] < b["start"] or b["end"] < a["start"]:
        return None
    if a["start"] <= b["start"] and b["end"] <= a["end"]:
        return f"{b['name']} is nested inside {a['name']}"
    if b["start"] <= a["start"] and a["end"] <= b["end"]:
        return f"{a['name']} is nested inside {b['name']}"
    ov = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    return f"{a['name']} and {b['name']} overlap by {ov} bp"


def part1_report(seq, truth, found, orphans):
    out = [f"Sequence length: {len(seq)} bp\n",
           "GROUND TRUTH (what the generator inserted)  - 1-based coordinates",
           f"  {'TE':<5}{'start':>7}{'end':>7}{'len':>6}   left IR      TSD"]
    for a in truth:
        out.append(f"  {a['name']:<5}{a['start']+1:>7}{a['end']+1:>7}{a['end']-a['start']+1:>6}"
                   f"   {a['ir']}   {a['tsd'] or '-'}")
    out.append("\nDETECTED (known inverted repeats, bracket pairing)")
    out.append(f"  {'TE':<5}{'start':>7}{'end':>7}{'len':>6}   correct?  intact copy?")
    tmap = {a["name"]: a for a in truth}
    for f in found:
        t = tmap.get(f["name"])
        correct = "yes" if t and (t["start"], t["end"]) == (f["start"], f["end"]) else "NO"
        intact = "yes" if t and seq[f["start"]:f["end"] + 1] == t["orig"] else \
                 "no - interrupted by another TE"
        out.append(f"  {f['name']:<5}{f['start']+1:>7}{f['end']+1:>7}{f['length']:>6}   "
                   f"{correct:<9} {intact}")
    for o in orphans:
        out.append(f"  ! {o[0]}: {o[1]} at {o[2]+1}")
    missed = [a["name"] for a in truth if a["name"] not in {f["name"] for f in found}]
    if missed:
        out.append(f"  MISSED: {', '.join(missed)}")

    out.append("\nINTERSECTIONS")
    rels = [relation(a, b) for i, a in enumerate(found) for b in found[i + 1:]]
    rels = [r for r in rels if r]
    out += [f"  - {r}" for r in rels] or ["  none"]
    out.append("\nObservation: a search for the COMPLETE element sequence fails for any TE that"
               "\nhas another TE inside it; the inverted-repeat search still finds both ends.")
    return "\n".join(out)


# ==========================================================================
# PART 2 - DE NOVO SEARCH (inverted repeats UNKNOWN)
# ==========================================================================
ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def find_inverted_repeats_naive(seq, k=10, min_gap=100, max_gap=None):
    """
    The literal s1/s2 algorithm from the lab sheet. O(n^2) - only for short sequences!
      s1 = seq[i:i+k] slides with step 1
      s2 slides from (end of s1 + min_gap) to the end, looking for revcomp(s1)
    Returns seed pairs (i, j): seq[j:j+k] == revcomp(seq[i:i+k]).
    """
    n, seeds = len(seq), []
    for i in range(n - k + 1):
        s1 = seq[i:i + k]
        target = revcomp(s1)                       # s1 inverted + complemented
        j_end = n - k if max_gap is None else min(n - k, i + k + max_gap)
        for j in range(i + k + min_gap, j_end + 1):
            if seq[j:j + k] == target:
                seeds.append((i, j))
    return seeds


def find_seeds_indexed(seq, k=12, min_gap=100, max_gap=10000, progress=None):
    """
    Same result as the naive version, but s2 does not scan base by base:
    every k-mer position is stored in an index (a linked list per k-mer), so
    'where does revcomp(s1) occur downstream?' is answered by a lookup.
    k-mers are encoded as integers, 2 bits per base (A=0 C=1 G=2 T=3).
    """
    n = len(seq)
    m = n - k + 1
    if m <= 0:
        return []
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    fwd = array("i", [-1]) * m      # code of seq[i:i+k]
    rev = array("i", [-1]) * m      # code of revcomp(seq[i:i+k])
    code = rc = valid = 0
    for idx, ch in enumerate(seq):
        b = ENC.get(ch, -1)
        if b < 0:                    # N or other symbol -> restart the window
            valid = code = rc = 0
            continue
        code = ((code << 2) | b) & mask
        rc = (rc >> 2) | ((3 - b) << shift)
        valid += 1
        if valid >= k:
            p = idx - k + 1
            fwd[p] = code
            rev[p] = rc
        if progress and idx % 500000 == 0:
            progress(0.3 * idx / n)

    # head[code] = first position of that k-mer, nxt[p] = next position with the same k-mer
    head = array("i", [-1]) * (1 << (2 * k))
    nxt = array("i", [-1]) * m
    for p in range(m - 1, -1, -1):
        c = fwd[p]
        if c >= 0:
            nxt[p] = head[c]
            head[c] = p

    seeds = []
    for i in range(m):               # s1 slides with step 1
        target = rev[i]
        if target < 0:
            continue
        lo = i + k + min_gap         # s2 must start at least min_gap after s1 ends
        hi = n - k if max_gap is None else i + k + max_gap
        p = head[target]
        while p != -1 and p < lo:    # skip occurrences that are too close (lo only grows,
            p = nxt[p]               # so we can permanently advance the pointer)
        head[target] = p
        while p != -1 and p <= hi:
            seeds.append((i, p))
            p = nxt[p]
        if progress and i % 500000 == 0:
            progress(0.3 + 0.65 * i / m)
    return seeds


def merge_seeds(seq, seeds, k, min_ir=12, max_join=3):
    """
    Consecutive seeds of the same inverted repeat lie on one 'anti-diagonal':
    (i, j) -> (i+1, j-1) -> ... so i + j is constant. Seeds with the same i+j and
    close i are glued into one longer IR (gaps up to max_join bases = mismatches).
    """
    open_runs, runs = {}, []
    for i, j in seeds:                 # seeds come sorted by i
        s = i + j
        r = open_runs.get(s)
        if r and i - r[1] <= max_join + 1:
            r[1] = i
        else:
            if r:
                runs.append((s, r[0], r[1]))
            open_runs[s] = [i, i]
    runs += [(s, r[0], r[1]) for s, r in open_runs.items()]

    out = []
    n = len(seq)
    for s, i0, i1 in runs:
        L = i1 - i0 + k                         # IR arm length
        if L < min_ir:
            continue
        left = (i0, i0 + L)                     # [start, end) of left arm
        right = (s - i1, s - i0 + k)            # [start, end) of right arm
        left_arm = seq[left[0]:left[1]]
        right_arm = seq[right[0]:right[1]]
        mism = sum(1 for a, b in zip(left_arm, revcomp(right_arm)) if a != b)
        top = max(left_arm.count(c) for c in "ACGT") / L
        start, end = left[0], right[1] - 1      # whole element, 0-based inclusive
        tsd = ""
        for t in range(12, 2, -1):              # target site duplication flanking the element?
            if start - t >= 0 and end + 1 + t <= n and seq[start - t:start] == seq[end + 1:end + 1 + t]:
                tsd = seq[start - t:start]
                break
        out.append({"start": start, "end": end, "length": end - start + 1,
                    "ir_len": L, "mismatches": mism, "left_ir": left_arm,
                    "right_ir": right_arm, "gap": right[0] - left[1],
                    "tsd": tsd, "low_complexity": top > 0.7})
    out.sort(key=lambda d: (-d["ir_len"], d["mismatches"], d["start"]))
    return out


def find_transposon_candidates(seq, k=12, min_gap=100, max_gap=10000, min_ir=16,
                               max_join=3, skip_low_complexity=True, progress=None):
    t0 = time.time()
    seeds = find_seeds_indexed(seq, k, min_gap, max_gap, progress)
    cands = merge_seeds(seq, seeds, k, min_ir, max_join)
    if skip_low_complexity:
        cands = [c for c in cands if not c["low_complexity"]]
    # how many elements in the genome end with this IR? (each element = 1 left + 1 right arm,
    # on either strand). Real transposon families usually have several copies.
    for c in cands:
        probe = c["left_ir"][:min(20, c["ir_len"])]
        c["copies"] = max(1, (seq.count(probe) + seq.count(revcomp(probe))) // 2)
    if progress:
        progress(1.0)
    return cands, len(seeds), time.time() - t0


def candidates_report(cands, limit=100):
    out = [f"{'#':>4} {'start':>10} {'end':>10} {'length':>8} {'IR':>4} {'mm':>3} {'copies':>6} {'TSD':<12} left IR"]
    for i, c in enumerate(cands[:limit], 1):
        out.append(f"{i:>4} {c['start']+1:>10,} {c['end']+1:>10,} {c['length']:>8,} "
                   f"{c['ir_len']:>4} {c['mismatches']:>3} {c.get('copies', 1):>6} {c['tsd'] or '-':<12} {c['left_ir'][:40]}")
    if len(cands) > limit:
        out.append(f"... {len(cands) - limit:,} more")
    return "\n".join(out)


# ==========================================================================
# CONSOLE MODE
# ==========================================================================
def cli_part1(scenario="nested3", seed=1):
    seq, truth = generate_sequence(scenario, seed)
    print(f"Scenario: {scenario} - {SCENARIOS[scenario]}")
    print(">artificial\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    found, orphans = detect_known(seq, {a["name"]: a["ir"] for a in truth})
    print(part1_report(seq, truth, found, orphans))


def cli_genome(path, k=12, min_gap=100, max_gap=10000, min_ir=16):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    text, nrec = first_record(text)
    if nrec > 1:
        print(f"(file has {nrec} records - using the first one)")
    seq = clean_sequence(text)
    print(f"{fasta_header(text)} | {len(seq):,} bp | k={k} min_gap={min_gap} max_gap={max_gap} min_IR={min_ir}")
    cands, nseeds, secs = find_transposon_candidates(
        seq, k, min_gap, max_gap, min_ir,
        progress=lambda x: print(f"\r  {100*x:5.1f}%", end="", flush=True))
    print(f"\n{nseeds:,} seed pairs -> {len(cands):,} candidate elements in {secs:.1f} s\n")
    print(candidates_report(cands, 40))


# ==========================================================================
# GUI
# ==========================================================================
def run_gui():
    import threading
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox, scrolledtext

    COLORS = ["#e6550d", "#3182bd", "#31a354", "#756bb1", "#d6616b", "#8c6d31"]

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("Transposable Element Finder")
            self.geometry("1280x820")
            self.seq1, self.truth, self.found = "", [], []
            self.genome, self.genome_name, self.cands = "", "", []
            nb = ttk.Notebook(self)
            nb.pack(fill="both", expand=True)
            self.nb = nb
            self._build_part1(ttk.Frame(nb, padding=6))
            self._build_part2(ttk.Frame(nb, padding=6))

        # ---------------- PART 1 ----------------
        def _build_part1(self, f):
            self.nb.add(f, text="1. Artificial sequence (known IRs)")
            top = ttk.Frame(f)
            top.pack(fill="x")
            ttk.Label(top, text="Scenario:").pack(side="left")
            self.scenario = tk.StringVar(value="nested3")
            cb = ttk.Combobox(top, textvariable=self.scenario, values=list(SCENARIOS),
                              state="readonly", width=10)
            cb.pack(side="left")
            self.scen_desc = ttk.Label(top, text=SCENARIOS["nested3"], foreground="gray")
            cb.bind("<<ComboboxSelected>>",
                    lambda e: self.scen_desc.config(text=SCENARIOS[self.scenario.get()]))
            ttk.Label(top, text="  seed:").pack(side="left")
            self.seed = ttk.Entry(top, width=6)
            self.seed.insert(0, "1")
            self.seed.pack(side="left")
            ttk.Button(top, text="Generate", command=self.generate).pack(side="left", padx=4)
            ttk.Button(top, text="Detect (known IRs)", command=self.detect1).pack(side="left")
            ttk.Button(top, text="Save FASTA...", command=self.save_fasta1).pack(side="left", padx=4)
            ttk.Button(top, text="Search it with Part 2 ->", command=self.send_to_part2).pack(side="left")
            self.scen_desc.pack(side="left", padx=8)

            self.seq_box = tk.Text(f, height=8, font=("Courier", 10), wrap="char")
            self.seq_box.pack(fill="x", pady=4)
            self.map1 = tk.Canvas(f, height=170, bg="white", highlightthickness=1,
                                  highlightbackground="#ccc")
            self.map1.pack(fill="x")
            self.map1.bind("<Configure>", lambda e: self.draw_map1())
            self.out1 = scrolledtext.ScrolledText(f, font=("Courier", 9))
            self.out1.pack(fill="both", expand=True, pady=4)

        def generate(self):
            try:
                seed = int(self.seed.get())
            except ValueError:
                seed = None
            self.seq1, self.truth = generate_sequence(self.scenario.get(), seed)
            self.scen_desc.config(text=SCENARIOS[self.scenario.get()])
            self.found = []
            self.seq_box.delete("1.0", "end")
            self.seq_box.insert("1.0", self.seq1)
            for t in self.seq_box.tag_names():
                self.seq_box.tag_delete(t)
            for i, a in enumerate(self.truth):          # colour the IRs in the text
                col = COLORS[i % len(COLORS)]
                tag = f"ir{i}"
                self.seq_box.tag_config(tag, background=col, foreground="white")
                L = len(a["ir"])
                for s in (a["start"], a["end"] - L + 1):
                    self.seq_box.tag_add(tag, f"1.0+{s}c", f"1.0+{s + L}c")
            self.out1.delete("1.0", "end")
            self.out1.insert("end", f"Generated {len(self.seq1)} bp - {SCENARIOS[self.scenario.get()]}\n"
                                    "Coloured = inverted repeats. Press 'Detect'.\n")
            self.draw_map1()

        def detect1(self):
            if not self.seq1:
                self.generate()
            # the user may have edited the sequence in the box - always read it back
            self.seq1 = clean_sequence(self.seq_box.get("1.0", "end"))
            fam = {a["name"]: a["ir"] for a in self.truth}
            self.found, orphans = detect_known(self.seq1, fam)
            self.out1.delete("1.0", "end")
            self.out1.insert("end", part1_report(self.seq1, self.truth, self.found, orphans))
            self.draw_map1()

        def draw_map1(self):
            c = self.map1
            c.delete("all")
            if not self.seq1:
                return
            W = c.winfo_width()
            n = len(self.seq1)
            x0, x1 = 60, W - 20
            X = lambda p: x0 + (x1 - x0) * p / n
            c.create_line(x0, 20, x1, 20, width=3, fill="#555")
            for p in range(0, n + 1, 50):
                c.create_line(X(p), 16, X(p), 24, fill="#555")
                c.create_text(X(p), 10, text=str(p), font=("Arial", 7))

            def bars(items, y0, label, colmap):
                c.create_text(5, y0 + 20, text=label, anchor="w", font=("Arial", 8, "bold"))
                levels = []                              # stack overlapping bars on levels
                for a in sorted(items, key=lambda d: (d["start"], -d["end"])):
                    lvl = next((li for li, e in enumerate(levels) if e < a["start"]), len(levels))
                    if lvl == len(levels):
                        levels.append(-1)
                    levels[lvl] = a["end"]
                    y = y0 + lvl * 15
                    col = colmap.get(a["name"], "#999")
                    c.create_rectangle(X(a["start"]), y, X(a["end"] + 1), y + 11, fill=col, outline="")
                    c.create_text(X(a["start"]) + 3, y + 5, text=a["name"], anchor="w",
                                  fill="white", font=("Arial", 7, "bold"))

            colmap = {a["name"]: COLORS[i % len(COLORS)] for i, a in enumerate(self.truth)}
            bars(self.truth, 32, "truth", colmap)
            if self.found:
                bars(self.found, 100, "found", colmap)

        def save_fasta1(self):
            if not self.seq1:
                return
            path = filedialog.asksaveasfilename(defaultextension=".fasta")
            if path:
                with open(path, "w") as fh:
                    fh.write(f">artificial_{self.scenario.get()}\n")
                    fh.write("\n".join(self.seq1[i:i + 70] for i in range(0, len(self.seq1), 70)) + "\n")

        def send_to_part2(self):
            if not self.seq1:
                self.generate()
            self._set_genome(">artificial sequence from part 1\n" + self.seq1, "artificial")
            self.k.set(8)
            self.min_gap.set(10)
            self.max_gap.set(400)
            self.min_ir.set(9)
            self.nb.select(1)
            messagebox.showinfo("Part 2", "Loaded the artificial sequence. Parameters lowered "
                                          "(k=8, min gap=10, min IR=9) because these TEs are short.\n"
                                          "Press 'Find candidates'.")

        # ---------------- PART 2 ----------------
        def _build_part2(self, f):
            self.nb.add(f, text="2. Genome search (unknown IRs)")
            row = ttk.Frame(f)
            row.pack(fill="x")
            ttk.Button(row, text="Load genome FASTA...", command=self.load_genome).pack(side="left")
            self.gen_info = ttk.Label(row, text="  no genome loaded", foreground="gray")
            self.gen_info.pack(side="left")

            par = ttk.LabelFrame(f, text="Parameters", padding=6)
            par.pack(fill="x", pady=4)
            self.k = tk.IntVar(value=12)
            self.min_gap = tk.IntVar(value=100)
            self.max_gap = tk.IntVar(value=10000)
            self.min_ir = tk.IntVar(value=16)
            self.max_join = tk.IntVar(value=3)
            self.skip_lc = tk.BooleanVar(value=True)
            fields = [("s1 window k (6-12)", self.k), ("min gap s1->s2", self.min_gap),
                      ("max gap (0 = to genome end)", self.max_gap),
                      ("min IR length", self.min_ir), ("mismatch gap join", self.max_join)]
            for col, (lbl, var) in enumerate(fields):
                ttk.Label(par, text=lbl).grid(row=0, column=col, padx=6, sticky="w")
                ttk.Entry(par, textvariable=var, width=8).grid(row=1, column=col, padx=6, sticky="w")
            ttk.Checkbutton(par, text="skip low-complexity (AAAA..)", variable=self.skip_lc).grid(
                row=1, column=len(fields), padx=6)
            ttk.Button(par, text="Find candidates", command=self.run_search).grid(
                row=1, column=len(fields) + 1, padx=6)
            ttk.Button(par, text="Export CSV...", command=self.export).grid(
                row=1, column=len(fields) + 2)
            self.pbar = ttk.Progressbar(f, maximum=1.0)
            self.pbar.pack(fill="x")
            self.status = ttk.Label(f, text="")
            self.status.pack(anchor="w")

            cols = ("start", "end", "length", "ir", "mm", "copies", "tsd", "left", "right")
            heads = ("start", "end", "length", "IR len", "mism.", "IR copies", "TSD", "left IR (s1)", "right IR (s2)")
            widths = (90, 90, 70, 55, 50, 70, 90, 250, 250)
            frame = ttk.Frame(f)
            frame.pack(fill="both", expand=True)
            self.tree = ttk.Treeview(frame, columns=cols, show="headings")
            for c, h, w in zip(cols, heads, widths):
                self.tree.heading(c, text=h, command=lambda c=c: self.sort_tree(c))
                self.tree.column(c, width=w, anchor="e" if w < 100 else "w")
            sb = ttk.Scrollbar(frame, command=self.tree.yview)
            self.tree.configure(yscrollcommand=sb.set)
            self.tree.pack(side="left", fill="both", expand=True)
            sb.pack(side="left", fill="y")
            self.map2 = tk.Canvas(f, height=80, bg="white", highlightthickness=1,
                                  highlightbackground="#ccc")
            self.map2.pack(fill="x", pady=4)
            self.map2.bind("<Configure>", lambda e: self.draw_map2())

        def _set_genome(self, text, name):
            text, nrec = first_record(text)
            if nrec > 1:
                messagebox.showinfo("FASTA", f"The file has {nrec} records - using the first one "
                                             "(normally the chromosome).")
            self.genome = clean_sequence(text)
            self.genome_name = fasta_header(text) or name
            gc = 100 * (self.genome.count("G") + self.genome.count("C")) / max(1, len(self.genome))
            self.gen_info.config(text=f"  {self.genome_name[:80]} | {len(self.genome):,} bp | GC {gc:.1f}%",
                                 foreground="black")

        def load_genome(self):
            path = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta *.fa *.fna *.txt"), ("All", "*.*")])
            if path:
                with open(path, encoding="utf-8", errors="replace") as fh:
                    self._set_genome(fh.read(), path.split("/")[-1])

        def run_search(self):
            if not self.genome:
                messagebox.showwarning("Search", "Load a genome FASTA first.")
                return
            k = self.k.get()
            if not 6 <= k <= 12:
                messagebox.showwarning("Search", "k must be 6..12 (the index has 4^k slots).")
                return
            args = dict(k=k, min_gap=self.min_gap.get(),
                        max_gap=self.max_gap.get() or None, min_ir=max(k, self.min_ir.get()),
                        max_join=self.max_join.get(), skip_low_complexity=self.skip_lc.get())
            self.status.config(text="searching...")
            self.config(cursor="watch")

            def prog(x):
                self.after(0, lambda: self.pbar.config(value=x))

            def worker():
                try:
                    res = find_transposon_candidates(self.genome, progress=prog, **args)
                    self.after(0, lambda: self.show_results(*res))
                except MemoryError:
                    self.after(0, lambda: messagebox.showerror("Search", "Out of memory - lower max gap."))
                finally:
                    self.after(0, lambda: self.config(cursor=""))
            threading.Thread(target=worker, daemon=True).start()

        def show_results(self, cands, nseeds, secs):
            self.cands = cands
            self.tree.delete(*self.tree.get_children())
            for c in cands[:5000]:
                self.tree.insert("", "end", values=(f"{c['start']+1:,}", f"{c['end']+1:,}", f"{c['length']:,}",
                                                    c["ir_len"], c["mismatches"], c["copies"], c["tsd"] or "-",
                                                    c["left_ir"], c["right_ir"]))
            self.status.config(text=f"{nseeds:,} seed pairs (s1/s2 hits) -> {len(cands):,} candidate "
                                    f"elements in {secs:.1f} s" +
                                    ("  (showing the first 5,000; Export for all)" if len(cands) > 5000 else ""))
            self.draw_map2()

        def sort_tree(self, col):
            rows = [(self.tree.set(i, col), i) for i in self.tree.get_children("")]
            def key(v):
                try:
                    return float(v[0].replace(",", ""))
                except ValueError:
                    return v[0]
            rows.sort(key=key, reverse=getattr(self, "_rev", False))
            self._rev = not getattr(self, "_rev", False)
            for idx, (_, i) in enumerate(rows):
                self.tree.move(i, "", idx)

        def draw_map2(self):
            c = self.map2
            c.delete("all")
            if not self.cands or not self.genome:
                return
            W, n = c.winfo_width(), len(self.genome)
            x0, x1 = 20, W - 20
            c.create_line(x0, 40, x1, 40, width=3, fill="#555")
            c.create_text(x0, 70, text="1", anchor="w", font=("Arial", 7))
            c.create_text(x1, 70, text=f"{n:,} bp", anchor="e", font=("Arial", 7))
            best = self.cands[0]["ir_len"]
            for cd in self.cands[:2000]:
                x = x0 + (x1 - x0) * cd["start"] / n
                h = 6 + 24 * cd["ir_len"] / best
                col = "#d62728" if cd["tsd"] else "#1f77b4"
                c.create_line(x, 40 - h, x, 40, fill=col)
            c.create_text(W / 2, 70, text="candidate positions (bar height = IR length, red = has TSD)",
                          font=("Arial", 7), fill="#555")

        def export(self):
            if not self.cands:
                return
            path = filedialog.asksaveasfilename(defaultextension=".csv")
            if path:
                with open(path, "w") as fh:
                    fh.write("start,end,length,ir_len,mismatches,ir_copies,gap,tsd,left_ir,right_ir\n")
                    for c in self.cands:
                        fh.write(f"{c['start']+1},{c['end']+1},{c['length']},{c['ir_len']},{c['mismatches']},"
                                 f"{c['copies']},{c['gap']},{c['tsd']},{c['left_ir']},{c['right_ir']}\n")

    App().mainloop()


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--cli-part1":
        cli_part1(a[1] if len(a) > 1 else "nested3", int(a[2]) if len(a) > 2 else 1)
    elif a and a[0] == "--cli-genome":
        nums = [int(x) for x in a[2:6]]
        cli_genome(a[1], *nums)
    else:
        run_gui()