import tkinter as tk
import tkinter.font as tkfont
from dataclasses import dataclass

TSD_LENGTH = 5
IR_MIN = 8
ELEMENT_MIN = 30
LINE_WIDTH = 60
COLORS = ["#ffc9c9", "#c9e4ff", "#cdf2c6", "#ffeaa8", "#e3cfff"]
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

HOST = ("CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACC"
        "TGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCT"
        "CATGTAGCCAGAAGGCTGC")

TRANSPOSONS = {
    "T1": ("GATTACAGCA", "AACTCATCGACTCTATGTAGTGACCGCGTC"),
    "T2": ("TTCAGGCATG", "GATGTCAAACCCCGGGGGGA"),
    "T3": ("AGGTCCTTAC", "GCTCAGATATCCGATACAGGGATGAA"),
    "T4": ("CCATTGGACA", "GAAATAACCTCATCCCATTGGT"),
}

SEPARATE = [("T1", "host", 30), ("T2", "host", 70), ("T3", "host", 110), ("T4", "host", 150)]
INTERSECTING = [("T1", "host", 30), ("T2", "T1", 22), ("T3", "host", 120), ("T4", "T3", 3)]

INTRO = """Exercise 1 shows the artificial DNA sequence written in this program: a host DNA with
4 transposons. Tick "Intersecting transposons" to use the same transposons, but inserted
into each other.
Exercise 2 finds the transposons again, using only the DNA sequence.
"""

STRUCTURE = f"""Each transposon = left inverted repeat (IR) + inner part + right IR (reverse complement
of the left IR). It is inserted after a {TSD_LENGTH} bp target site, which is copied on its other
side, so it ends up between two direct repeats (TSD = target site duplication).
"""


@dataclass
class Transposon:
    name: str
    start: int
    end: int
    ir: str = ""


def reverse_complement(dna: str) -> str:
    return dna.translate(COMPLEMENT)[::-1]


def insert_transposon(sequence: str, transposons: list[Transposon], name: str, cut: int) -> str:
    ir, inner = TRANSPOSONS[name]
    element = ir + inner + reverse_complement(ir)
    shift = len(element) + TSD_LENGTH
    for t in transposons:
        if t.start >= cut:
            t.start += shift
        if t.end > cut:
            t.end += shift
    transposons.append(Transposon(name, cut, cut + len(element), ir))
    return sequence[:cut] + element + sequence[cut - TSD_LENGTH:cut] + sequence[cut:]


def build(layout: list[tuple[str, str, int]]) -> tuple[str, list[Transposon]]:
    sequence = HOST
    transposons: list[Transposon] = []
    into_host = sorted((step for step in layout if step[1] == "host"), key=lambda step: -step[2])
    into_transposons = [step for step in layout if step[1] != "host"]
    for name, target, offset in into_host + into_transposons:
        cut = offset if target == "host" else next(t.start for t in transposons if t.name == target) + offset
        sequence = insert_transposon(sequence, transposons, name, cut)
    transposons.sort(key=lambda t: (t.start, -t.end))
    return sequence, transposons


def recipe(layout: list[tuple[str, str, int]]) -> list[str]:
    return [f"  {name} inserted " + (f"after base {offset} of the host DNA" if target == "host"
                                     else f"{offset} bases into {target}")
            for name, target, offset in layout]


def find_transposons(sequence: str) -> list[Transposon]:
    found = []
    for start in range(TSD_LENGTH, len(sequence) - ELEMENT_MIN - TSD_LENGTH + 1):
        right_ir = reverse_complement(sequence[start:start + IR_MIN])
        hit = sequence.find(right_ir, start + ELEMENT_MIN - IR_MIN)
        while hit != -1:
            end = hit + IR_MIN
            if sequence[start - TSD_LENGTH:start] == sequence[end:end + TSD_LENGTH]:
                found.append(Transposon("", start, end))
            hit = sequence.find(right_ir, hit + 1)
    found.sort(key=lambda t: (t.start, -t.end))
    for number, t in enumerate(found, 1):
        t.name = f"D{number}"
    return found


def inverted_repeat_length(sequence: str, t: Transposon) -> int:
    length = IR_MIN
    while (t.start + length < t.end - length - 1
           and sequence[t.start + length] == reverse_complement(sequence[t.end - length - 1])):
        length += 1
    return length


def parent_of(t: Transposon, transposons: list[Transposon]) -> str:
    parents = [p for p in transposons if p.start < t.start and t.end < p.end]
    return min(parents, key=lambda p: p.end - p.start).name if parents else "-"


def why_missed(sequence: str, t: Transposon) -> str:
    ir_ok = sequence[t.start:t.start + IR_MIN] == reverse_complement(sequence[t.end - IR_MIN:t.end])
    tsd_ok = sequence[t.start - TSD_LENGTH:t.start] == sequence[t.end:t.end + TSD_LENGTH]
    if not ir_ok and not tsd_ok:
        return "its ends are no longer inverted repeats and its TSDs differ"
    if not ir_ok:
        return "another transposon cut into one of its inverted repeats"
    if not tsd_ok:
        return "another transposon cut into one of its TSDs"
    return "its signature is intact, but it was not reported"


def match_name(d: Transposon, inserted: list[Transposon]) -> str:
    for t in inserted:
        if (t.start, t.end) == (d.start, d.end):
            return t.name
    for t in inserted:
        if t.start <= d.start and d.end <= t.end and (t.start == d.start or t.end == d.end):
            return f"part of {t.name}"
    return "no"


def describe_extra(d: Transposon, inserted: list[Transposon]) -> str:
    for t in inserted:
        if t.start <= d.start and d.end <= t.end and t.end == d.end:
            return f"is the part of {t.name} after the place where another transposon cut it"
        if t.start <= d.start and d.end <= t.end and t.start == d.start:
            return f"is the part of {t.name} before the place where another transposon cut it"
    return "was not inserted, this part of the sequence matches the rule by chance"


def inserted_table(sequence: str, transposons: list[Transposon]) -> list[str]:
    lines = [f"  {'#':<4}{'start':>6}{'end':>6}{'length':>8}   {'left IR':<14}{'right IR':<14}{'TSD':<8}inside"]
    for t in transposons:
        lines.append(f"  {t.name:<4}{t.start + 1:>6}{t.end:>6}{t.end - t.start:>8}   {t.ir:<14}"
                     f"{reverse_complement(t.ir):<14}{sequence[t.start - TSD_LENGTH:t.start]:<8}"
                     f"{parent_of(t, transposons)}")
    return lines


def detected_table(sequence: str, found: list[Transposon], inserted: list[Transposon]) -> list[str]:
    lines = [f"  {'#':<4}{'start':>6}{'end':>6}{'length':>8}   {'left IR':<14}{'right IR':<14}{'TSD':<8}"
             f"{'inside':<8}inserted as"]
    for d in found:
        size = min(inverted_repeat_length(sequence, d), 13)
        lines.append(f"  {d.name:<4}{d.start + 1:>6}{d.end:>6}{d.end - d.start:>8}   "
                     f"{sequence[d.start:d.start + size]:<14}{sequence[d.end - size:d.end]:<14}"
                     f"{sequence[d.start - TSD_LENGTH:d.start]:<8}{parent_of(d, found):<8}"
                     f"{match_name(d, inserted)}")
    return lines


def comparison(sequence: str, found: list[Transposon], inserted: list[Transposon]) -> list[str]:
    detected = {(d.start, d.end) for d in found}
    missed = [t for t in inserted if (t.start, t.end) not in detected]
    extra = [d for d in found if (d.start, d.end) not in {(t.start, t.end) for t in inserted}]
    lines = [f"Found {len(inserted) - len(missed)} of {len(inserted)} inserted transposons"
             + (f", plus {len(extra)} other detection(s)." if extra else ", no false detections.")]
    lines += [f"  {t.name} ({t.start + 1}-{t.end}) was not found: {why_missed(sequence, t)}." for t in missed]
    lines += [f"  {d.name} ({d.start + 1}-{d.end}) {describe_extra(d, inserted)}." for d in extra]
    return lines


class App:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.layout = SEPARATE
        self.sequence = ""
        self.inserted: list[Transposon] = []
        self.view = 0
        root.title("Transposable elements")

        top = tk.Frame(root)
        top.pack(fill="x", padx=8, pady=(8, 4))
        tk.Button(top, text="Exercise 1 - Artificial sequence", command=self.exercise1).pack(side="left")
        self.intersect = tk.BooleanVar(value=False)
        tk.Checkbutton(top, text="Intersecting transposons", variable=self.intersect,
                       command=self.refresh).pack(side="left", padx=(4, 16))
        tk.Button(top, text="Exercise 2 - Detect transposons", command=self.exercise2).pack(side="left")
        self.info = tk.Label(root, anchor="w")
        self.info.pack(fill="x", padx=8)

        main = tk.Frame(root)
        main.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        scrollbar = tk.Scrollbar(main)
        scrollbar.pack(side="right", fill="y")
        family = "Consolas" if "Consolas" in tkfont.families() else tkfont.nametofont("TkFixedFont").actual("family")
        self.output = tk.Text(main, width=92, height=38, font=(family, 10), wrap="none", state="disabled",
                              yscrollcommand=scrollbar.set)
        self.output.pack(side="left", fill="both", expand=True)
        scrollbar.configure(command=self.output.yview)
        self.output.tag_configure("tsd", underline=True)
        for number, color in enumerate(COLORS):
            self.output.tag_configure(f"t{number}", background=color)
        self.output.tag_raise("tsd")
        self.select()
        self.write(INTRO)

    def select(self) -> None:
        self.layout = INTERSECTING if self.intersect.get() else SEPARATE
        self.sequence, self.inserted = build(self.layout)
        kind = "intersecting" if self.intersect.get() else "separate"
        self.info.configure(text=f"Artificial sequence with {kind} transposons:  {len(self.sequence)} bp, "
                                 f"{len(self.inserted)} transposons")

    def refresh(self) -> None:
        if self.view == 1:
            self.exercise1()
        elif self.view == 2:
            self.exercise2()
        else:
            self.select()

    def clear(self) -> None:
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.configure(state="disabled")

    def write(self, text: str = "") -> None:
        self.output.configure(state="normal")
        self.output.insert("end", text + "\n")
        self.output.configure(state="disabled")

    def write_sequence(self, transposons: list[Transposon]) -> None:
        first_line = int(self.output.index("end-1c").split(".")[0])
        for line_start in range(0, len(self.sequence), LINE_WIDTH):
            chunk = self.sequence[line_start:line_start + LINE_WIDTH]
            self.write(f"{line_start + 1:>6} " + " ".join(chunk[i:i + 10] for i in range(0, len(chunk), 10)))
        self.output.configure(state="normal")
        for number, t in enumerate(transposons):
            self.mark(first_line, f"t{number % len(COLORS)}", t.start, t.end)
            self.mark(first_line, "tsd", t.start - TSD_LENGTH, t.start)
            self.mark(first_line, "tsd", t.end, t.end + TSD_LENGTH)
        self.output.configure(state="disabled")

    def mark(self, first_line: int, tag: str, start: int, end: int) -> None:
        for position in range(max(start, 0), min(end, len(self.sequence))):
            line = first_line + position // LINE_WIDTH
            column = 7 + position % LINE_WIDTH + position % LINE_WIDTH // 10
            joined = position + 1 < end and (position + 1) % LINE_WIDTH != 0 and (position + 1) % 10 == 0
            self.output.tag_add(tag, f"{line}.{column}", f"{line}.{column + (2 if joined else 1)}")

    def exercise1(self) -> None:
        self.select()
        self.view = 1
        self.clear()
        kind = "intersecting" if self.intersect.get() else "separate"
        self.write(f"EXERCISE 1 - Artificial DNA sequence with {kind} transposons ({len(self.sequence)} bp)")
        self.write()
        self.write(STRUCTURE)
        self.write(f"How it was made, starting from a host DNA of {len(HOST)} bp:")
        self.write("\n".join(recipe(self.layout)))
        self.write()
        self.write("Transposons in the final sequence (1-based positions, without the TSDs):")
        self.write("\n".join(inserted_table(self.sequence, self.inserted)))
        self.write()
        self.write("Sequence (transposons highlighted, TSDs underlined):")
        self.write_sequence(self.inserted)

    def exercise2(self) -> None:
        self.select()
        self.view = 2
        found = find_transposons(self.sequence)
        self.clear()
        kind = "intersecting" if self.intersect.get() else "separate"
        self.write(f"EXERCISE 2 - Detection in the sequence with {kind} transposons ({len(self.sequence)} bp)")
        self.write(f"Rule: the first {IR_MIN} bases of the element are the reverse complement of its last {IR_MIN} bases"
                   f"\n      (inverted repeats) and the {TSD_LENGTH} bases before it equal the {TSD_LENGTH} bases after it (TSD).")
        self.write()
        self.write("Inserted (Exercise 1):")
        self.write("\n".join(inserted_table(self.sequence, self.inserted)))
        self.write()
        self.write("Detected:")
        self.write("\n".join(detected_table(self.sequence, found, self.inserted)) if found else "  nothing")
        self.write()
        self.write("\n".join(comparison(self.sequence, found, self.inserted)))
        self.write()
        self.write("Sequence (detected transposons highlighted, TSDs underlined):")
        self.write_sequence(found)


def main() -> None:
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
