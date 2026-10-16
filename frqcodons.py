from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk


TOP_N = 10
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}
MIN_ORF_CODONS = 100


_BASES = "UCAG"
_RNA_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_AMINO_ACIDS_IN_TABLE_ORDER = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
GENETIC_CODE = dict(zip(_RNA_CODONS, _AMINO_ACIDS_IN_TABLE_ORDER))

AMINO_ACID_INFO = {
    "A": ("Ala", "Alanine"),        "R": ("Arg", "Arginine"),
    "N": ("Asn", "Asparagine"),     "D": ("Asp", "Aspartic acid"),
    "C": ("Cys", "Cysteine"),       "Q": ("Gln", "Glutamine"),
    "E": ("Glu", "Glutamic acid"),  "G": ("Gly", "Glycine"),
    "H": ("His", "Histidine"),      "I": ("Ile", "Isoleucine"),
    "L": ("Leu", "Leucine"),        "K": ("Lys", "Lysine"),
    "M": ("Met", "Methionine"),     "F": ("Phe", "Phenylalanine"),
    "P": ("Pro", "Proline"),        "S": ("Ser", "Serine"),
    "T": ("Thr", "Threonine"),      "W": ("Trp", "Tryptophan"),
    "Y": ("Tyr", "Tyrosine"),       "V": ("Val", "Valine"),
    "*": ("Stop", "Stop codon"),
}

VALID_DNA_CODONS = {codon.replace("U", "T") for codon in GENETIC_CODE}


def transcribe(dna_codon: str) -> str:
    return dna_codon.replace("T", "U")


def translate(rna_codon: str) -> tuple[str, str, str]:
    one_letter = GENETIC_CODE[rna_codon]
    three_letter, name = AMINO_ACID_INFO[one_letter]
    return one_letter, three_letter, name


@dataclass
class FastaRecord:
    header: str
    sequence: str


def _clean(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("U", "T")


def read_fasta(path: str | Path) -> list[FastaRecord]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, _clean("".join(chunks))))
                header, chunks = line[1:].strip(), []
            else:
                if header is None:
                    header = "(no header)"
                chunks.append(line)

    if header is not None:
        records.append(FastaRecord(header, _clean("".join(chunks))))

    records = [r for r in records if r.sequence]
    if not records:
        raise ValueError("No sequence was found. Is this a FASTA file?")
    return records


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass
class Orf:
    record: int
    strand: str
    start: int
    end: int

    @property
    def length_codons(self) -> int:
        return (self.end - self.start) // 3


def read_until_stop(sequence: str, start: int) -> int | None:
    for j in range(start, len(sequence) - 2, 3):
        codon = sequence[j:j + 3]
        if codon in STOP_CODONS:
            return j + 3
        if codon not in VALID_DNA_CODONS:
            return None
    return None


def find_orfs(sequence: str, min_codons: int = MIN_ORF_CODONS,
              strand: str = "+", record: int = 0) -> list[Orf]:
    orfs: list[Orf] = []
    i = 0
    while i <= len(sequence) - 3:
        if sequence[i:i + 3] == START_CODON:
            end = read_until_stop(sequence, i)
            if end is not None and (end - i) // 3 >= min_codons:
                orfs.append(Orf(record, strand, i, end))
                i = end
                continue
        i += 1
    return orfs


def count_orf_codons(sequence: str, orfs: list[Orf]) -> Counter:
    counts = Counter()
    for orf in orfs:
        counts.update(sequence[j:j + 3] for j in range(orf.start, orf.end, 3))
    return counts


def count_windows(sequence: str) -> tuple[Counter, int]:
    windows = Counter(sequence[i:i + 3] for i in range(len(sequence) - 2))
    valid = Counter({codon: n for codon, n in windows.items()
                     if codon in VALID_DNA_CODONS})
    skipped = sum(windows.values()) - sum(valid.values())
    return valid, skipped


@dataclass
class CodonStat:
    rank: int
    dna: str
    rna: str
    count: int
    frequency: float
    aa_one: str
    aa_three: str
    aa_name: str


@dataclass
class AnalysisResult:
    path: Path
    method: str
    records: list[tuple[str, int]]
    genome_length: int
    codons_counted: int
    counts: Counter
    top: list[CodonStat]
    min_codons: int = 0
    both_strands: bool = False
    orfs: list[Orf] = field(default_factory=list)
    longest_orf: Orf | None = None
    windows_skipped: int = 0


def analyze(path: str | Path, method: str = "orf",
            min_codons: int = MIN_ORF_CODONS, both_strands: bool = False,
            top_n: int = TOP_N) -> AnalysisResult:
    if method not in ("orf", "window"):
        raise ValueError(f"Unknown method: {method!r}")
    if min_codons < 1:
        raise ValueError("The minimum ORF length must be at least 1 codon.")

    records = read_fasta(path)
    genome_length = sum(len(r.sequence) for r in records)

    acgt = sum(r.sequence.count(base) for r in records for base in "ACGT")
    if acgt / genome_length < 0.9:
        raise ValueError(
            f"Only {100 * acgt / genome_length:.0f}% of the letters are A/C/G/T. "
            "This does not look like a DNA/RNA genome (is it a protein FASTA?)."
        )

    total = Counter()
    orfs: list[Orf] = []
    strand_sequences: dict[tuple[int, str], str] = {}
    skipped = 0
    for index, record in enumerate(records):
        if method == "window":
            counts, bad = count_windows(record.sequence)
            total.update(counts)
            skipped += bad
            continue

        strands = {"+": record.sequence}
        if both_strands:
            strands["-"] = reverse_complement(record.sequence)
        for strand, sequence in strands.items():
            strand_sequences[(index, strand)] = sequence
            orfs.extend(find_orfs(sequence, min_codons, strand, index))

    longest_orf = None
    if method == "orf":
        if not orfs:
            hint = "" if both_strands else ", or also scan the reverse strand"
            raise ValueError(f"No ORF of at least {min_codons} codons was found. "
                             f"Try a smaller minimum length{hint}.")
        longest_orf = max(orfs, key=lambda orf: orf.length_codons)
        sequence = strand_sequences[(longest_orf.record, longest_orf.strand)]
        total = count_orf_codons(sequence, [longest_orf])

    codons_counted = sum(total.values())
    if codons_counted == 0:
        raise ValueError("The sequence is too short to contain a single codon.")

    ranked = sorted(total.items(), key=lambda item: (-item[1], item[0]))[:top_n]

    top = []
    for rank, (dna, count) in enumerate(ranked, start=1):
        rna = transcribe(dna)
        one, three, name = translate(rna)
        top.append(CodonStat(rank, dna, rna, count,
                             100 * count / codons_counted, one, three, name))

    return AnalysisResult(
        path=Path(path), method=method,
        records=[(r.header, len(r.sequence)) for r in records],
        genome_length=genome_length, codons_counted=codons_counted,
        counts=total, top=top, min_codons=min_codons,
        both_strands=both_strands, orfs=orfs, longest_orf=longest_orf,
        windows_skipped=skipped,
    )


def describe_method(result: AnalysisResult) -> str:
    if result.method == "window":
        return "sliding window of 3, step 1"
    strands = "both strands" if result.both_strands else "+ strand only"
    return f"ORF finder, longest ORF (ORFs of at least {result.min_codons} codons, {strands})"


def describe_longest_orf(result: AnalysisResult) -> str:
    orf = result.longest_orf
    if orf is None:
        return ""
    if orf.strand == "+":
        first, last = orf.start + 1, orf.end
        where = "+ strand"
    else:
        length = result.records[orf.record][1]
        first, last = length - orf.end + 1, length - orf.start
        where = "- strand"
    record = f", record {orf.record + 1}" if len(result.records) > 1 else ""
    return f"{first}..{last} ({where}{record}), {orf.length_codons:,} codons"


def describe_counts(result: AnalysisResult) -> str:
    if result.method == "window":
        return (f"Windows counted: {result.codons_counted:,}  "
                f"(skipped {result.windows_skipped:,} with non-ACGT letters)")
    return (f"ORFs found: {len(result.orfs):,}  |  "
            f"Longest ORF: {describe_longest_orf(result)}")


def describe_amino_acid(stat: CodonStat) -> str:
    if stat.aa_one == "*":
        return "Stop codon"
    return f"{stat.aa_three} ({stat.aa_one}) - {stat.aa_name}"


def describe_records(result: AnalysisResult) -> str:
    header = result.records[0][0]
    extra = len(result.records) - 1
    if extra:
        header += f"  (+{extra} more record{'s' if extra > 1 else ''})"
    return header


def format_report(result: AnalysisResult) -> str:
    lines = [
        f"File:          {result.path.name}",
        f"Sequence:      {describe_records(result)}",
        f"Genome length: {result.genome_length:,} nt",
        f"Method:        {describe_method(result)}",
        describe_counts(result),
        "",
        f"Top {len(result.top)} codons "
        + ("inside the longest ORF:" if result.method == "orf" else "(sliding window):"),
        f"{'#':>3}  {'DNA':<5} {'RNA':<5} {'Amino acid':<26} {'Count':>10} {'Freq':>8}",
        "-" * 63,
    ]
    for s in result.top:
        lines.append(
            f"{s.rank:>3}  {s.dna:<5} {s.rna:<5} {describe_amino_acid(s):<26} "
            f"{s.count:>10,} {s.frequency:>7.2f}%"
        )
    return "\n".join(lines)


class CodonApp:
    COLUMNS = (
        ("rank", "#", 40),
        ("dna", "Codon (DNA)", 100),
        ("rna", "Codon (RNA)", 100),
        ("aa", "Amino acid", 240),
        ("count", "Count", 110),
        ("freq", "Frequency", 110),
    )

    def __init__(self, root: tk.Tk):
        self.root = root
        self.current_path: str | None = None
        root.title("Codon Frequency Analyzer")
        root.geometry("980x800")
        root.minsize(760, 600)

        bar = ttk.Frame(root, padding=(10, 10, 10, 4))
        bar.pack(fill="x")
        ttk.Button(bar, text="Open FASTA file...", command=self.choose_file).pack(side="left")
        self.file_var = tk.StringVar(value="No file loaded")
        ttk.Label(bar, textvariable=self.file_var).pack(side="left", padx=10)

        opts = ttk.Frame(root, padding=(10, 4, 10, 8))
        opts.pack(fill="x")
        self.method_var = tk.StringVar(value="orf")
        ttk.Label(opts, text="Method:").pack(side="left")
        ttk.Radiobutton(opts, text="ORF finder", value="orf", variable=self.method_var,
                        command=self.on_settings_changed).pack(side="left", padx=(6, 0))
        ttk.Radiobutton(opts, text="Sliding window (all triplets)", value="window",
                        variable=self.method_var,
                        command=self.on_settings_changed).pack(side="left", padx=(10, 20))

        ttk.Label(opts, text="Min ORF length (codons):").pack(side="left")
        self.min_var = tk.StringVar(value=str(MIN_ORF_CODONS))
        self.min_spin = ttk.Spinbox(opts, from_=1, to=100000, increment=10, width=7,
                                    textvariable=self.min_var,
                                    command=self.on_settings_changed)
        self.min_spin.pack(side="left", padx=(6, 20))
        self.min_spin.bind("<Return>", lambda _event: self.on_settings_changed())

        self.strand_var = tk.BooleanVar(value=False)
        self.strand_check = ttk.Checkbutton(opts, text="Also scan the reverse strand",
                                            variable=self.strand_var,
                                            command=self.on_settings_changed)
        self.strand_check.pack(side="left")

        self.info_var = tk.StringVar(
            value="Choose a virus or bacteria genome in FASTA format "
                  "(.fasta, .fa, .fna, or the same gzipped as .gz)."
        )
        ttk.Label(root, textvariable=self.info_var, padding=(10, 0),
                  wraplength=940, justify="left").pack(fill="x")

        self.title_var = tk.StringVar(value=f"Top {TOP_N} most frequent codons")
        ttk.Label(root, textvariable=self.title_var, font=("TkDefaultFont", 11, "bold"),
                  padding=(10, 12, 10, 4)).pack(anchor="w")
        self.table = ttk.Treeview(root, columns=[c[0] for c in self.COLUMNS],
                                  show="headings", height=TOP_N)
        for key, title, width in self.COLUMNS:
            self.table.heading(key, text=title)
            self.table.column(key, width=width, anchor="center")
        self.table.column("aa", anchor="w")
        self.table.tag_configure("stop", foreground="#b00020")
        self.table.pack(fill="x", padx=10)

        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            from matplotlib.figure import Figure
        except ImportError:
            self.has_chart = False
            ttk.Label(root, padding=10,
                      text="Install matplotlib (pip install matplotlib) "
                           "to also see a bar chart.").pack()
        else:
            self.has_chart = True
            self.figure = Figure(figsize=(9.4, 3.4), dpi=100)
            self.ax = self.figure.add_subplot(111)
            self.canvas = FigureCanvasTkAgg(self.figure, master=root)
            self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

    def on_settings_changed(self):
        orf_mode = self.method_var.get() == "orf"
        state = "normal" if orf_mode else "disabled"
        self.min_spin.config(state=state)
        self.strand_check.config(state=state)
        if self.current_path:
            self.load_file(self.current_path)

    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select a genome FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.fas *.gz"),
                       ("All files", "*.*")],
        )
        if path:
            self.load_file(path)

    def read_min_codons(self) -> int | None:
        try:
            value = int(self.min_var.get())
            if value < 1:
                raise ValueError
            return value
        except ValueError:
            messagebox.showerror("Invalid setting",
                                 "The minimum ORF length must be a whole number of at least 1.")
            self.min_var.set(str(MIN_ORF_CODONS))
            return None

    def load_file(self, path: str | Path):
        min_codons = self.read_min_codons()
        if min_codons is None:
            return
        self.root.config(cursor="watch")
        self.root.update_idletasks()
        try:
            result = analyze(path, method=self.method_var.get(),
                             min_codons=min_codons, both_strands=self.strand_var.get())
        except (OSError, ValueError) as exc:
            messagebox.showerror("Could not analyse the file", str(exc))
            return
        finally:
            self.root.config(cursor="")
        self.current_path = str(path)
        self.show(result)

    def show(self, result: AnalysisResult):
        self.file_var.set(str(result.path.name))
        self.info_var.set(
            f"{describe_records(result)}\n"
            f"Genome length: {result.genome_length:,} nt    |    {describe_counts(result)}"
        )
        where = "inside the longest ORF" if result.method == "orf" else "(sliding window)"
        self.title_var.set(f"Top {len(result.top)} most frequent codons {where}")

        self.table.delete(*self.table.get_children())
        for s in result.top:
            self.table.insert(
                "", "end",
                values=(s.rank, s.dna, s.rna, describe_amino_acid(s),
                        f"{s.count:,}", f"{s.frequency:.2f} %"),
                tags=("stop",) if s.aa_one == "*" else (),
            )

        if self.has_chart:
            self.draw_chart(result)

    def draw_chart(self, result: AnalysisResult):
        ax = self.ax
        ax.clear()
        labels = [f"{s.dna}\n{s.aa_three}" for s in result.top]
        values = [s.frequency for s in result.top]
        colors = ["#c0392b" if s.aa_one == "*" else "#2e86de" for s in result.top]
        bars = ax.bar(labels, values, color=colors)
        ax.bar_label(bars, fmt="%.2f", fontsize=8, padding=2)
        ax.set_ylabel("Frequency (%)")
        ax.set_ylim(0, max(values) * 1.15)
        if result.method == "orf":
            title = f"{result.path.name} - longest ORF: {describe_longest_orf(result)}"
        else:
            title = f"{result.path.name} - {describe_method(result)}"
        ax.set_title(title, fontsize=10)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        self.figure.tight_layout()
        self.canvas.draw()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Codon frequencies of a virus or bacteria genome (FASTA).")
    parser.add_argument("fasta", nargs="?",
                        help="FASTA file to analyse; leave it out to open the window")
    parser.add_argument("--method", choices=["orf", "window"], default="orf",
                        help="orf = codons inside ORFs (default); "
                             "window = every triplet, sliding window of 3, step 1")
    parser.add_argument("--min-codons", type=int, default=MIN_ORF_CODONS,
                        help=f"minimum ORF length in codons (default {MIN_ORF_CODONS})")
    parser.add_argument("--both-strands", action="store_true",
                        help="also look for ORFs on the reverse strand")
    args = parser.parse_args(argv)

    if args.fasta:
        try:
            result = analyze(args.fasta, method=args.method,
                             min_codons=args.min_codons, both_strands=args.both_strands)
        except (OSError, ValueError) as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        print(format_report(result))
        return 0

    root = tk.Tk()
    CodonApp(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())