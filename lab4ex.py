import math
import random
import tkinter as tk
import tkinter.font as tkfont
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog

SEQUENCE_MIN, SEQUENCE_MAX = 1000, 3000
DNA_LETTERS = set("ACGTRYSWKMBDHVN")

SAMPLE_COUNT = 10
SAMPLE_MIN, SAMPLE_MAX = 100, 3000
LADDER = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000]
GEL_LENGTH_MM = 80
RUN_TIME_MIN = 45
MOBILITY = 0.6
SIZE_THAT_DOES_NOT_MOVE = 50_000

PX_PER_MM = 6
LANE_WIDTH = 40
BAND_WIDTH = 28
GEL_LEFT = 70
WELL_Y = 52
FRAMES = 60
FRAME_MS = 40

UNIT_MIN, UNIT_MAX = 3, 5
REPEAT_MIN, REPEAT_MAX = 6, 10
FLANK = 6
LINE_WIDTH = 60

INSTRUCTIONS = f"""How to get a DNA sequence from NCBI:
  1. Search https://www.ncbi.nlm.nih.gov/nuccore
     Tip: add  AND {SEQUENCE_MIN}:{SEQUENCE_MAX}[SLEN]  to get sequences of {SEQUENCE_MIN}-{SEQUENCE_MAX} nt.
  2. Open a record, then Send to > Complete Record > File
     > Format: FASTA > Create File.
  3. Click "Open FASTA file..." and choose the downloaded file.
"""

Repeat = tuple[int, str, int]
Lane = tuple[str, list[int]]


@dataclass
class Band:
    tag: str
    length: int
    y: float = 0.0


def parse_fasta(text: str) -> tuple[str, str]:
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("this is not a FASTA file (the first line must start with '>')")
    sequence = []
    for line in lines[1:]:
        if line.startswith(">"):
            break
        sequence.append(line.strip())
    return lines[0][1:].strip(), "".join(sequence).upper()


def read_fasta_file(path: str) -> tuple[str, str]:
    try:
        text = Path(path).read_text(encoding="utf-8")
    except UnicodeDecodeError:
        raise ValueError("this is not a text (FASTA) file") from None
    header, sequence = parse_fasta(text)
    wrong = sorted(set(sequence) - DNA_LETTERS)
    if wrong:
        raise ValueError(f"the sequence has letters that are not DNA bases: {' '.join(wrong[:10])}")
    if not SEQUENCE_MIN <= len(sequence) <= SEQUENCE_MAX:
        raise ValueError(f"the sequence has {len(sequence)} nucleotides, the exercises need {SEQUENCE_MIN}-{SEQUENCE_MAX}")
    return header, sequence


def take_samples(sequence: str, count: int) -> list[tuple[int, str]]:
    longest = min(SAMPLE_MAX, len(sequence))
    samples = []
    for _ in range(count):
        length = random.randint(SAMPLE_MIN, longest)
        start = random.randint(0, len(sequence) - length)
        samples.append((start, sequence[start:start + length]))
    return samples


def migration_distance(length: int, minutes: float) -> float:
    speed = MOBILITY * math.log10(SIZE_THAT_DOES_NOT_MOVE / length)
    return min(max(speed, 0.0) * minutes, GEL_LENGTH_MM)


def band_color(length: int) -> str:
    level = 120 + 135 * min(length, 3000) // 3000
    return f"#{level:02x}{level:02x}{level:02x}"


def samples_report(samples: list[tuple[int, str]]) -> str:
    lines = [f"{len(samples)} random samples (positions are 1-based), distance after {RUN_TIME_MIN} min:", ""]
    for number, (start, fragment) in enumerate(samples, 1):
        distance = migration_distance(len(fragment), RUN_TIME_MIN)
        lines.append(f"  S{number:<3}{start + 1:>6}-{start + len(fragment):<6}{len(fragment):>6} bp"
                     f"{distance:7.1f} mm   {fragment[:20]}...")
    order = sorted(range(len(samples)), key=lambda i: len(samples[i][1]))
    lines.append("")
    lines.append("Fastest to slowest: " + " ".join(f"S{i + 1}" for i in order))
    return "\n".join(lines)


def is_primitive(unit: str) -> bool:
    size = len(unit)
    return all(unit != unit[:part] * (size // part) for part in range(1, size) if size % part == 0)


def find_repeats(sequence: str) -> tuple[list[Repeat], int]:
    repeats: list[Repeat] = []
    too_long = 0
    for size in range(UNIT_MIN, UNIT_MAX + 1):
        for start in range(len(sequence) - 2 * size + 1):
            if start > 0 and sequence[start - 1] == sequence[start - 1 + size]:
                continue
            unit = sequence[start:start + size]
            copies = 1
            while sequence[start + copies * size:start + (copies + 1) * size] == unit:
                copies += 1
            if copies < 2 or not is_primitive(unit) or set(unit) - set("ACGT"):
                continue
            if copies * size > REPEAT_MAX:
                too_long += 1
            elif copies * size >= REPEAT_MIN:
                repeats.append((start, unit, copies))
    repeats.sort()
    return repeats, too_long


def repeats_report(sequence: str, repeats: list[Repeat], too_long: int) -> str:
    lines = [f"{'#':>3}  {'position':<11}{'unit':<7}{'copies':>6}{'length':>8}   repeat with {FLANK} bases on each side"]
    for number, (start, unit, copies) in enumerate(repeats, 1):
        end = start + len(unit) * copies
        left = sequence[max(start - FLANK, 0):start].lower()
        right = sequence[end:end + FLANK].lower()
        position = f"{start + 1}-{end}"
        lines.append(f"{number:>3}  {position:<11}{unit:<7}{copies:>6}{end - start:>8}   "
                     f"{left}[{sequence[start:end]}]{right}")
    sizes = Counter(len(unit) for _, unit, _ in repeats)
    units = Counter(unit for _, unit, _ in repeats)
    lines.append(f"\nFound {len(repeats)} repeats of {REPEAT_MIN}-{REPEAT_MAX} bases:")
    lines += [f"  unit of {size} bases: {sizes[size]}" for size in range(UNIT_MIN, UNIT_MAX + 1)]
    if units:
        lines.append("Most frequent units: " + ", ".join(f"{unit} ({count}x)" for unit, count in units.most_common(5)))
    if too_long:
        lines.append(f"Skipped {too_long} tandem repeat(s) longer than {REPEAT_MAX} bases.")
    return "\n".join(lines)


def sequence_map(sequence: str, repeats: list[Repeat]) -> str:
    marked = [False] * len(sequence)
    for start, unit, copies in repeats:
        for position in range(start, start + len(unit) * copies):
            marked[position] = True
    text = "".join(base if mark else base.lower() for base, mark in zip(sequence, marked))
    lines = ["Sequence map (repeats in UPPERCASE):"]
    for line_start in range(0, len(text), LINE_WIDTH):
        line = text[line_start:line_start + LINE_WIDTH]
        lines.append(f"{line_start + 1:>6} " + " ".join(line[i:i + 10] for i in range(0, len(line), 10)))
    return "\n".join(lines)


class App:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.header = ""
        self.sequence = ""
        self.bands: list[Band] = []
        self.clock = 0
        self.job: str | None = None
        root.title("Bioinformatics lab")

        top = tk.Frame(root)
        top.pack(fill="x", padx=8, pady=(8, 4))
        tk.Button(top, text="Open FASTA file...", command=self.open_file).pack(side="left")
        self.info = tk.Label(top, text="No file loaded", anchor="w")
        self.info.pack(side="left", padx=10)

        buttons = tk.Frame(root)
        buttons.pack(fill="x", padx=8, pady=4)
        tk.Button(buttons, text="Exercise 1 - Gel electrophoresis", command=self.exercise1).pack(side="left")
        tk.Button(buttons, text="Exercise 2 - Find repeats", command=self.exercise2).pack(side="left", padx=6)

        main = tk.Frame(root)
        main.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        self.canvas = tk.Canvas(main, width=GEL_LEFT + (SAMPLE_COUNT + 2) * LANE_WIDTH + 10,
                                height=WELL_Y + GEL_LENGTH_MM * PX_PER_MM + 20, highlightthickness=0)
        self.canvas.pack(side="right", anchor="n", padx=(8, 0))
        self.canvas.create_text(GEL_LEFT, WELL_Y, text="The gel appears here (Exercise 1)", anchor="w")
        scrollbar = tk.Scrollbar(main)
        scrollbar.pack(side="right", fill="y")
        family = "Consolas" if "Consolas" in tkfont.families() else tkfont.nametofont("TkFixedFont").actual("family")
        self.output = tk.Text(main, width=76, height=34, font=(family, 10), wrap="none", state="disabled",
                              yscrollcommand=scrollbar.set)
        self.output.pack(side="left", fill="both", expand=True)
        scrollbar.configure(command=self.output.yview)
        self.show(INSTRUCTIONS)

    def show(self, text: str) -> None:
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("end", text)
        self.output.configure(state="disabled")

    def open_file(self) -> None:
        path = filedialog.askopenfilename(title="Open FASTA file",
                                          filetypes=[("FASTA files", "*.fasta *.fa *.fna *.fas *.fsa *.txt"),
                                                     ("All files", "*.*")])
        if path:
            self.load_file(path)

    def load_file(self, path: str) -> None:
        name = Path(path).name
        try:
            self.header, self.sequence = read_fasta_file(path)
        except (OSError, ValueError) as error:
            self.header, self.sequence = "", ""
            self.info.configure(text=f"Could not load {name}")
            self.show(f"Could not load {name}:\n{error}\n\n{INSTRUCTIONS}")
            return
        self.info.configure(text=f"{name}:  {self.header[:60]}  ({len(self.sequence)} nt)")
        self.show(f"Loaded {name}\n{self.header}\nLength: {len(self.sequence)} nt\n\nChoose an exercise.")

    def exercise1(self) -> None:
        if not self.sequence:
            self.show("Open a FASTA file first.\n\n" + INSTRUCTIONS)
            return
        samples = take_samples(self.sequence, SAMPLE_COUNT)
        lengths = [len(fragment) for _, fragment in samples]
        self.show("EXERCISE 1 - Gel electrophoresis\n\n" + samples_report(samples))
        lanes = [("Ladder", LADDER)] + [(f"S{n}", [length]) for n, length in enumerate(lengths, 1)] + [("Mix", lengths)]
        if self.job is not None:
            self.root.after_cancel(self.job)
        self.draw_gel(lanes)
        self.animate(0)

    def draw_gel(self, lanes: list[Lane]) -> None:
        c = self.canvas
        c.delete("all")
        right = GEL_LEFT + LANE_WIDTH * len(lanes)
        bottom = WELL_Y + GEL_LENGTH_MM * PX_PER_MM
        c.create_rectangle(GEL_LEFT, WELL_Y - 12, right, bottom + 10, fill="black")
        self.clock = c.create_text(GEL_LEFT, 10, anchor="w")
        c.create_text(4, WELL_Y - 6, text="(-)", anchor="w")
        c.create_text(4, bottom + 4, text="(+)", anchor="w")
        self.bands = []
        for index, (name, lengths) in enumerate(lanes):
            x = GEL_LEFT + LANE_WIDTH * index + LANE_WIDTH / 2
            c.create_text(x, WELL_Y - 22, text=name)
            c.create_rectangle(x - BAND_WIDTH / 2, WELL_Y - 8, x + BAND_WIDTH / 2, WELL_Y - 3, outline="gray50")
            for length in lengths:
                tag = f"band{len(self.bands)}"
                c.create_rectangle(x - BAND_WIDTH / 2, -2, x + BAND_WIDTH / 2, 2, fill=band_color(length),
                                   outline="", tags=tag)
                if name == "Ladder":
                    c.create_text(GEL_LEFT - 6, 0, text=f"{length} bp", anchor="e", tags=tag)
                elif name != "Mix":
                    c.create_text(x, 4, text=str(length), anchor="n", fill="gray70", tags=tag)
                self.bands.append(Band(tag, length))

    def animate(self, frame: int) -> None:
        minutes = RUN_TIME_MIN * frame / FRAMES
        for band in self.bands:
            y = WELL_Y + migration_distance(band.length, minutes) * PX_PER_MM
            self.canvas.move(band.tag, 0, y - band.y)
            band.y = y
        self.canvas.itemconfigure(self.clock, text=f"Migration: t = {minutes:.1f} / {RUN_TIME_MIN} min")
        self.job = self.root.after(FRAME_MS, self.animate, frame + 1) if frame < FRAMES else None

    def exercise2(self) -> None:
        if not self.sequence:
            self.show("Open a FASTA file first.\n\n" + INSTRUCTIONS)
            return
        repeats, too_long = find_repeats(self.sequence)
        self.show(f"EXERCISE 2 - Repeats of {REPEAT_MIN}-{REPEAT_MAX} bases made of a {UNIT_MIN}-{UNIT_MAX} base unit\n\n"
                  + repeats_report(self.sequence, repeats, too_long) + "\n\n" + sequence_map(self.sequence, repeats))


def main() -> None:
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
