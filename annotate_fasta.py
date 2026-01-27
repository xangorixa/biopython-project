import sys

import os

from collections import Counter

# allowed character sets for basic type detection
DNA_Chars = set("ACGTN")
RNA_CHARS = set("ACGUN")
Protein_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZ*")

def parse_fasta(path):
    """ returns a list of (header, sequence) tuples.header: text after '>' sequence: concatenated sequence lines, uppercased,whitespace removed
    """
    records = []
    header = None
    seq_parts = []

    with open(path, "r", encoding="utf-8") as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    records.append((header, seq))

                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_num}") 
                seq_parts = []

            else:
                if header is None:
                    raise ValueError(f"Found sequence befroe first header at line {line_num}")
                seq_parts.append("".join(line.split()))

    if header is not None:
        seq = "".join(seq_parts).upper()
        records.append((header, seq))

    if not records:
        raise ValueError("No FASTA records found. File must contain a line starting with'>'")

    return records

def guess_sequence_type(seq):
    letters = set(seq)

    if letters.issubset(DNA_Chars):
        return "DNA"

    if letters.issubset(RNA_CHARS):
        return "RNA"

    if letters.issubset(Protein_CHARS):
        return "protein"

    return "unknown"

def counts(seq):
    return Counter(seq)

def gc_percent(seq):
    g = seq.count("G")
    c = seq.count("C")
    total = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
    return None if total == 0 else 100 * (g + c) / total

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 annotate_fasta.py <fasta_file>")
        return

    path = sys.argv[1]
    records = parse_fasta(path)

    # results file name based on input file name
    base = os.path.basename(path)
    out_path = f"{base}.results.txt"   # e.g., QSOX1.txt.results.txt

    with open(out_path, "w", encoding="utf-8") as r:
        def out(*args):
            line = " ".join(str(a) for a in args)
            print(line)         # still show on screen
            r.write(line + "\n")  # also save to file

        out("Records:", len(records))
        for header, seq in records:
            out("Header:", header)
            out("Length:", len(seq))

            seq_type = guess_sequence_type(seq)
            out("Type:", seq_type)

            c = counts(seq)
            out("Counts (top):", dict(c.most_common(10)))

            if seq_type == "DNA":
                gc = gc_percent(seq)
                out("GC%:", round(gc, 2) if gc is not None else "n/a")

    print(f"Wrote results to: {out_path}")

if __name__ == "__main__":
    main()
