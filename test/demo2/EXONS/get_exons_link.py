#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

# IUPAC DNA complement map
COMPLEMENT = str.maketrans(
    {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "a": "t", "t": "a", "c": "g", "g": "c",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N",
        "r": "y", "y": "r", "s": "s", "w": "w",
        "k": "m", "m": "k", "b": "v", "v": "b",
        "d": "h", "h": "d", "n": "n"
    }
)

def revcomp(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]

def read_first_fasta_seq(path: Path) -> str:
    """Read the first sequence from a FASTA file."""
    seq_lines = []
    with path.open() as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if seq_lines:
                    break
                else:
                    continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)

def load_exons(exon_dir: Path) -> Dict[int, str]:
    """
    Load exon sequences from files named exonN.fa (N is integer).
    Returns a dict {N: sequence}.
    """
    exons: Dict[int, str] = {}
    for p in sorted(exon_dir.iterdir()):
        m = re.fullmatch(r"exon(\d+)\.fa", p.name)
        if m and p.is_file():
            idx = int(m.group(1))
            seq = read_first_fasta_seq(p)
            if not seq:
                raise ValueError(f"No sequence found in {p}")
            exons[idx] = seq
    if not exons:
        raise ValueError(f"No exonN.fa files found in directory: {exon_dir}")
    return exons

def parse_scheme(s: str) -> List[Tuple[int, str]]:
    """
    Parse a scheme string into a list of (exon_index, strand) tuples.

    Supported formats:
      - With prefix:  exon_1+2-3+4+
      - Without:      1+2-3+4+

    Each exon is indicated by a number followed by + or -.
    """
    s = s.strip().replace(" ", "")
    if s.startswith("exon_"):
        s = s[len("exon_"):]
    tokens = re.findall(r"(\d+)([+-])", s)
    if not tokens:
        raise ValueError(f"Cannot parse scheme: {s}")
    return [(int(i), strand) for i, strand in tokens]

def build_transcript(scheme: str, exon_seqs: Dict[int, str]) -> str:
    """Build a transcript sequence according to the scheme."""
    pairs = parse_scheme(scheme)
    pieces = []
    for exon_id, strand in pairs:
        if exon_id not in exon_seqs:
            raise KeyError(f"Exon exon{exon_id}.fa referenced in scheme not found: {scheme}")
        seq = exon_seqs[exon_id]
        if strand == "-":
            seq = revcomp(seq)
        pieces.append(seq)
    return "".join(pieces)

def main():
    ap = argparse.ArgumentParser(
        description="Build transcript sequences by concatenating exon sequences with specified strand per exon."
    )
    ap.add_argument("--exon-dir", required=True, help="Directory containing exonN.fa files, e.g., ../")
    ap.add_argument("--scheme", action="append", required=True,
                    help="Transcript scheme, e.g., 'exon_1+2+3+4-5+6+9+' or '1+2-3+'. "
                         "Provide this option multiple times for multiple transcripts.")
    ap.add_argument("--out", required=True, help="Output FASTA file path")
    args = ap.parse_args()

    exon_dir = Path(args.exon_dir)
    if not exon_dir.is_dir():
        raise SystemExit(f"Directory does not exist: {exon_dir}")

    exon_seqs = load_exons(exon_dir)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as out_fh:
        for i, scheme in enumerate(args.scheme, 1):
            seq = build_transcript(scheme, exon_seqs)
            out_fh.write(f">tx{i} {scheme}\n")
            out_fh.write(seq + "\n")  # single-line sequence, no wrapping

    print(f"Done: wrote {len(args.scheme)} sequence(s) to {out_path}")

if __name__ == "__main__":
    main()