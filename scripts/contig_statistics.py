#!/usr/bin/env python3
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import GC

def main(fasta_txt, output_path):
    fasta_txt = Path(fasta_txt)
    output_path = Path(output_path)

    # Read FASTA paths from the .txt file
    with open(fasta_txt) as f:
        fasta_files = [line.strip() for line in f if line.strip()]

    with open(output_path, 'w') as out:
        out.write("contig_id\tlength\tgc_percent\n")
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_len = len(record.seq)
                gc_pct = round(GC(record.seq), 2)
                out.write(f"{record.id}\t{seq_len}\t{gc_pct}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <fasta_paths.txt> <output_tsv>")
        sys.exit(1)
    fasta_txt, out_file = sys.argv[1:]
    main(fasta_txt, out_file)
