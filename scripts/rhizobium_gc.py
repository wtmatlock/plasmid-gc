from Bio import SeqIO

# four-fold degenerate codons (for GC4)
FOUR_FOLD = {
    'GCT','GCC','GCA','GCG',  # Ala
    'GGT','GGC','GGA','GGG',  # Gly
    'CCT','CCC','CCA','CCG',  # Pro
    'ACT','ACC','ACA','ACG',  # Thr
    'GTT','GTC','GTA','GTG',  # Val
    'TCT','TCC','TCA','TCG',  # Ser (part)
    'CGT','CGC','CGA','CGG',  # Arg (part)
    'CTT','CTC','CTA','CTG'   # Leu (part)
}

input_file = "ncbi_dataset/data/GCF_000009265.1/genomic.gbff"
output_file = "rhizobium_gc.tsv"

records = []

for record in SeqIO.parse(input_file, "genbank"):
    # get plasmid name from source feature if present
    plasmid_name = None
    for feature in record.features:
        if feature.type == "source":
            plasmid_name = feature.qualifiers.get("plasmid", [None])[0]
            break
    replicon = plasmid_name if plasmid_name else record.name

    for feature in record.features:
        if feature.type == "CDS":
            locus = feature.qualifiers.get("old_locus_tag", [""])[0]
            if not locus:
                locus = feature.qualifiers.get("locus_tag", [""])[0]

            seq = feature.extract(record.seq)
            seq_str = str(seq).upper()
            gc = seq_str.count("G") + seq_str.count("C")
            total = len(seq_str)
            gc_content = (gc / total * 100) if total > 0 else 0

            # GC1, GC2, GC3, GC4
            gc1 = gc2 = gc3 = gc4 = 'NA'
            n_codons = fourfold_count = 'NA'
            if total % 3 == 0 and total > 0:
                n_codons = total // 3
                counts = [0, 0, 0]
                fourfold_count = 0
                gc4_count = 0
                for i in range(n_codons):
                    codon = seq_str[3*i:3*i+3]
                    for pos in range(3):
                        if codon[pos] in ('G','C'):
                            counts[pos] += 1
                    if codon in FOUR_FOLD:
                        fourfold_count += 1
                        if codon[2] in ('G','C'):
                            gc4_count += 1
                if n_codons:
                    gc1 = (counts[0] / n_codons) * 100
                    gc2 = (counts[1] / n_codons) * 100
                    gc3 = (counts[2] / n_codons) * 100
                if fourfold_count > 0:
                    gc4 = (gc4_count / fourfold_count) * 100

            records.append({
                "old_locus_tag": locus,
                "replicon": replicon,
                "length": total,
                "gc_content": f"{gc_content:.2f}",
                "gc1": f"{gc1:.2f}" if gc1 != 'NA' else 'NA',
                "gc2": f"{gc2:.2f}" if gc2 != 'NA' else 'NA',
                "gc3": f"{gc3:.2f}" if gc3 != 'NA' else 'NA',
                "gc4": f"{gc4:.2f}" if gc4 != 'NA' else 'NA'
            })

with open(output_file, "w") as out:
    out.write("old_locus_tag\treplicon\tlength\tgc_content\tgc1\tgc2\tgc3\tgc4\n")
    for rec in records:
        out.write(f"{rec['old_locus_tag']}\t{rec['replicon']}\t{rec['length']}\t{rec['gc_content']}\t{rec['gc1']}\t{rec['gc2']}\t{rec['gc3']}\t{rec['gc4']}\n")

print(f"Results written to {output_file}")