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

# core region: drop start codon + first 11 codons (Lewin et al. 2023 5' window)
# and the last 10 codons + stop codon
CORE_5P = 12
CORE_3P = 11

input_file = "ncbi_dataset/data/GCF_000009265.1/genomic.gbff"
output_file = "rhizobium_gc.tsv"


def codon_gc(codons):
    """Return gc1, gc2, gc3, gc4 (percent, or 'NA') for a list of codons."""
    if not codons:
        return 'NA', 'NA', 'NA', 'NA'
    counts = [0, 0, 0]
    fourfold_count = 0
    gc4_count = 0
    for codon in codons:
        for pos in range(3):
            if codon[pos] in ('G','C'):
                counts[pos] += 1
        if codon in FOUR_FOLD:
            fourfold_count += 1
            if codon[2] in ('G','C'):
                gc4_count += 1
    n = len(codons)
    gc4 = (gc4_count / fourfold_count) * 100 if fourfold_count > 0 else 'NA'
    return counts[0] / n * 100, counts[1] / n * 100, counts[2] / n * 100, gc4


def fmt(x):
    return f"{x:.2f}" if x != 'NA' else 'NA'


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

            # GC1, GC2, GC3, GC4 (whole gene, and core region)
            gc1 = gc2 = gc3 = gc4 = 'NA'
            gc1_core = gc2_core = gc3_core = gc4_core = 'NA'
            core_codons = 0
            if total % 3 == 0 and total > 0:
                codons = [seq_str[3*i:3*i+3] for i in range(total // 3)]
                gc1, gc2, gc3, gc4 = codon_gc(codons)
                core = codons[CORE_5P:len(codons) - CORE_3P]
                core_codons = len(core)
                gc1_core, gc2_core, gc3_core, gc4_core = codon_gc(core)

            records.append({
                "old_locus_tag": locus,
                "replicon": replicon,
                "length": total,
                "gc_content": f"{gc_content:.2f}",
                "gc1": fmt(gc1),
                "gc2": fmt(gc2),
                "gc3": fmt(gc3),
                "gc4": fmt(gc4),
                "core_codons": core_codons,
                "gc1_core": fmt(gc1_core),
                "gc2_core": fmt(gc2_core),
                "gc3_core": fmt(gc3_core),
                "gc4_core": fmt(gc4_core)
            })

columns = ["old_locus_tag", "replicon", "length", "gc_content",
           "gc1", "gc2", "gc3", "gc4",
           "core_codons", "gc1_core", "gc2_core", "gc3_core", "gc4_core"]

with open(output_file, "w") as out:
    out.write("\t".join(columns) + "\n")
    for rec in records:
        out.write("\t".join(str(rec[c]) for c in columns) + "\n")

print(f"Results written to {output_file}")
