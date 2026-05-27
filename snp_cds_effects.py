from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict

# Load genome
genome = SeqIO.to_dict(SeqIO.parse("genome.fa", "fasta"))

# Load CDS annotations
cds_by_transcript = defaultdict(list)
with open("AnoSag2.1_cds.gtf") as gtf:
    for line in gtf:
        if line.startswith("#") or '\tCDS\t' not in line:
            continue
        chrom, _, _, start, end, _, strand, _, attrs = line.strip().split('\t')
        start, end = int(start), int(end)
        transcript_id = [x.split('"')[1] for x in attrs.split(';') if "transcript_id" in x][0]
        cds_by_transcript[transcript_id].append((chrom, start, end, strand))

# Load SNPs
snps = pd.read_csv("snp_genome_coords.tsv", sep="\t", header=None,
                   names=["chrom", "start", "end", "radtag", "rad_pos", "snp_id", "ref", "alt", "strand"])

# Step: Classify SNPs
results = []

for i, row in snps.iterrows():
    chrom, pos, ref, alt = row["chrom"], row["start"], row["ref"], row["alt"]

    for tid, cds_parts in cds_by_transcript.items():
        parts = [p for p in cds_parts if p[0] == chrom]
        if not parts:
            continue
        strand = parts[0][3]

        # Flatten and sort CDS regions
        cds_coords = []
        for _, s, e, _ in sorted(parts, key=lambda x: x[1]):
            cds_coords.extend(range(s, e + 1))
        if strand == "-":
            cds_coords = cds_coords[::-1]  # Reverse for negative strand

        if pos not in cds_coords:
            continue

        cds_index = cds_coords.index(pos)
        codon_index = cds_index // 3
        codon_start = codon_index * 3

        if codon_start + 2 >= len(cds_coords):
            continue  # incomplete codon

        codon_positions = cds_coords[codon_start:codon_start + 3]
        codon_seq = "".join([genome[chrom].seq[p - 1] for p in codon_positions])

        if strand == "-":
            codon_seq = str(Seq(codon_seq).reverse_complement())

        # Mutate the codon
        snp_in_codon_pos = cds_index % 3
        if strand == "-":
            snp_in_codon_pos = 2 - snp_in_codon_pos
            alt = str(Seq(alt).reverse_complement())

        codon_list = list(codon_seq)
        codon_list[snp_in_codon_pos] = alt.upper()
        mutated_codon = "".join(codon_list)

        # Translate and compare
        original_aa = Seq(codon_seq).translate()
        mutated_aa = Seq(mutated_codon).translate()

        classification = "synonymous" if original_aa == mutated_aa else "non-synonymous"
        results.append([
            row["snp_id"], chrom, pos, ref, alt,
            codon_seq, mutated_codon, original_aa, mutated_aa, classification
        ])

# Output
out_df = pd.DataFrame(results, columns=[
    "snp_id", "chrom", "pos", "ref", "alt",
    "original_codon", "mutated_codon", "original_aa", "mutated_aa", "classification"
])
out_df.to_csv("snp_codon_effects.tsv", sep="\t", index=False)
