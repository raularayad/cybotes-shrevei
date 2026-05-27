import pandas as pd

def parse_attributes(attr_str):
    attrs = {}
    for attr in attr_str.strip().split(";"):
        if attr.strip():
            key, val = attr.strip().split(" ", 1)
            attrs[key] = val.strip('"')
    return attrs

# --- Load SNPs ---
snps_file = "snp_genome_coords.tsv"
snps_cols = ["chrom", "start", "end", "RAD_tag", "RAD_pos", "snp_id", "ref", "alt", "strand"]
snps = pd.read_csv(snps_file, sep="\t", names=snps_cols, skiprows=1, dtype={"start": int, "end": int})

# Ensure 'start' is int (in case file has quotes or weird formatting)
snps["start"] = snps["start"].astype(int)

# --- Load GTF ---
gtf_file = "AnoSag2.1_cds.gtf"
gtf_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=gtf_cols)

# Convert start/end to int
gtf["start"] = gtf["start"].astype(int)
gtf["end"] = gtf["end"].astype(int)

# Extract transcript_id
gtf["transcript_id"] = gtf["attribute"].apply(lambda x: parse_attributes(x).get("transcript_id"))

# Filter only exons
exons = gtf[gtf["feature"] == "exon"]

# --- Group exon intervals and compute introns ---
intron_regions = []
grouped = exons.groupby("transcript_id")

for tid, group in grouped:
    strand = group["strand"].iloc[0]
    chrom = group["chrom"].iloc[0]
    
    # Sort exons by genomic coordinate
    sorted_exons = group.sort_values("start")
    exon_intervals = list(zip(sorted_exons["start"], sorted_exons["end"]))
    
    # Compute introns as gaps between exons
    for i in range(len(exon_intervals) - 1):
        intron_start = exon_intervals[i][1] + 1
        intron_end = exon_intervals[i + 1][0] - 1
        if intron_start <= intron_end:
            intron_regions.append((chrom, intron_start, intron_end, strand, tid))

# --- Check SNPs in introns ---
results = []
for idx, snp in snps.iterrows():
    chrom = snp["chrom"]
    pos = int(snp["start"])
    for region in intron_regions:
        r_chrom, r_start, r_end, r_strand, r_tid = region
        if chrom == r_chrom and r_start <= pos <= r_end:
            results.append({
                "snp_id": snp["snp_id"],
                "chrom": chrom,
                "pos": pos,
                "ref": snp["ref"],
                "alt": snp["alt"],
                "strand": r_strand,
                "transcript_id": r_tid
            })
            break  # Found a match

# --- Save results ---
df_results = pd.DataFrame(results)
df_results.to_csv("snps_in_introns.tsv", sep="\t", index=False)

print(f"Identified {len(df_results)} SNPs in intron regions.")