import pandas as pd

def parse_attributes(attr_str):
    attrs = {}
    for attr in attr_str.strip().split(";"):
        if attr.strip():
            key, val = attr.strip().split(" ", 1)
            attrs[key] = val.strip('"')
    return attrs

def merge_intervals(intervals):
    if not intervals:
        return []
    sorted_ints = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_ints[0]]
    for current in sorted_ints[1:]:
        last = merged[-1]
        if current[0] <= last[1] + 1:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def subtract_intervals(minuend, subtrahend):
    """Subtract intervals in subtrahend from minuend intervals."""
    result = []
    for start, end in minuend:
        temp_start = start
        for s_start, s_end in subtrahend:
            if s_end < temp_start:
                continue
            if s_start > end:
                break
            if s_start <= temp_start <= s_end:
                temp_start = s_end + 1
            elif temp_start < s_start <= end:
                result.append((temp_start, s_start - 1))
                temp_start = s_end + 1
        if temp_start <= end:
            result.append((temp_start, end))
    return result

# --- Load SNPs ---

snps_file = "snp_genome_coords.tsv"
snps_cols = ["chrom", "start", "end", "RAD_tag", "RAD_pos", "snp_id", "ref", "alt", "strand"]
snps = pd.read_csv(snps_file, sep="\t", names=snps_cols)

# --- Load GTF ---

gtf_file = "AnoSag2.1_cds.gtf"
gtf_cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=gtf_cols)

# Extract transcript_id from attribute column
gtf["transcript_id"] = gtf["attribute"].apply(lambda x: parse_attributes(x).get("transcript_id"))

# Filter exons and CDS
exons = gtf[gtf["feature"] == "exon"]
cds = gtf[gtf["feature"] == "CDS"]

# Group intervals by transcript_id
exons_by_transcript = {}
for tid, group in exons.groupby("transcript_id"):
    intervals = list(zip(group["start"], group["end"]))
    exons_by_transcript[tid] = merge_intervals(intervals)

cds_by_transcript = {}
for tid, group in cds.groupby("transcript_id"):
    intervals = list(zip(group["start"], group["end"]))
    cds_by_transcript[tid] = merge_intervals(intervals)

# Get strand info per transcript (assumed consistent)
strand_by_transcript = gtf.groupby("transcript_id")["strand"].first().to_dict()
chrom_by_transcript = gtf.groupby("transcript_id")["chrom"].first().to_dict()

# Compute UTRs (exon minus CDS) per transcript
utr_by_transcript = {}
for tid in exons_by_transcript:
    exon_intervals = exons_by_transcript[tid]
    cds_intervals = cds_by_transcript.get(tid, [])
    utr_intervals = subtract_intervals(exon_intervals, cds_intervals)
    utr_by_transcript[tid] = utr_intervals

# Determine 5' UTR intervals per transcript based on strand
five_utr_regions = []
for tid, utrs in utr_by_transcript.items():
    strand = strand_by_transcript.get(tid)
    chrom = chrom_by_transcript.get(tid)
    cds_intervals = cds_by_transcript.get(tid, [])

    if not strand or not cds_intervals:
        continue

    # Determine the 5' end of CDS for this transcript
    if strand == "+":
        cds_5prime = min(e[0] for e in cds_intervals)
        # 5' UTR is UTR region **before** CDS start (lower coordinates)
        five_prime_utrs = [(start, end) for start, end in utrs if end < cds_5prime]
    else:  # strand == "-"
        cds_5prime = max(e[1] for e in cds_intervals)
        # 5' UTR is UTR region **after** CDS end (higher coordinates)
        five_prime_utrs = [(start, end) for start, end in utrs if start > cds_5prime]

    for start, end in five_prime_utrs:
        five_utr_regions.append((chrom, start, end, strand, tid))

# Now check which SNPs fall in 5' UTR regions

results = []
for idx, snp in snps.iterrows():
    chrom = snp["chrom"]
    pos = snp["start"]
    for region in five_utr_regions:
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
            break  # SNP assigned to first matching region

# Save results
df_results = pd.DataFrame(results)
df_results.to_csv("snps_in_5UTR.tsv", sep="\t", index=False)

print(f"Identified {len(df_results)} SNPs in 5' UTR regions.")