import pandas as pd

# === 1. Load SNP data ===
snp_file = "snp_genome_coords.tsv"
snp_cols = ["chrom", "start", "end", "RAD_tag", "RAD_pos", "snp_id", "ref", "alt", "strand"]
snps = pd.read_csv(snp_file, sep="\t", names=snp_cols)

# === 2. Load GTF file and define columns ===
gtf_file = "AnoSag2.1_cds.gtf"
gtf = pd.read_csv(gtf_file, sep="\t", header=None, comment="#")
gtf.columns = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# === 3. Extract transcript_id ===
def extract_transcript_id(attr):
    for part in attr.split(";"):
        if "transcript_id" in part:
            return part.split('"')[1]
    return None

gtf["transcript_id"] = gtf["attribute"].apply(extract_transcript_id)

# === 4. Filter CDS entries ===
cds_regions = gtf[gtf["feature"] == "CDS"].copy()

# === 5. Identify SNPs in CDS regions ===
snp_in_cds = []
for _, snp in snps.iterrows():
    matches = cds_regions[
        (cds_regions["chrom"] == snp["chrom"]) &
        (cds_regions["start"] <= snp["start"]) &
        (cds_regions["end"] >= snp["start"])
    ]
    for _, match in matches.iterrows():
        snp_in_cds.append([
            snp["chrom"], snp["start"], snp["end"], snp["snp_id"],
            snp["ref"], snp["alt"], snp["strand"], match["transcript_id"]
        ])

# === 6. Output results ===
snp_cds_df = pd.DataFrame(snp_in_cds, columns=[
    "chrom", "start", "end", "snp_id", "ref", "alt", "strand", "transcript_id"
])
snp_cds_df.drop_duplicates(inplace=True)
snp_cds_df.to_csv("snp_in_CDS.tsv", sep="\t", index=False)