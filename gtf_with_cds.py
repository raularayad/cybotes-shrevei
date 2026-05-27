#Creates gtf including CDS data

from collections import defaultdict

input_gtf = "AnoSag2.1.gtf"
output_gtf = "AnoSag2.1_cds.gtf"

transcripts = defaultdict(lambda: {
    'strand': None,
    'chrom': None,
    'exons': [],
    'start_codon': None,
    'stop_codon': None,
    'attrs': None,
})

# Step 1: Parse original GTF
with open(input_gtf) as f:
    gtf_lines = f.readlines()

for line in gtf_lines:
    if line.startswith("#") or not line.strip():
        continue
    fields = line.strip().split('\t')
    chrom, source, feature, start, end, score, strand, phase, attr_str = fields
    start, end = int(start), int(end)

    attrs = dict(
        kv.strip().replace('"', '').split(' ', 1)
        for kv in attr_str.strip().split(';') if kv.strip()
    )
    transcript_id = attrs.get('transcript_id')
    if not transcript_id:
        continue

    tx = transcripts[transcript_id]
    tx['chrom'] = chrom
    tx['strand'] = strand
    tx['attrs'] = attrs

    if feature == 'exon':
        tx['exons'].append((start, end))
    elif feature == 'start_codon':
        tx['start_codon'] = start if strand == '+' else end
    elif feature == 'stop_codon':
        tx['stop_codon'] = end if strand == '+' else start

# Step 2: Write original GTF + inferred CDS entries
with open(output_gtf, "w") as out:
    out.writelines(gtf_lines)  # Write original content

    for tid, data in transcripts.items():
        if not data['start_codon'] or not data['stop_codon']:
            continue  # Skip incomplete

        chrom = data['chrom']
        strand = data['strand']
        exons = sorted(data['exons'])

        cds_start = data['start_codon'] if strand == '+' else data['stop_codon']
        cds_end = data['stop_codon'] if strand == '+' else data['start_codon']

        for exon_start, exon_end in exons:
            overlap_start = max(cds_start, exon_start)
            overlap_end = min(cds_end, exon_end)
            if overlap_start <= overlap_end:
                attrs = data['attrs']
                attr_string = ' '.join([f'{k} "{v}";' for k, v in attrs.items()])
                out.write(f"{chrom}\tAUGUSTUS\tCDS\t{overlap_start}\t{overlap_end}\t.\t{strand}\t0\t{attr_string}\n")