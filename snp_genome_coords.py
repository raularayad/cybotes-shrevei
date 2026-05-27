from Bio import SearchIO

# STEP 1: Load SNPs into a dictionary {RAD_tag: [positions]}
snp_dict = {}
with open("variants_noprefix.tsv") as f:
    next(f)  # skip header
    for line in f:
        rad_tag, pos, var_id, ref, alt = line.strip().split('\t')[:5]
        snp_dict.setdefault(rad_tag, []).append((int(pos), var_id, ref, alt))

# STEP 2: Parse BLAST XML and map SNPs to genomic coordinates
output = []

for qresult in SearchIO.parse("blast_results.xml", "blast-xml"):
    qid = qresult.id
    if qid not in snp_dict:
        continue
    if not qresult.hits:
        print(f"# Warning: No hits for {qid}")
        continue
    hit = qresult.hits[0]  # Best hit
    hsp = hit.hsps[0]      # First HSP (alignment)
    
    query_seq = hsp.query.seq
    hit_seq = hsp.hit.seq

    q_pos = hsp.query_start  # 0-based
    h_pos = hsp.hit_start    # 0-based
    strand = hsp.hit_strand  # +1 or -1

    query_to_genome = {}  # {query_pos (1-based): genome_pos (0-based)}

    q_counter = q_pos
    h_counter = h_pos

    for i in range(len(query_seq)):
        q_char = query_seq[i]
        h_char = hit_seq[i]

        if q_char != '-':
            q_real = q_counter + 1  # convert to 1-based for SNP input
            if h_char != '-':
                query_to_genome[q_real] = h_counter
                h_counter += 1
            else:
                query_to_genome[q_real] = None  # gap in genome
            q_counter += 1
        elif h_char != '-':
            h_counter += 1

    for snp_pos, var_id, ref, alt in snp_dict[qid]:
        genome_coord = query_to_genome.get(snp_pos)
        if genome_coord is not None:
            chrom = hit.id
            if strand == -1:
                # Reverse strand: return reverse-complemented alleles if needed
                # genome_coord = real position; still valid
                output.append((chrom, genome_coord, genome_coord + 1, qid, snp_pos, var_id, ref, alt, "-"))
            else:
                output.append((chrom, genome_coord, genome_coord + 1, qid, snp_pos, var_id, ref, alt, "+"))
        else:
            print(f"# Warning: SNP {var_id} at pos {snp_pos} in {qid} maps to a gap")

# STEP 3: Write output
with open("snp_genome_coords.tsv", "w") as out:
    out.write("chrom\tstart\tend\tRAD_tag\tRAD_pos\tID\tREF\tALT\tstrand\n")
    for row in output:
        out.write('\t'.join(map(str, row)) + '\n')