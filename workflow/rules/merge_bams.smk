rule merge_bams:
    input:
        reads = expand('{outdir}/tmp/bams/{sample}_filtered.bam', outdir=OUTDIR, sample=READ_PREFIXES),
        genome_index = rules.samtools_index_assembly.output,
    conda: envs.arima
    output:
        merged_bam = f"{OUTDIR}/tmp/bams/{PREFIX}_merged.bam"
    params: 
        merge_script = scripts.merge
    shell: "perl {params.merge_script} {input.reads} samtools 10 | \
            samtools view -bS -t {input.genome_index} - | \
            samtools sort -@ {threads} -o {output.merged_bam}"
    