rule yahs_scaffold:
    input:
        bam = BAM_TO_SCAFFOLD,
        assembly = GENOME,
    conda: envs.yahs
    threads: workflow.cores
    output: 
        scaffolds = f"{OUTDIR}/yahs/{PREFIX}_scaffolds_final.fa",
        scaffolds_agp = f"{OUTDIR}/yahs/{PREFIX}_scaffolds_final.agp",
        alignment_bin = f"{OUTDIR}/yahs/{PREFIX}.bin",
    params:
        prefix = PREFIX,
        yahs_dir = directory(f"{OUTDIR}/yahs")
    shell: "yahs {input.assembly} {input.bam} -q 40 \
    -r 1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000,10000000,20000000,50000000,100000000,200000000,500000000 \
    -o {params.yahs_dir}/{params.prefix}"

rule scaffold_index:
    input:
        scaffolds = rules.yahs_scaffold.output.scaffolds,
    conda: envs.arima
    output:
        scaffolds_index = f"{OUTDIR}/yahs/{PREFIX}_scaffolds_final.fa.fai",
        scaffolds_chrom_sizes = f"{OUTDIR}/yahs/{PREFIX}_scaffolds_final.chrom_sizes",
    shell: """
    samtools faidx {input.scaffolds}

    awk -F'\\t' '{{print $1"\\t"$2}}' {output.scaffolds_index} > {output.scaffolds_chrom_sizes}
    """