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
    shell: "yahs {input.assembly} {input.bam} -o {params.yahs_dir}/{params.prefix}"

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