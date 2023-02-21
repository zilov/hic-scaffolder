rule yahs:
    input:
        bam = BAM_TO_SCAFFOLD,
        assembly = GENOME,
    conda: envs.yahs
    threads: workflow.cores
    output: 
        f"{OUTDIR}/yahs/{PREFIX}_scaffolds_final.fa"
    params:
        prefix = PREFIX,
        yahs_dir = directory(f"{OUTDIR}/yahs")
    shell: "yahs {input.assembly} {input.bam} -o {params.yahs_dir}/{params.prefix}"