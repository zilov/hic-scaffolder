import glob

rule review_to_fasta: 
    input:
        review_assembly = REVIEWED_ASSEMBLY,
        run_folder = RUN_FOLDER
    output:
        reviewed_fasta = f"{OUTDIR}/{PREFIX}_reviewed.FINAL.fa",
    threads: workflow.cores
    conda: envs.yahs
    params:
        liftover = glob.glob(f"{RUN_FOLDER}/juicer_JBAT/*.liftover.agp")[0],
        assembly = glob.glob(f"{RUN_FOLDER}/tmp/genome/*.fasta")[0],
        prefix = f"{OUTDIR}/{PREFIX}_reviewed",
        juicer_dir = directory(f"{OUTDIR}/juicer_JBAT"),
    shell:"""
    juicer post -o {params.prefix} {input.review_assembly} {params.liftover} {params.assembly} 
    """