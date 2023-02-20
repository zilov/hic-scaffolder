rule mark_duplicates:
    input:
        bam_with_readgroup = rules.add_read_group.output
    conda: envs.arima
    output:
        bam_with_duplicate = f'{OUTDIR}/tmp/bams/{PREFIX}_marked_duplicates.bam',
        metrics = f'{OUTDIR}/tmp/bams/{PREFIX}_marked_duplicates.metrics'
    params:
        tmpdir = f"{OUTDIR}/tmp/bams/tmp"
    shell: "picard MarkDuplicates INPUT={input} OUTPUT={output.bam_with_duplicate} METRICS_FILE={output.metrics} \
            TMP_DIR={params.tmpdir} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE"