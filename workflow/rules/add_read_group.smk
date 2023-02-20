rule add_read_group:
    input: 
        merged_bam = rules.merge_bams.output
    conda: envs.arima
    output: 
        bam_with_readgroup = f'{OUTDIR}/tmp/bams/{PREFIX}_with_readgroup.bam'
    params: 
        prefix = PREFIX
    shell: "picard AddOrReplaceReadGroups INPUT={input} OUTPUT={output} \
            ID={params.prefix} LB={params.prefix} SM={params.prefix} PL=ILLUMINA PU=none"
    