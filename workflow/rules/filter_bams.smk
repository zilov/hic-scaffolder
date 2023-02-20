rule filter_bams:
    input: 
        bam_read = rules.map_reads.output,
    conda: envs.arima
    output:
        bam_filtered = '%s/tmp/bams/{sample}_filtered.bam' % OUTDIR,
    params:
        filter_script = scripts.filtering
    shell: "samtools view -h {input} | perl {params.filter_script} | samtools view -Sb - > {output}"
    