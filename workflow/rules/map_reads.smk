rule samtools_index_assembly:
    input: GENOME
    conda: envs.arima
    threads: workflow.cores * 0.5
    output: 
        samtools_index = f"{GENOME}.fai"
    shell: "samtools faidx {input}"


rule bwa_index_assembly:
    input: GENOME
    conda: envs.arima
    threads: workflow.cores * 0.5
    output:
        bwa_index = f"{GENOME}.ann"
    shell: "bwa index -a bwtsw {input}"

rule map_reads:
    input:
        read_fq = '%s/tmp/reads/{sample}.fq' % OUTDIR,
        genome = GENOME,
        samtools_index = rules.samtools_index_assembly.output.samtools_index,
        bwa_index_assembly = rules.bwa_index_assembly.output.bwa_index,
    conda: envs.arima
    threads: workflow.cores
    output:
        read_bam = '%s/tmp/bams/{sample}.bam' % OUTDIR
    params: 
        outdir = f'{OUTDIR}/tmp/bams/'
    shell: "bwa mem -t {threads} {input.genome} {input.read_fq} | samtools view -@ {threads} -Sb - > {output.read_bam}"
    

