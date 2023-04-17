rule juicer_alignment: 
    input:
        alignment_bin = rules.yahs_scaffold.output.alignment_bin,
        scaffolds_agp = rules.yahs_scaffold.output.scaffolds_agp,
        contigs_index = rules.samtools_index_assembly.output.samtools_index,
    output:
        alignments_sorted = f"{OUTDIR}/juicer/{PREFIX}_alignment_sorted.txt",
    threads: workflow.cores
    conda: envs.yahs
    params:
        juicer_dir = directory(f"{OUTDIR}/juicer"),
        tmp_alignments = f"{OUTDIR}/juicer/{PREFIX}_alignment_sorted.txt.tmp"
    shell:"""
    (juicer pre {input.alignment_bin} {input.scaffolds_agp} {input.contigs_index} |\
      sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > {params.tmp_alignments}) && \
      (mv {params.tmp_alignments} {output.alignments_sorted})
    """

rule juicer_hic:
    input: 
        alignment = rules.juicer_alignment.output.alignments_sorted,
        scaffolds_chrom_sizes = rules.scaffold_index.output.scaffolds_chrom_sizes,
    output:
        hic_map = f"{OUTDIR}/juicer/{PREFIX}.hic"
    conda: envs.yahs
    params:
        juicer_tools = scripts.juicer_tools,
        hic_map_tmp = f"{OUTDIR}/juicer/{PREFIX}.hic.tmp"

    shell: """
    ({params.juicer_tools} pre {input.alignment} {params.hic_map_tmp} {input.scaffolds_chrom_sizes}) && \
    (mv {params.hic_map_tmp} {output.hic_map})
    """


rule juicer_JBAT: 
    input:
        alignment_bin = rules.yahs_scaffold.output.alignment_bin,
        scaffolds_agp = rules.yahs_scaffold.output.scaffolds_agp,
        contigs_index = rules.samtools_index_assembly.output.samtools_index,
    output:
        bed_JBAT = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.txt",
        assembly_JBAT = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.assembly",
        log_JBAT = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.log",
    threads: workflow.cores
    conda: envs.yahs
    params:
        prefix = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT",
        juicer_dir = directory(f"{OUTDIR}/juicer_JBAT"),
    shell:"""
    juicer pre -a -o {params.prefix} {input.alignment_bin} {input.scaffolds_agp} \
    {input.contigs_index} >{output.log_JBAT} 2>&1
    """

rule juicer_JBAT_hic:
    input: 
        bed_JBAT = rules.juicer_JBAT.output.bed_JBAT,
        log_JBAT = rules.juicer_JBAT.output.log_JBAT
    output:
        hic_map_JBAT = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.hic"
    conda: envs.yahs
    params:
        juicer_tools = scripts.juicer_tools,
        hic_map_tmp = f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.hic.tmp"
    shell: """
    ({params.juicer_tools} pre {input.bed_JBAT} {params.hic_map_tmp} \
    <(cat {input.log_JBAT}  | grep PRE_C_SIZE | awk '{{print $2" "$3}}')) && \
    (mv {params.hic_map_tmp} {output.hic_map_JBAT})
    """
