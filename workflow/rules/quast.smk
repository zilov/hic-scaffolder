rule quast:
    input:
        SCAFFOLDS,
    conda:
        envs.quast,
    output:
        quast_output = "{OUTDIR}/quast/report.txt"
    params:
        directory("{OUTDIR}/quast")
    shell: "quast -o {params} {input}"