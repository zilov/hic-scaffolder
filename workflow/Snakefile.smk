from pathlib import Path, PurePath

OUTDIR = config["outdir"]
GENOME = config["assembly_fasta"]
FORWARD_READ = config["forward_hic_read"]
REVERSE_READ = config["reverse_hic_read"]
BAM = config["bam"]
MODE = config["mode"]
GENOME_DIR = Path(GENOME).parent.resolve()
EXE_FOLDER = config["execution_folder"]
SCAFFOLDER = config["scaffolder"]
REVIEWED_ASSEMBLY = config["reviewed_assembly"]
RUN_FOLDER = config['run_folder']
PREFIX = config["prefix"]

def fasta_reader_yield(path_to_fasta_file):
    header = None
    with open(path_to_fasta_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header,"".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header,"".join(seq)


if MODE == 'assembly_to_hic_map':
    ## copy genome to tmp folder without extra-short contigs < 500bp

    if not Path(f"{OUTDIR}/tmp/genome").exists():
        Path(f"{OUTDIR}/tmp/genome").mkdir(parents=True)
    if not Path(f"{OUTDIR}/tmp/genome/{PREFIX}.fasta").exists():
        print("Extracting contigs with length > 500bp")
        with open(f"{OUTDIR}/tmp/genome/{PREFIX}.fasta", 'w') as fw:
            for h, s in fasta_reader_yield(GENOME):
                if len(s) > 500:
                    fw.write(f">{h}\n{s}\n")
        print('Done!')
    GENOME = f"{OUTDIR}/tmp/genome/{PREFIX}.fasta"


    if not Path(f"{OUTDIR}/tmp/reads").exists():
        Path(f"{OUTDIR}/tmp/reads").mkdir()
    if not Path(f"{OUTDIR}/tmp/reads/{PREFIX}_1.fq").exists():
        Path(f"{OUTDIR}/tmp/reads/{PREFIX}_1.fq").symlink_to(FORWARD_READ)
    if not Path(f"{OUTDIR}/tmp/reads/{PREFIX}_2.fq").exists():
        Path(f"{OUTDIR}/tmp/reads/{PREFIX}_2.fq").symlink_to(REVERSE_READ)
    READ_PREFIXES = [f"{PREFIX}_1", f"{PREFIX}_2"]

if MODE == "assembly_to_hic_map":
    if SCAFFOLDER == "yahs":
        OUTPUTS = [f"{OUTDIR}/juicer/{PREFIX}.hic", f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.hic"]
elif MODE == "reviewed_to_fasta":
    OUTPUTS = [f"{OUTDIR}/{PREFIX}_reviewed.FINAL.fa"]

OUTPUTS.append(f"{OUTDIR}/quast/report.txt")

rule all:
    input:
         OUTPUTS

rule envs:
    params:
        arima = "../envs/arima_pipeline.yaml",
        yahs = "../envs/yahs.yaml",
        quast = "../envs/quast.yaml"
    
envs = rules.envs.params

rule scripts:
    params:
        filtering = f"{EXE_FOLDER}/scripts/mapping_pipeline/filter_five_end.pl",
        merge = f"{EXE_FOLDER}/scripts/mapping_pipeline/two_read_bam_combiner.pl",
        stats = f"{EXE_FOLDER}/scripts/mapping_pipeline/get_stats.pl",
        juicer_tools = f"{EXE_FOLDER}/scripts/juicer_tools.jar"
    
scripts = rules.scripts.params

if MODE == "assembly_to_hic_map":
    include: "./rules/map_reads.smk"
    include: "./rules/filter_bams.smk"
    include: "./rules/merge_bams.smk"
    include: "./rules/add_read_group.smk"
    include: "./rules/mark_duplicates.smk"
    BAM_TO_SCAFFOLD = rules.mark_duplicates.output.bam_with_duplicate
    if SCAFFOLDER == 'yahs':
        include: "./rules/yahs.smk"
        SCAFFOLDS = rules.yahs_scaffold.output.scaffolds
        include: "./rules/juicer_maps.smk"
    elif SCAFFOLDER == 'salsa':
        include: "./rules/salsa.smk"
        SCAFFOLDS = rules.salsa_scaffold.output.scaffolds
        
elif MODE == "reviewed_to_fasta":
    include: "./rules/review_to_fasta.smk"
    SCAFFOLDS = rules.review_to_fasta.output.reviewed_fasta

include: "./rules/quast.smk"
