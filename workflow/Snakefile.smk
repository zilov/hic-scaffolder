from pathlib import Path, PurePath

OUTDIR = config["outdir"]
GENOME = config["assembly"]
FORWARD_READ = config["forward_hic_read"]
REVERSE_READ = config["reverse_hic_read"]
BAM = config["bam"]
REPLICATES = config["replicates"]
MODE = config["mode"]
GENOME_DIR = Path(GENOME).parent.resolve()
EXE_FOLDER = config["execution_folder"]

PREFIX = PurePath(GENOME).stem

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

if MODE == 'align_hic_pair':
    if not Path(f"{OUTDIR}/tmp/reads").exists():
        Path(f"{OUTDIR}/tmp/reads").mkdir()
    if not Path(f"{OUTDIR}/tmp/reads/{PREFIX}_1.fq").exists():
        Path(f"{OUTDIR}/tmp/reads/{PREFIX}_1.fq").symlink_to(FORWARD_READ)
    if not Path(f"{OUTDIR}/tmp/reads/{PREFIX}_2.fq").exists():
        Path(f"{OUTDIR}/tmp/reads/{PREFIX}_2.fq").symlink_to(REVERSE_READ)
    READ_PREFIXES = [f"{PREFIX}_1", f"{PREFIX}_2"]

rule all:
    input:
         f"{OUTDIR}/quast/report.txt",
         f"{OUTDIR}/juicer/{PREFIX}.hic",
         f"{OUTDIR}/juicer_JBAT/{PREFIX}_JBAT.hic",

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

include: "./rules/map_reads.smk"

include: "./rules/filter_bams.smk"

include: "./rules/merge_bams.smk"

include: "./rules/add_read_group.smk"

include: "./rules/mark_duplicates.smk"

BAM_TO_SCAFFOLD = rules.mark_duplicates.output.bam_with_duplicate

include: "./rules/yahs.smk"

SCAFFOLDS = rules.yahs_scaffold.output.scaffolds

include: "./rules/quast.smk"

include: "./rules/juicer_maps.smk"