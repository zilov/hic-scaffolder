from pathlib import Path, PurePath

OUTDIR = config["outdir"]
GENOME = config["assembly"]
FORWARD_READ = config["forward_hic_read"]
REVERSE_READ = config["reverse_hic_read"]
BAM = config["bam"]
REPLICATES = config["replicates"]
MODE = config["mode"]
THREADS = ["threads"]
GENOME_DIR = Path(GENOME).parent.resolve()

PREFIX = PurePath(GENOME).stem

## create symlink to genome to not create index files in genome folder

if not Path(f"{OUTDIR}/tmp/genome").exists():
    Path(f"{OUTDIR}/tmp/genome").mkdir()
if not Path(f"{OUTDIR}/tmp/genome/{PREFIX}.fasta").exists():
    Path(f"{OUTDIR}/tmp/genome/{PREFIX}.fasta").symlink_to(GENOME)
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
        expand('{OUTDIR}/tmp/bams/{sample}_filtered.bam', sample=READ_PREFIXES, OUTDIR=OUTDIR)

rule envs:
    params:
        arima = "../envs/arima_pipeline.yaml",
        yahs = "../envs/yahs.yaml"
    
envs = rules.envs.params

rule scripts:
    params:
        filtering = Path("./scripts/mapping_pipeline/filter_five_end.pl").resolve(),
        merge = Path("./scripts/mapping_pipeline/two_read_bam_combiner.pl").resolve(),
        stats = Path("./scripts/mapping_pipeline/get_stats.pl").resolve(),
    
scripts = rules.scripts.params

include: "./rules/map_reads.smk"

include: "./rules/filter_bams.smk"