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

def parse_replicates(table):
    result = {}
    with open(table) as fh:
        header = fh.readline().strip().split(",")
        for row in fh:
            row = row.strip().split(',')
            sample = row[0]
            biological_replicate = row[1]
            technical_replicate = row[2]
            value1 = row[3]
            value2 = row[4]

            if sample not in result:
                result[sample] = {}

            if biological_replicate not in result[sample]:
                result[sample][biological_replicate] = {}
                
            if technical_replicate not in result[sample][biological_replicate]:
                result[sample][biological_replicate][technical_replicate] = (value1, value2)
        
        return result

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

elif MODE == 'align_hic_replicates':
    SAMPLES = []
    BIO_REPS = []
    TECH_REPS = []
    REPLICATES = parse_replicates(REPLICATES)
    for sample, bio_repl in REPLICATES.items():
        SAMPLES.append(sample) if sample not in SAMPLES
        for bio_repl_n, tech_repl in bio_repl.items():
            BIO_REPS.append(f"{sample}_{bio_repl_n}")
            for tech_repl_n, reads in tech_repl.items():
                TECH_REPS.append(f"{sample}_{bio_repl_n}_{tech_repl_n}")
                Path(f"{OUTDIR}/tmp/reads/{sample}_{bio_repl_n}_{tech_repl_n}_1.fq").symlink_to(reads[0])
                Path(f"{OUTDIR}/tmp/reads/{sample}_{bio_repl_n}_{tech_repl_n}_2.fq").symlink_to(reads[1])

rule all:
    input:
       f'{OUTDIR}/tmp/bams/{PREFIX}_marked_duplicates.bam'

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

include: "./rules/merge_bams.smk"

include: "./rules/add_read_group.smk"

include: "./rules/mark_duplicates.smk"