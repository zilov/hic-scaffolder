# hic-scaffolder: A Snakemake workflow for Hi-C scaffolding and analysis
This repository provides a Snakemake workflow for processing Hi-C data, scaffolding genome assemblies, and generating Hi-C interaction maps.

## Key Features:
**Two Modes of Operation:** 

1. *Assembly to Hi-C map:* This mode takes an assembly FASTA file and paired-end Hi-C reads (FASTQ or BAM format) as input and performs the following steps:
- **Read mapping and filtering:** Aligns Hi-C reads to the assembly and filters out low-quality alignments and duplicates using [bwa-mem](https://github.com/bwa-mem2/bwa-mem2), [samtols](https://github.com/samtools/samtools) and [picard](https://github.com/broadinstitute/picard) ([Arima mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline) standart).
- **Scaffolding:** Utilizes either [YaHS](https://github.com/c-zhou/yahs) or [Salsa](https://github.com/marbl/SALSA) (TODO) for scaffolding based on Hi-C interaction information.
- **Hi-C map generation:** Creates Hi-C interaction maps using [Juicer tools](https://github.com/aidenlab/JuicerTools) ready to edit in [Juicebox](https://github.com/aidenlab/Juicebox). Juicebox version 1.11.08 is recommended.
2. *Reviewed to FASTA:* This mode converts reviewed assembly files from previous hic-scaffolder runs into FASTA format for downstream analysis.
   
**Supported Scaffolders:**
- YaHS (default)
- Salsa (TODO)

**Quality Control:** Employs [quast](https://github.com/ablab/quast) for assembly quality assessment.

**Reproducibility:** Leverages conda environments for dependency management ensuring consistent and reproducible results. Dont necessary to install all tools, only snakemake environment requeired.

## Getting Started:
1. Clone the repository:
```git clone https://github.com/your-username/hic-scaffolder.git```
2. Install dependencies:
- Install conda or mamba.
- Create snakemake conda environment `conda/mamba create -n snakemake -c bioconda -c conda forge -c defaults snakemake`
3. Run the pipeline:
- Assembly to Hi-C map mode, to get scaffolds:

```python hic-scaffolder.py -a assembly.fasta -1 forward_reads.fq.gz -2 reverse_reads.fq.gz -o output_directory```

`.hic` and `.assembly` files for Juicebox manual curation of genome lays in `output_directory/juicer_JBAT` folder. After editing of scaffolds use `*_reviewed.assembly` file to get fasta of resulted assembly. 
- Reviewed to FASTA mode:
   
```python hic-scaffolder.py -m reviewed_to_fasta -o output_directory -f hic-scaffolder_run_folder -r reviewed_assembly.txt```

## Options:
```
| Option                | Description                                                            | Required                                            |
|-----------------------|------------------------------------------------------------------------|-----------------------------------------------------|
| -o/--outdir           | Path to the output directory                                           | Yes                                                 |
| -m/--mode             | Mode of operation (assembly_to_hic_map or reviewed_to_fasta)           | Yes                                                 |
| -a/--assembly         | Path to the assembly FASTA file                                        | assembly_to_hic_map mode                            |
| -1/--forward_hic_read | Path to the forward Hi-C reads file                                    | assembly_to_hic_map mode (if not using -b)          |
| -2/--reverse_hic_read | Path to the reverse Hi-C reads file                                    | assembly_to_hic_map mode (if not using -b)          |
| -b/--bam              | Path to a BAM file containing aligned Hi-C reads                       | assembly_to_hic_map mode (alternative to -1 and -2) |
| -r/--reviewed         | Path to the reviewed assembly file                                     | reviewed_to_fasta mode                              |
| -f/--run_folder       | Path to the hic-scaffolder run folder containing the reviewed assembly | reviewed_to_fasta mode                              |
| -t/--threads          | Number of threads to use                                               | No (default: 8)                                     |
| -d/--debug            | Activate debug mode                                                    | No                                                  |
| -s/--scaffolder       | Scaffolding tool (yahs or salsa)                                       | No (default: yahs)                                  |
```