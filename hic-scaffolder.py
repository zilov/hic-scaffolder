#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.02.2023
#@author: Danil Zilov
#@contact: zilov.d@gmail.com

import argparse
import os
import os.path
from inspect import getsourcefile
from datetime import datetime
import string
import random

def config_maker(settings, config_file):
    config = f"""
    "outdir" : "{settings["outdir"]}"
    "assembly" : "{settings["assembly_fasta"]}"
    "forward_read" : "{settings["fr"]}"
    "reverse_read" : "{settings["rr"]}"
    "bams": "{settings["bams"]}"
    "mode" : "{settings["mode"]}"
    "threads" : "{settings["threads"]}"
    """

    if not os.path.exists(os.path.dirname(config_file)):
        os.mkdir(os.path.dirname(config_file))


    with open(config_file, "w") as fw:
        fw.write(config)
        print(f"CONFIG IS CREATED! {config_file}")
      

def main(settings):
    ''' Function description.
    '''
        
    if settings["debug"]:
        snake_debug = "-n"
    else:
        snake_debug = ""

    #Snakemake
    command = f"""
    snakemake --snakefile {settings["execution_folder"]}/workflow/Snakefile \
              --configfile {settings["config_file"]} \
              --cores {settings["threads"]} \
              --use-conda
              --conda-frontend mamba
              --conda-prefix {settings["execution_folder"]}/.conda_envs
              {snake_debug}"""
    print(command)
    os.system(command)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parrots pipeline - a tool for parrots genome annotation project.')
    parser.add_argument('-m','--mode', help="mode to use [default = fasta]", 
                        choices=["align_hic", "merge_tech_reps", "merge_bio_reps", "yahs"], default="align_hic")
    parser.add_argument('-a','--assembly', help="path to asssembly fasta file", type=argparse.FileType('r'), default="")
    parser.add_argument('-1','--forward_hix_read', help="path to forward hic read", type=argparse.FileType('r'), default="")
    parser.add_argument('-2','--reverse_hic_read', help="path to reverse hic read", type=argparse.FileType('r'), default="")
    parser.add_argument('--bams-tech-rep', 
                        help='technical replicates bam files to merge, space-separated list',
                        type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--bio-rep-bams',
                        help='technical replicates bam files to merge, space-separated list',
                        type=argparse.FileType('r'), nargs='+')
    parser.add_argument('-b', '--bam', help='path to bam pile to run YaHS', type=argparse.FileType('r'), default='')
    parser.add_argument('-o','--outdir', help='output directory', required=True)
    parser.add_argument('-t','--threads', help='number of threads [default == 8]', default = "8")
    parser.add_argument('-d','--debug', help='debug mode', action='store_true')
    args = vars(parser.parse_args())

    assembly_fasta = os.path.abspath(args["assembly"])
    threads = args["threads"]
    debug = args["debug"]
    mode = args["mode"]
    forward_hic_read = args["forward_hic_read"]
    reverse_hic_read = args["reverse_hic_read"]
    bams_tech_rep = args['bams-tech-rep']
    bams_bio_rep = args['bams-bio-rep']
    bam = args['bam']
    outdir = os.path.abspath(args["outdir"])
    
    execution_folder = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    execution_time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    config_random_letters = "".join([random.choice(string.ascii_letters) for n in range(3)])
    config_file = os.path.join(execution_folder, f"config/config_{config_random_letters}{execution_time}.yaml")
    # default braker_file for fasta mode
    braker_file = os.path.join(execution_folder,"rules/braker_fasta.smk")

    if mode == "align_hic":
        if not (forward_hic_read or reverse_hic_read):
            parser.error("\nalign_hic mode requires -1 {path_to_forward_read} and -2 {path_to_reverse_read}!")
        else:
            forward_hic_read = os.path.abspath(forward_hic_read)
            reverse_hic_read = os.path.abspath(reverse_hic_read)    
    elif mode == "merge_tech_reps":
        if not bams_tech_rep:
            parser.error("\nmerge_tech_reps mode requires --bams-tech-rep tech_rep1.bam tech_rep2.bam")
        else:
            bams_tech_rep = [os.path.abspath(x) for x in bams_tech_rep]
    elif mode == "merge_bio_reps":
        if not bams_bio_rep:
            parser.error("\nmerge_bio_reps mode requires --bams-tech-rep tech_rep1.bam tech_rep2.bam")
        else:
            bams_bio_rep = [os.path.abspath(x) for x in bams_bio_rep]
    elif mode == "yahs":
        if not (bam or assembly_fasta):
            parser.error("\nmerge_tech_reps mode requires --bams-tech-rep tech_rep1.bam tech_rep2.bam")
        else:
            bam = os.path.abspath(bam)
            assembly_fasta = os.path.abspath(assembly_fasta)
        
    settings = {
        "assembly_fasta" : assembly_fasta,
        "fr" : forward_hic_read, 
        "rr" : reverse_hic_read,
        "bam" : bam,
        "bam_bio_rep" : bams_bio_rep,
        "bam_tech_rep" : bams_tech_rep,
        "outdir" : outdir,
        "threads" : threads,
        "mode" : mode,
        "execution_folder" : execution_folder,
        "debug" : debug,
        "config_file" : config_file,
        "braker_mode" : braker_file,
    }
    
    config_maker(settings, config_file)
    main(settings)