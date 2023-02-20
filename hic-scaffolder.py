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
    "replicates": "{settings["replicates"]}"
    "forward_hic_read": "{settings["fr"]}"
    "reverse_hic_read": "{settings["rr"]}"
    "mode" : "{settings["mode"]}"
    "bam": "{settings["bam"]}"
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
    snakemake   --snakefile {settings["execution_folder"]}/workflow/Snakefile.smk \
                --configfile {settings["config_file"]} \
                --cores {settings["threads"]} \
                --use-conda \
                --conda-frontend mamba \
                {snake_debug}"""
    print(command)
    os.system(command)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parrots pipeline - a tool for parrots genome annotation project.')
    parser.add_argument('-m','--mode', help="mode to use [default = align_hic_pair]", 
                        choices=["align_hic_pair", "align_hic_replicates", "yahs"], default="align_hic_pair")
    parser.add_argument('-a','--assembly', help="path to asssembly fasta file", default="", required=True)
    parser.add_argument('-1','--forward_hic_read', help="path to forward hic read", default="")
    parser.add_argument('-2','--reverse_hic_read', help="path to reverse hic read", default="")
    parser.add_argument('-r','--replicates', help="path to replicates CSV file", default="")
    parser.add_argument('-b', '--bam', help='path to bam pile to run YaHS', default='')
    parser.add_argument('-o','--outdir', help='output directory', required=True)
    parser.add_argument('-t','--threads', help='number of threads [default == 8]', default = "8")
    parser.add_argument('-d','--debug', help='debug mode', action='store_true')
    args = vars(parser.parse_args())

    mode = args["mode"]
    assembly_fasta = os.path.abspath(args["assembly"])
    threads = args["threads"]
    debug = args["debug"]
    forward_hic_read = args["forward_hic_read"]
    reverse_hic_read = args["reverse_hic_read"]
    replicates = args["replicates"]
    bam = args['bam']
    outdir = os.path.abspath(args["outdir"])
    
    execution_folder = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    execution_time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    config_random_letters = "".join([random.choice(string.ascii_letters) for n in range(3)])
    config_file = os.path.join(execution_folder, f"config/config_{config_random_letters}{execution_time}.yaml")
    # default braker_file for fasta mode
    braker_file = os.path.join(execution_folder,"rules/braker_fasta.smk")

    if mode == "align_hic_pair":
        if not (forward_hic_read or reverse_hic_read):
            parser.error("\nalign_hic_pair mode requires -1 {path_to_forward_read} and -2 {path_to_reverse_read}!")
        else:
            forward_hic_read = os.path.abspath(forward_hic_read)
            reverse_hic_read = os.path.abspath(reverse_hic_read)
    elif mode == "align_hic_replicates":
        if not (replicates):
            parser.error("\nalign_hic_pair mode requires -r {replicates.csv}!")
        else:
            replicates = os.path.abspath(replicates)
    elif mode == "yahs":
        if not (bam or assembly_fasta):
            parser.error("\nyahs mode requires --bam {sample.bam} -a {assembly.fasta}")
        else:
            bam = os.path.abspath(bam)
            assembly_fasta = os.path.abspath(assembly_fasta)
        
    settings = {
        "assembly_fasta" : assembly_fasta,
        "fr" : forward_hic_read, 
        "rr" : reverse_hic_read,
        "bam" : bam,
        "replicates": replicates,
        "outdir" : outdir,
        "threads" : threads,
        "mode" : mode,
        "execution_folder" : execution_folder,
        "debug" : debug,
        "config_file" : config_file,
    }
    
    config_maker(settings, config_file)
    main(settings)