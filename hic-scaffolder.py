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
import yaml

def config_maker(settings, config_file):
    
    if not os.path.exists(os.path.dirname(config_file)):
        os.mkdir(os.path.dirname(config_file))

    with open(config_file, "w") as f:
        yaml.dump(settings, f)
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

    parser = argparse.ArgumentParser(description='hic-scaffolder - a pipleine to align HiC reads to genome, scaffold it with YaHS adn get ready to curation HiC maps')
    parser.add_argument('-m','--mode', help="mode to use [default = assembly_to_hic_map]", 
                        choices=["assembly_to_hic_map", "reviewed_to_fasta"], default="assembly_to_hic_map")
    parser.add_argument('-s','--scaffolder', help="scaffolder to use [default = yahs]", 
                        choices=["yahs", "salsa"], default="yahs")
    parser.add_argument('-a','--assembly', help="path to asssembly fasta file", default="")
    parser.add_argument('-r','--reviewed', help="path to reviewed.asssembly file", default="")
    parser.add_argument('-f', '--run_folder', help="path to hic-scaffolder run folder to get reviewed assembly fasta", default="")
    parser.add_argument('-1','--forward_hic_read', help="path to forward hic read", default="")
    parser.add_argument('-2','--reverse_hic_read', help="path to reverse hic read", default="")
    parser.add_argument('-b', '--bam', help='path to bam pile to run YaHS', default='')
    parser.add_argument('-o','--outdir', help='output directory', required=True)
    parser.add_argument('-t','--threads', help='number of threads [default == 8]', default = "8")
    parser.add_argument('-d','--debug', help='debug mode', action='store_true')
    args = vars(parser.parse_args())

    mode = args["mode"]
    assembly_fasta = os.path.abspath(args["assembly"])
    threads = args["threads"]
    debug = args["debug"]
    forward_hic_read = os.path.abspath(args["forward_hic_read"])
    reverse_hic_read = os.path.abspath(args["reverse_hic_read"])
    bam = os.path.abspath(args['bam'])
    scaffolder = args['scaffolder']
    reviewed_assembly = os.path.abspath(args['reviewed'])
    run_folder = os.path.abspath(args['run_folder'])
    outdir = os.path.abspath(args["outdir"])
    
    execution_folder = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
    print(execution_folder)
    execution_time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    config_random_letters = "".join([random.choice(string.ascii_letters) for n in range(3)])
    config_file = os.path.join(execution_folder, f"config/config_{config_random_letters}{execution_time}.yaml")
    # default braker_file for fasta mode
    braker_file = os.path.join(execution_folder,"rules/braker_fasta.smk")
    os.chdir(execution_folder)

    if mode == "assembly_to_hic_map":
        if not ((forward_hic_read and reverse_hic_read) or bam):
            parser.error("\nassembly_to_hic_map mode requires -1 {path_to_forward_read} and -2 {path_to_reverse_read} or {path_to_bam}!")
        else:
            forward_hic_read = os.path.abspath(forward_hic_read)
            reverse_hic_read = os.path.abspath(reverse_hic_read)
            prefix = os.path.splitext(os.path.split(assembly_fasta)[1])[0]
    elif mode == "reviewed_to_fasta":
        if not (reviewed_assembly or run_folder):
            parser.error("\nreviewed_to_fasta mode requires -r {path_to_reviewed_assembly} and -f {path_to_hic-scaffolder_folder}!")
        else:
            reviewed_assembly = os.path.abspath(reviewed_assembly)
            prefix = os.path.splitext(os.path.split(reviewed_assembly)[1])[0]

        
    settings = {
        "scaffolder" : scaffolder,
        "assembly_fasta" : assembly_fasta,
        "reviewed_assembly" : reviewed_assembly,
        "run_folder": run_folder,
        "forward_hic_read" : forward_hic_read, 
        "reverse_hic_read" : reverse_hic_read,
        "bam" : bam,
        "outdir" : outdir,
        "threads" : threads,
        "mode" : mode,
        "execution_folder" : execution_folder,
        "debug" : debug,
        "config_file" : config_file,
        "prefix": prefix
    }
    
    config_maker(settings, config_file)
    main(settings)