"""
FILE: varsect.py
AUTHOR: Caitlin Falconer
USAGE: python3 varsect_batch.py -r reference.fa -o path/to/outdir -s samples.txt -t 16 [-E recomb.gff] [-A]
REQUIRES:
-r Reference fasta file and write access to reference sequence directory,
-o path to output directory
-s file containing a line separated list of paired fastq sample files in the
   order of sampleA_R1.fastq.gz, sampleA_R2.fastq.gz, sampleB_R1.fastq.gz,
   sampleB_R2.fastq.gz
-t number of threads to use

Run Varsect using Slurm batch scripts. Creates all the batch scripts for each of
the steps in Varsect. Creates the directories for the final output files and all
necessary indexes for the reference fasta file.
"""

import argparse
import os
from lib.write_batch_scripts import *

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-r', type=str, required=True, metavar='[reference.fa]', help="Reference file in .fasta format")
    parser.add_argument('-o', type=str, required=True, metavar='[path/to/outdir]', help="Output directory")
    parser.add_argument('-s', type=str, required=True, metavar='[samples.txt]', help="Read filenames to analyse")
    parser.add_argument('-t', type=int, required=True, metavar='[threads]', help='Number of threads')
    parser.add_argument('-E', type=str, metavar='[recomb.gff]', help="Remove recombination regions specified in a gff file")
    parser.add_argument('-A', action='store_true', help='Run all steps')
    parser.add_argument('-M', action='store_true', help='Perform Mapping')
    parser.add_argument('-D', action='store_true', help='Run Delly')
    parser.add_argument('-F', action='store_true', help='Run Freebayes')
    parser.add_argument('-G', action='store_true', help='Run GATK-HaplotypeCaller')
    parser.add_argument('-P', action='store_true', help='Run Pindel')
    parser.add_argument('-C', action='store_true', help='Collate data')
    parser.add_argument('-I', action='store_true', help='Make Samplot Images')
    parser.add_argument('-R', action='store_true', help='Run RaxML')
    parser.add_argument('-B', action='store_true', help='Run MrBayes')

    return parser.parse_args()

def read_sample_file(args):
    """ Extract sample names from input file """
    file_sets = {} # key = sample_name, value = list of paired files
    with open(args.s, "r") as input:
        for line in input:
            sample = line.strip().split("/")[-1].split("_")[0]
            if file_sets.get(sample) is None:
                file_sets[sample] = [line.strip()]
            else:
                file_sets[sample].append(line.strip())

    fwd_file =  open("{}/forward.reads".format(args.o), "w")
    rev_file = open("{}/reverse.reads".format(args.o), "w")
    names_file = open("{}/sample_names.txt".format(args.o), "w")
    fwd = ""
    rev = ""
    names = ""
    try:
        for sample in file_sets.keys():
            fwd += "{}\n".format(file_sets[sample][0])
            rev += "{}\n".format(file_sets[sample][1])
            names += "{}\n".format(sample)
        fwd_file.write(fwd)
        rev_file.write(rev)
        names_file.write(names)
    except IndexError:
        print("Error: Samples are not paired (specifically {}). Check input sample file".format(sample))
        exit(1)
    return file_sets

def make_directories(args):
    dir = "{}/batch_output/".format(args.o)
    if not os.path.exists(os.path.dirname(dir)):
        os.makedirs(dir)
    dir = "{}/1_Mapping/".format(args.o)
    if not os.path.exists(os.path.dirname(dir)):
        os.mkdir(dir)
    dir = "{}/2_SVs/Final".format(args.o)
    if not os.path.exists(os.path.dirname(dir)):
        os.makedirs(dir)
    dir = "{}/3_Trees/".format(args.o)
    if not os.path.exists(os.path.dirname(dir)):
        os.mkdir(dir)
    dir = "{}/4_Images/".format(args.o)
    if not os.path.exists(os.path.dirname(dir)):
        os.mkdir(dir)
    if args.D:
        dir = "{}/2_SVs/Delly/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.makedirs(dir)
    if args.F:
        dir = "{}/2_SVs/Freebayes/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.makedirs(dir)
    if args.G:
        dir = "{}/2_SVs/GATK/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.makedirs(dir)
    if args.P:
        dir = "{}/2_SVs/Pindel/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.makedirs(dir)
    if args.R:
        dir = "{}/3_Trees/RaxML/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.mkdir(dir)
    if args.B:
        dir = "{}/3_Trees/MrBayes/".format(args.o)
        if not os.path.exists(os.path.dirname(dir)):
            os.mkdir(dir)

def index_reference(args):
    """ Check that reference is properly indexed
    Requires:
        samtools faidx (fasta.fai)
        bwa index (fasta.amb, fasta.ann, fasta.bwt, fasta.pac, fasta.sa
        gatk CreateSequenceDictionary (.dict))"""

    dir = args.r + '.fai'
    if not os.path.exists(dir):
        print("Building samtools index")
        os.system("samtools faidx {}".format(args.r))

    index = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    for i in index:
        dir = args.r + i
        if not os.path.exists(dir):
            print("Building bwa index")
            os.system("bwa index {}".format(args.r))

    dir = args.r.split(".fa")[0] + '.dict'
    if not os.path.exists(dir):
        print("Building gatk index")
        os.system("gatk-launch CreateSequenceDictionary -R {}".format(args.r))


def main():
    args = parse_args()
    if args.A:
        args.M = True
        args.D = True
        args.F = True
        args.G = True
        args.P = True
        args.C = True
        args.I = True
        args.R = True
        args.B = True

    make_directories(args)
    file_sets = read_sample_file(args) # get sample names and read file names
    write_batch_files(args, file_sets)

    index_reference(args) # make sure reference is indexed properly


if __name__ == '__main__':
    main()
