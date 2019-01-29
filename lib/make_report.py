"""
Makes SNP reports from Varsect output for each isolate that is suitable for use
in Gubbins to predict recombination regions.
"""

import argparse
from pysam import VariantFile

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-i', type=str, required=True, help="Input vcf")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-s', type=str, required=True, help="Name of sample")
    parser.add_argument('-r', type=str, help="Name of the reference of interest in the vcf for cases where plasmids were included")

    return parser.parse_args()

def main():
    args = parse_args()

    outfile = open("{}/{}_report.txt".format(args.o, args.s), "w")
    outfile.write("Sequence\tPosition in reference\tChange type\tOld\tNew\tEvidence\tConsequences\n")
    bcf_in = VariantFile(args.i)
    for rec in bcf_in.fetch():
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0]
        SVtype = rec.info["SVTYPE"]
        if SVtype == "SNP":
            if args.r:
                if chrom == args.r:
                    outfile.write("{}\t{}\tsubstitution\t{}\t{}\t.\n".format(chrom, pos, ref, alt))
            else:
                outfile.write("{}\t{}\tsubstitution\t{}\t{}\t.\n".format(chrom, pos, ref, alt))

if __name__ == '__main__':
    main()
