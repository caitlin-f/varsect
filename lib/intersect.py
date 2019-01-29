"""
Find the intersection of variant calls from Freebayes and GATK HaplotypeCaller
(SNPs and small indels), and Delly and Pindel (large variants). Writes results
to a vcf file for each sample in Final/sample_name.vcf.

The output vcf file is used by matrix_fast.py
"""
import argparse
import sys
import datetime
from pipeline_collate_classes import *

"""
TODO: Make Final directory
"""

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-r', type=str, required=True, help="Reference file in .fasta format")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-s', type=str, required=True, help="Name of sample")
    parser.add_argument('-S', action='store_true', help="Include SNPs")
    parser.add_argument('-I', action='store_true', help="Include Indels")
    parser.add_argument('-D', action='store_true', help='Run Delly')
    parser.add_argument('-F', action='store_true', help='Run Freebayes')
    parser.add_argument('-G', action='store_true', help='Run GATK-HaplotypeCaller')
    parser.add_argument('-P', action='store_true', help='Run Pindel')

    return parser.parse_args()


def main():

    print("Finding variant insersection")
    args = parse_args()

    all_trees = AllTrees(args, [args.s])

    low_cov = getLowCov(args, [args.s])

    print("Building variant matrices...")

    # column headers are of class Variant
    mtx = VariantMatrix(args, all_trees.delly_trees, \
        all_trees.freebayes_snp_trees, all_trees.freebayes_indel_trees, \
        all_trees.gatk_snp_trees, all_trees.gatk_indel_trees, \
        all_trees.pindel_trees, \
        all_trees.fail_snp_trees, all_trees.fail_indel_trees, \
        [args.s], low_cov)

    # var_dif = mtx.getDifference()
    #
    # for sample in var_dif.keys():
    #     with open("{}/2_SVs/{}_var_dif.txt".format(args.o, sample), "w") as dif_out:
    #         dif_out.write("#CHROM\tSTART\tEND\tTYPE\tTOOL\tMTX\n")
    #         for line in sorted(var_dif[sample]):
    #             dif_out.write(line)

    #
    # if args.S:
    #     mtx.snp_mtx.to_csv("{}/2_SVs/unique_snps.csv".format(args.o))
    # if args.I:
    #     mtx.indel_mtx.to_csv("{}/2_SVs/unique_indels.csv".format(args.o))


    if args.S and args.I:
        common_mtx = pd.concat([mtx.common_snps, mtx.common_indels, mtx.common_snps_low_cov, mtx.common_indels_low_cov], axis=1)
        # unique_mtx = pd.concat([mtx.snp_mtx, mtx.indel_mtx], axis=1)
    elif args.S:
        common_mtx = pd.concat([mtx.common_snps, mtx.common_snps_low_cov], axis=1)
        # unique_mtx = mtx.snp_mtx
    elif args.I:
        common_mtx = pd.concat([mtx.common_indels, mtx.common_indels_low_cov], axis=1)
        # unique_mtx = mtx.indel_mtx

    # write_vcf(args, common_mtx, "Final/{}_common.vcf", sample_names)
    write_vcf(args, common_mtx, "Final/{}.vcf".format(args.s), [args.s], "tmp_{}.vcf".format(args.s))


    # write_phylo(args, mtx, sample_name)


if __name__ == '__main__':
    main()
