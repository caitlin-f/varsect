"""
Builds sample variant matrices based off previously intersected variant vcf
files for each individual isolate and populate using bcf files from individual
tools. Outputs the final common.vcf and unique calls in coregenome.vcf and
pangenome.vcf. Writes out the matrices as csv files. Also creates the necessary
alignment files (both phylip and nexus) for phylogenomic analysis in tools such
as RAxML and MrBayes.
"""

import argparse
import sys
import datetime
from pipeline_collate_classes import *


def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-r', type=str, required=True, help="Reference file in .fasta format")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-s', type=str, required=True, help="File of sample names")
    parser.add_argument('-E', type=str, help='Gubbins recombination gff file')
    parser.add_argument('-S', action='store_true', help="Include SNPs")
    parser.add_argument('-I', action='store_true', help="Include Indels")
    parser.add_argument('-D', action='store_true', help='Run Delly')
    parser.add_argument('-F', action='store_true', help='Run Freebayes')
    parser.add_argument('-G', action='store_true', help='Run GATK-HaplotypeCaller')
    parser.add_argument('-P', action='store_true', help='Run Pindel')

    return parser.parse_args()


class CombinedTrees:
    """ Build BST for each sample from single vcf (from intersect) """
    def __init__(self, args, sample_names):
        self.args = args
        self.sample_names = sample_names
        self.snp_trees = {}
        self.indel_trees = {}
        self.fail_snp_trees = {}
        self.fail_indel_trees = {}
        print("Building BST")

        self.buildAllTrees()


    def buildAllTrees(self):
        for sample in self.sample_names:
            print("Working on sample {} ({}/{})".format(sample, self.sample_names.index(sample)+1, len(self.sample_names)))
            # store 'failed' SV calls for each tool for each SV type (snp and indel/other)
            snp_fail_nodes = {}
            indel_fail_nodes = {}

            snp_nodes, indel_nodes, translocations = ReadBCF("Combined").read("{}/2_SVs/Final/{}.vcf".format(self.args.o, sample))
            if self.args.S:
                self.snp_trees[sample] = self.buildTrees(snp_nodes)
            if self.args.I:
                self.indel_trees[sample] = self.buildTrees(indel_nodes)

            if self.args.D and self.args.I: # Delly
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Delly").read("{0}/2_SVs/Delly/{1}_Results/{1}.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.F: # Freebayes
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Freebayes").read("{0}/2_SVs/Freebayes/{1}_Results/{1}.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.G: # gatk
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("GATK").read("{0}/2_SVs/GATK/{1}_Results/{1}.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.P and self.args.I: # Pindel
                tmp_indel_fail = {}
                tmp_end_fail = {}
                tmp_indel_fail, tmp_end_fail = ReadBCF("Pindel").read("{0}/2_SVs/Pindel/{1}_Results/{1}.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.S:
                self.fail_snp_trees[sample] = self.buildTrees(snp_fail_nodes)
            if self.args.I:
                self.fail_indel_trees[sample] = self.buildTrees(indel_fail_nodes)
                # remove false deletion signal in presence of translocation
                for trans_node in translocations:
                    try:
                        # perform twice for each Pindel and Delly failed nodes added (need to correct with proper range query)
                        trans_matches = self.fail_indel_trees[sample][trans_node.var.chrom].getNodes([], trans_node.var.end, 5, self.fail_indel_trees[sample][trans_node.var.chrom].root)
                        for match in trans_matches:
                            self.fail_indel_trees[sample][match.var.chrom].delete(match)
                        trans_matches = self.fail_indel_trees[sample][trans_node.var.chrom].getNodes([], trans_node.var.end, 5, self.fail_indel_trees[sample][trans_node.var.chrom].root)
                        for match in trans_matches:
                            self.fail_indel_trees[sample][match.var.chrom].delete(match)
                    except KeyError: # no SV calls on that chrom in failed call tree
                        pass


    def buildTrees(self, allNodes):
        """
        param:
            allNodes = dict{chrom:[nodes], chrom:[nodes], chrom:[nodes]}
        returns:
            dict{chom1:BST, chrom2:BST}
        """
        trees = {}
        for chrom in allNodes.keys():
            trees[chrom] = BinarySearchTree()
            trees[chrom].buildTree(allNodes[chrom])
        return trees

    def collateFailedNodes(self, all_fail_nodes, tmp_fail):
        """ Puts 'failed' call Nodes into dictionary for each chrom """
        for chrom in tmp_fail.keys():
            if all_fail_nodes.get(chrom) != None:
                all_fail_nodes[chrom].extend(tmp_fail[chrom])
            else:
                all_fail_nodes[chrom] = tmp_fail[chrom]
        return all_fail_nodes


class CombinedMatrix:
    def __init__(self, args, all_trees, samples, low_cov):
        self.args = args
        self.snp_trees = all_trees.snp_trees
        self.indel_trees = all_trees.indel_trees
        self.failed_snp = all_trees.fail_snp_trees
        self.failed_indel = all_trees.fail_indel_trees
        self.samples = samples
        self.low_cov = low_cov
        self.snp_positions = set()
        self.indel_positions = set()

        self.snp_mtx = None # coregenome snps
        self.indel_mtx = None # coregenome indels

        self.pan_snps = None
        self.pan_indels = None

        self.common_snps = None
        self.common_snps_low_cov = None

        self.common_indels = None
        self.common_indels_low_cov = None


        if self.args.S:
            self.buildSnpMtx()
            self.cleanSnpMtx()
        if self.args.I:
            self.buildIndelMtx()
            self.cleanIndelMtx()


    def buildSnpMtx(self):
        print("Building SNP Matrix")

        self.appendPositions(self.snp_positions, self.snp_trees)
        self.snp_positions = sorted(self.snp_positions)

        snp_mtx = np.zeros((len(self.samples), len(self.snp_positions)), dtype="int")

        self.__populateMtx(self.snp_positions, self.snp_trees, snp_mtx, 0)
        self.__populateMtx(self.snp_positions, self.failed_snp, snp_mtx, 0)
        self.snp_mtx = pd.DataFrame(snp_mtx, index=self.samples, columns=self.snp_positions)


    def buildIndelMtx(self):
        print("Building Indel Matrix")

        self.appendPositions(self.indel_positions, self.indel_trees)
        self.indel_positions = sorted(self.indel_positions)

        indel_mtx = np.zeros((len(self.samples), len(self.indel_positions)), dtype="int")

        self.__populateMtx(self.indel_positions, self.indel_trees, indel_mtx, 3)
        self.__populateMtx(self.indel_positions, self.failed_indel, indel_mtx, 3)

        self.indel_mtx = pd.DataFrame(indel_mtx, index=self.samples, columns=self.indel_positions)


    def __populateMtx(self, positions, trees, mtx, leniency):
        ("Populating Matrix")
        for r in range(len(self.samples)):
            sample = self.samples[r]
            print("Populating {} at {}".format(sample, datetime.datetime.time(datetime.datetime.now())))
            for c in range(len(positions)):
                var = positions[c]
                try:
                    if trees[sample][var.chrom].searchTree(var.pos, leniency, trees[sample][var.chrom].root):
                        mtx[r, c] = 1
                    elif "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and mtx[r, c] == 0:
                        mtx[r, c] = len(self.samples)*3
                except KeyError:
                    if "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and mtx[r, c] == 0:
                        mtx[r, c] = len(self.samples)*3


    def appendPositions(self, positions, tree):
        for sample in self.samples:
            for chrom in tree[sample].keys():
                node_list = tree[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)


    def cleanSnpMtx(self):
        print("Cleaning SNP matrix...")

        # remove cols where SNPs in all samples (shared SNPs)
        self.common_snps = self.snp_mtx[[col for col, val in self.snp_mtx.sum().iteritems() if val == len(self.samples)]]
        self.snp_mtx.drop([col for col, val in self.snp_mtx.sum().iteritems() if val == len(self.samples)], axis=1, inplace=True)
        # remove cols where either SNP call or low coverage
        self.common_snps_low_cov = self.snp_mtx[[col for col, val in self.snp_mtx.iteritems() if ((any((x == len(self.samples)*3) for x in val.values)) and len(set(val.values.tolist())) == 2)]]
        self.snp_mtx.drop([col for col, val in self.snp_mtx.iteritems() if ((any((x == len(self.samples)*3) for x in val.values)) \
            and len(set(val.values.tolist())) == 2)], axis=1, inplace=True)

        # remove cols where any sample has low coverage
        self.pan_snps = self.snp_mtx[[col for col, val in self.snp_mtx.sum().iteritems() if val > len(self.samples)]]
        self.snp_mtx.drop([col for col, val in self.snp_mtx.sum().iteritems() if val > len(self.samples)], axis=1, inplace=True)

        # reassign placeholder with N for no coverage
        self.snp_mtx[self.snp_mtx == len(self.samples)*3] = 'N'
        self.pan_snps[self.pan_snps == len(self.samples)*3] = 'N'
        self.common_snps[self.common_snps == len(self.samples)*3] = 'N'
        self.common_snps_low_cov[self.common_snps_low_cov == len(self.samples)*3] = 'N'

    def cleanIndelMtx(self):
        print("Cleaning Indel matrix...")

        # cols where INDELs in all samples (shared INDELs)
        self.common_indels = self.indel_mtx[[col for col, val in self.indel_mtx.sum().iteritems() if val == len(self.samples)]]

        # remove SVs appearing in all cols (either all samples have SV or low cov)
        col_sums = [] # store column sums where value is either 1 or '?' placeholder (i.e. len(sample_names)*3)
        for n in range(1, len(self.samples)+1):
            col_sums.append(len(self.samples)-n + (len(self.samples)*3)*n)
        # remove cols where sum of column is in the col_sums list
        self.common_indels_low_cov = self.indel_mtx[[col for col, val in self.indel_mtx.sum().iteritems() if val in col_sums]]
        self.indel_mtx.drop([col for col, val in self.indel_mtx.sum().iteritems() if (val == len(self.samples) or val in col_sums)], axis=1, inplace=True)

        # remove cols where any sample has low coverage (note can lose a bit of data taking just coregenome indels)
        self.pan_indels = self.indel_mtx[[col for col, val in self.indel_mtx.sum().iteritems() if val > len(self.samples)]]
        self.indel_mtx.drop([col for col, val in self.indel_mtx.sum().iteritems() if val > len(self.samples)], axis=1, inplace=True)

        # reassign placeholder with N for no coverage
        self.indel_mtx[self.indel_mtx == len(self.samples)*3] = 'N'
        self.pan_indels[self.pan_indels == len(self.samples)*3] = 'N'
        self.common_indels[self.common_indels == len(self.samples)*3] = 'N'
        self.common_indels_low_cov[self.common_indels_low_cov == len(self.samples)*3] = 'N'


def getIndex(alist, target):
    first = 0
    last = len(alist) -1
    found = -1
    while first <= last and found == -1:
        mid = (first + last) // 2
        if alist[mid] == target:
            found = mid
        else:
            if target < alist[mid]:
                last = mid - 1
            else:
                first = mid + 1
    return found



def main():

    args = parse_args()
    sample_names = get_sample_names(args)
    print("Building trees {}...".format(datetime.datetime.time(datetime.datetime.now())))
    all_trees = CombinedTrees(args, sample_names)
    print("Getting low coverage {}...".format(datetime.datetime.time(datetime.datetime.now())))
    low_cov = getLowCov(args, sample_names)

    # column headers are of class Variant
    print("Building matrices {}...".format(datetime.datetime.time(datetime.datetime.now())))

    mtx = CombinedMatrix(args, all_trees, sample_names, low_cov)

    # var_dif = mtx.getDifference()
    #
    # for sample in var_dif.keys():
    #     with open("{}/2_SVs/{}_var_dif.txt".format(args.o, sample), "w") as dif_out:
    #         dif_out.write("#CHROM\tSTART\tEND\tTYPE\tTOOL\tMTX\n")
    #         for line in sorted(var_dif[sample]):
    #             dif_out.write(line)


    if args.S: # core, pan and complete snps to csv
        mtx.snp_mtx.to_csv("{}/2_SVs/core_snps.csv".format(args.o))
        mtx.pan_snps.to_csv("{}/2_SVs/pan_snps.csv".format(args.o))
        complete_snps = pd.concat([mtx.snp_mtx, mtx.pan_snps], axis=1)
        complete_snps.to_csv("{}/2_SVs/complete_snps.csv".format(args.o))

    if args.I: # core, pan and complete indels to csv
        mtx.indel_mtx.to_csv("{}/2_SVs/core_indels.csv".format(args.o))
        mtx.pan_indels.to_csv("{}/2_SVs/pan_indels.csv".format(args.o))
        complete_indels = pd.concat([mtx.indel_mtx, mtx.pan_indels], axis=1)
        complete_indels.to_csv("{}/2_SVs/complete_indels.csv".format(args.o))

    if args.S and args.I: # concat snp and indel matrices
        common_mtx = pd.concat([mtx.common_snps, mtx.common_indels, mtx.common_snps_low_cov, mtx.common_indels_low_cov], axis=1)
        core = pd.concat([mtx.snp_mtx, mtx.indel_mtx], axis=1)
        pan = pd.concat([mtx.pan_snps, mtx.pan_indels], axis=1)
        complete = pd.concat([mtx.snp_mtx, mtx.indel_mtx, mtx.pan_snps, mtx.pan_indels], axis=1)
    elif args.S:
        common_mtx = pd.concat([mtx.common_snps, mtx.common_snps_low_cov], axis=1)
        core = mtx.snp_mtx
        pan = mtx.pan_snps
        complete = pd.concat([mtx.snp_mtx, mtx.pan_snps], axis=1)
    elif args.I:
        common_mtx = pd.concat([mtx.common_indels, mtx.common_indels_low_cov], axis=1)
        core = mtx.indel_mtx
        pan = mtx.pan_indels
        complete = pd.concat([mtx.indel_mtx, mtx.pan_indels], axis = 1)

    core.to_csv("{}/2_SVs/coregenome.csv".format(args.o))
    pan.to_csv("{}/2_SVs/pangenome.csv".format(args.o))
    complete.to_csv("{}/2_SVs/completegenome.csv".format(args.o))


    # common variants to vcf file (includes recombination regions)
    print("Writing common vcf {}...".format(datetime.datetime.time(datetime.datetime.now())))
    write_vcf(args, common_mtx, "common.vcf", sample_names, "tmp.vcf")


    if args.I:
        print("Writing unique indel vcf {}...".format(datetime.datetime.time(datetime.datetime.now())))
        write_vcf(args, mtx.indel_mtx, "core_indels.vcf", sample_names, "tmp.vcf")
        write_vcf(args, mtx.pan_indels, "pan_indels.vcf", sample_names, "tmp.vcf")

    if len(core.columns.values) < 100000:
        # write all variants to vcf (only if small number of snps ~ < 50k)
        print("Writing unique vcf {}...".format(datetime.datetime.time(datetime.datetime.now())))
        write_vcf(args, core, "coregenome.vcf", sample_names, "tmp.vcf")
        write_vcf(args, pan, "pangenome.vcf", sample_names, "tmp.vcf")

    # remove recombination regions before writing tree alignment files
    if args.E:
        os.system("cat {} | sort -k 4 -h > {}/2_SVs/recom_pred_sorted.gff".format(args.E, args.o))
        # for writing phylo files
        if args.S:
            remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), mtx.snp_mtx)
            remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), mtx.pan_snps)
        if args.I:
            remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), mtx.indel_mtx)
            remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), mtx.pan_indels)

        # for saving csv matrix minus recombination
        remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), core)
        core.to_csv("{}/2_SVs/coregenome_norec.csv".format(args.o))

        remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), pan)
        pan.to_csv("{}/2_SVs/pangenome_norec.csv".format(args.o))

        remove_recombination("{}/2_SVs/recom_pred_sorted.gff".format(args.o), complete)
        complete.to_csv("{}/2_SVs/completegenome_norec.csv".format(args.o))

    # if len(core.columns.values) < 100000:
        # write all variants to vcf (only if small number of snps ~ < 100k)
    print("Writing vcf without recombination {}...".format(datetime.datetime.time(datetime.datetime.now())))
    write_vcf(args, core, "coregenome_norec.vcf", sample_names, "tmp.vcf")
    write_vcf(args, pan, "pangenome_norec.vcf", sample_names, "tmp.vcf")

    print("Writing phylo files {}...".format(datetime.datetime.time(datetime.datetime.now())))
    write_phylo(args, sample_names, "core", snp_mtx = mtx.snp_mtx, indel_mtx = mtx.indel_mtx)

    if args.S:
        complete_snps = pd.concat([mtx.snp_mtx, mtx.pan_snps], axis=1)
    else:
        complete_snps = None

    if args.I:
        complete_indels = pd.concat([mtx.indel_mtx, mtx.pan_indels], axis = 1)
    else:
        complete_indels = None

    write_phylo(args, sample_names, "complete", snp_mtx = complete_snps, indel_mtx = complete_indels)


if __name__ == '__main__':
    main()
