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

        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building BST {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building BST")

        self.buildAllTrees()


    def buildAllTrees(self):
        for sample in self.sample_names:
            f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
            f.write("Working on sample {} ({}/{}) {}\n".format(sample, self.sample_names.index(sample)+1, len(self.sample_names), datetime.datetime.time(datetime.datetime.now())))
            print("Working on sample {} ({}/{})".format(sample, self.sample_names.index(sample)+1, len(self.sample_names)))
            # store 'failed' SV calls for each tool for each SV type (snp and indel/other)
            snp_fail_nodes = {}
            indel_fail_nodes = {}

            snp_nodes, indel_nodes = ReadBCF("Combined").read("{}/2_SVs/Final/{}.vcf".format(self.args.o, sample))
            if self.args.S:
                self.snp_trees[sample] = self.buildTrees(snp_nodes)
            if self.args.I:
                self.indel_trees[sample] = self.buildTrees(indel_nodes)

            if self.args.D and self.args.I: # Delly
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Delly").read("{0}/2_SVs/Delly/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.F: # Freebayes
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Freebayes").read("{0}/2_SVs/Freebayes/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.G: # gatk
                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("GATK").read("{0}/2_SVs/GATK/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.P and self.args.I: # Pindel
                tmp_indel_fail = {}
                tmp_end_fail = {}
                tmp_indel_fail, tmp_end_fail = ReadBCF("Pindel").read("{0}/2_SVs/Pindel/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.S:
                self.fail_snp_trees[sample] = self.buildTrees(snp_fail_nodes)
            if self.args.I:
                self.fail_indel_trees[sample] = self.buildTrees(indel_fail_nodes)

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
        self.snp_mtx = None
        self.indel_mtx = None
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
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building SNP Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building SNP Matrix")

        self.appendPositions(self.snp_positions, self.snp_trees)

        self.snp_mtx = pd.DataFrame(np.zeros((len(self.samples), len(self.snp_positions))), index=self.samples, columns=sorted(self.snp_positions), dtype="int")

        self.__populateMtx(self.snp_positions, self.snp_trees, self.snp_mtx, 0)
        self.__populateMtx(self.snp_positions, self.failed_snp, self.snp_mtx, 0)


    def buildIndelMtx(self):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building Indel Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building Indel Matrix")

        self.appendPositions(self.indel_positions, self.indel_trees)

        self.indel_mtx = pd.DataFrame(np.zeros((len(self.samples), len(self.indel_positions))), index=self.samples, columns=sorted(self.indel_positions), dtype="int")

        self.__populateMtx(self.indel_positions, self.indel_trees, self.indel_mtx, 3)
        self.__populateMtx(self.indel_positions, self.failed_indel, self.indel_mtx, 3)


    def __populateMtx(self, positions, trees, mtx, leniency):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Populating Matrix")
        for sample in self.samples:
            f.write("Populating {} at {}\n".format(sample, datetime.datetime.time(datetime.datetime.now())))
            print("Populating {} at {}".format(sample, datetime.datetime.time(datetime.datetime.now())))
            for var in positions:
                try:
                    if trees[sample][var.chrom].searchTree(var.pos, leniency, trees[sample][var.chrom].root):
                        mtx.loc[sample, var] = 1
                    elif "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and mtx.loc[sample, var] == 0:
                        mtx.loc[sample, var] = len(self.samples)*3
                except KeyError:
                    if "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and mtx.loc[sample, var] == 0:
                        mtx.loc[sample, var] = len(self.samples)*3


    def appendPositions(self, positions, tree):
        for sample in self.samples:
            for chrom in tree[sample].keys():
                node_list = tree[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)

    def cleanSnpMtx(self):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Cleaning SNP Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))

        print("Cleaning SNP matrix...")

        # remove cols where SNPs in all samples (shared SNPs)
        self.common_snps = self.snp_mtx[[col for col, val in self.snp_mtx.sum().iteritems() if val == len(self.samples)]]
        self.snp_mtx.drop([col for col, val in self.snp_mtx.sum().iteritems() if val == len(self.samples)], axis=1, inplace=True)
        # remove cols where either SNP call or '?'
        self.common_snps_low_cov = self.snp_mtx[[col for col, val in self.snp_mtx.iteritems() if ((any((x == len(self.samples)*3) for x in val.values)) and len(set(val.values.tolist())) == 2)]]
        self.snp_mtx.drop([col for col, val in self.snp_mtx.iteritems() if ((any((x == len(self.samples)*3) for x in val.values)) \
            and len(set(val.values.tolist())) == 2)], axis=1, inplace=True)

        self.snp_mtx[self.snp_mtx == len(self.samples)*3] = 'N'
        self.common_snps[self.common_snps == len(self.samples)*3] = 'N'
        self.common_snps_low_cov[self.common_snps_low_cov == len(self.samples)*3] = 'N'

    def cleanIndelMtx(self):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Cleaning Indel Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))

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

        self.indel_mtx[self.indel_mtx == len(self.samples)*3] = 'N'
        self.common_indels[self.common_indels == len(self.samples)*3] = 'N'
        self.common_indels_low_cov[self.common_indels_low_cov == len(self.samples)*3] = 'N'
        # # drop columns where ANY sample has low cov (generally lose unique SV data doing this)
        # indel_mtx.drop([col for col, val in indel_mtx.sum().iteritems() if val > len(sample_names)], axis=1, inplace=True)






def main():


    args = parse_args()
    f = open("{}/2_SVs/collate_status.txt".format(args.o), 'a')
    f.write("Building sample variant matrices {}\n".format(datetime.datetime.time(datetime.datetime.now())))

    print("Building variant matrices...")

    sample_names = get_sample_names(args)
    all_trees = CombinedTrees(args, sample_names)
    low_cov = getLowCov(args, sample_names)

    f.write("Instantiating CombinedMatrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
    f.close()

    # column headers are of class Variant
    mtx = CombinedMatrix(args, all_trees, sample_names, low_cov)

    # var_dif = mtx.getDifference()
    #
    # for sample in var_dif.keys():
    #     with open("{}/2_SVs/{}_var_dif.txt".format(args.o, sample), "w") as dif_out:
    #         dif_out.write("#CHROM\tSTART\tEND\tTYPE\tTOOL\tMTX\n")
    #         for line in sorted(var_dif[sample]):
    #             dif_out.write(line)

    #
    if args.S:
        mtx.snp_mtx.to_csv("{}/2_SVs/unique_snps.csv".format(args.o))
    if args.I:
        mtx.indel_mtx.to_csv("{}/2_SVs/unique_indels.csv".format(args.o))

    if args.S and args.I:
        common_mtx = pd.concat([mtx.common_snps, mtx.common_indels, mtx.common_snps_low_cov, mtx.common_indels_low_cov], axis=1)
        unique_mtx = pd.concat([mtx.snp_mtx, mtx.indel_mtx], axis=1)
        unique_mtx.to_csv("{}/2_SVs/unique.csv".format(args.o))
    elif args.S:
        common_mtx = pd.concat([mtx.common_snps, mtx.common_snps_low_cov], axis=1)
        unique_mtx = mtx.snp_mtx
    elif args.I:
        common_mtx = pd.concat([mtx.common_indels, mtx.common_indels_low_cov], axis=1)
        unique_mtx = mtx.indel_mtx

    write_vcf(args, common_mtx, "common.vcf", sample_names)
    write_vcf(args, unique_mtx, "unique.vcf", sample_names)

    write_phylo(args, mtx, sample_names)


if __name__ == '__main__':
    main()
