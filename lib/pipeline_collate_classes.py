"""
Miscellaneous classes used during the collation stage of Varsect
"""

import csv
import random
import sys
import os
import gc
import gff_parser
import re
import numpy as np
import pandas as pd
from pysam import VariantFile

import datetime

csv.field_size_limit(sys.maxsize)

LENIENCY_SNP = 0
LENIENCY_SMALL = 3 # leniency to compare positions of SVs
LENIENCY_LARGE = 10 # leniency when writing which SVs found to file

class Node:
    """ Node of the BST
    params:
       pos (int): variant start position
       var (Variant): Variant class for that variant
       left (Node): left child in BST
       right (Node): right child in BST
       parent (Node): parent in BST """
    def __init__(self, pos, var, left=None, right=None, parent=None):
        self.pos = int(pos)
        self.var = var # Class Variant
        self.leftChild = left
        self.rightChild = right
        self.parent = parent

    def isRoot(self):
        return not self.parent # not None

    def isLeaf(self):
        return not (self.leftChild or self.rightChild)

    def getAllChildren(self, childrenList):
        """ returns [Node, Node, ...] """
        childrenList.append(self)
        if self.leftChild is not None:
            self.leftChild.getAllChildren(childrenList)
        if self.rightChild is not None:
            self.rightChild.getAllChildren(childrenList)
        return childrenList

    def getAllPositions(self, positionList):
        """ returns ["chrom_pos", "chrom_pos", ...]"""
        positionList.append("{}:{}".format(self.var.chrom, self.pos))
        if self.leftChild is not None:
            self.leftChild.getAllPositions(positionList)
        if self.rightChild is not None:
            self.rightChild.getAllPositions(positionList)
        return positionList

    def __str__(self):
        return "TYPE:{}, START:{}, END:{}, LENGTH:{}".format(self.var.type,self.pos,self.var.end,self.var.length)

    def __repr__(self):
        return self.__str__()


class BinarySearchTree:
    """ BST data stucture
    params:
        root (Node): root of BST """
    def __init__(self, root=None):
        self.root = root

    def __setRoot(self, root):
        self.root = root # Class Node

    def insert(self, node):
        if self.root is None:
            self.__setRoot(node)
        else:
            self.__insertNode(self.root, node)

    def __insertNode(self, currentNode, newNode):
        if currentNode != newNode:
            if currentNode.pos >= newNode.pos: # new pos <= current pos
                if currentNode.leftChild != None:
                    self.__insertNode(currentNode.leftChild, newNode)
                else:
                    currentNode.leftChild = newNode
                    newNode.parent = currentNode
            else: # new pos > current pos
                if currentNode.rightChild != None:
                    self.__insertNode(currentNode.rightChild, newNode)
                else:
                    currentNode.rightChild = newNode
                    newNode.parent = currentNode

    def delete(self, node):
        if node.isLeaf(): # node has no children
            # remove parent link to node
            if node.parent.leftChild == node:
                node.parent.leftChild = None
            elif node.parent.rightChild == node:
                node.parent.rightChild = None

        elif node.rightChild == None: # no successor
            # link left child to parent
            if node.parent != None: # node is not root
                if node.parent.leftChild == node:
                    node.parent.leftChild = node.leftChild
                if node.parent.rightChild == node:
                    node.parent.rightChild = node.leftChild
                node.leftChild.parent = node.parent
            else: # node is root
                node.leftChild.parent = None
                self.root = node.leftChild

        else:
            successor = self.find_min(node.rightChild)

            # 1. node.parent.child - successor
            if node.parent != None: # node is not root
                if node.parent.leftChild == node:
                    node.parent.leftChild = successor
                elif node.parent.rightChild == node:
                    node.parent.rightChild = successor

            # 2. successor.parent.child = successor.child
            if successor.parent.leftChild == successor:
                successor.parent.leftChild = successor.rightChild
            elif successor.parent.rightChild == successor: # occurs if successor parent is node being deleted
                successor.parent.rightChild = successor.rightChild

            # 3. successor.rightChild = successor.parent
            if successor.rightChild != None: # if successor has any children (should only have right child)
                successor.rightChild.parent = successor.parent

            # 4. successor.parent = node.parent
            successor.parent = node.parent
            if node.parent == None: # node is root
                self.root = successor

            # 5. successor.children = node.children
            successor.leftChild = node.leftChild
            if node.leftChild != None:
                node.leftChild.parent = successor
            successor.rightChild = node.rightChild
            if node.rightChild != None:
                node.rightChild.parent = successor


    def find_min(self, currentNode):
        if currentNode.leftChild == None:
            return currentNode
        else:
            return self.find_min(currentNode.leftChild)

    def find_max(self, currentNode):
        if currentNode.rightChild is None:
            return currentNode
        else:
            return self.find_max(currentNode.rightChild)

    def searchTree(self, pos, leniency, currentNode):
        if currentNode is None: # no more nodes
            return False
        elif pos in range(currentNode.pos-leniency, currentNode.pos+leniency+1): # found match in range
            return True
        elif pos <= currentNode.pos-leniency: # pos smaller, search left
            return self.searchTree(pos, leniency, currentNode.leftChild)
        else: # pos greater, search right
            return self.searchTree(pos, leniency, currentNode.rightChild)

    def __str__(self):
        children = []
        return(str(self.root.getAllChildren(children)))

    def __repr__(self):
        return self.__str__()


    def getHeight(self, node):
        if node is None:
            return 0
        else:
            leftDepth = self.getHeight(node.leftChild)
            rightDepth = self.getHeight(node.rightChild)

        if (leftDepth > rightDepth):
            return leftDepth + 1
        else:
            return rightDepth + 1

    def buildTree(self, nodes):
        """
        param nodes = [node, node, ...]
        """
        random.shuffle(nodes) # randomise list to build more balanced BST
        for node in nodes:
            self.insert(node)

    def getAllPositions(self):
        """ returns ["chrom_pos", "chrom_pos", ...]"""
        positions = []
        return self.root.getAllPositions(positions)

    def getAllNodes(self):
        nodes = []
        return self.root.getAllChildren(nodes)

    def getNodes(self, nodeList, pos, leniency, currentNode):
        """ Return all nodes that match within a range of a position -
        Not returning nodes properly if both left and right child fall within
        the range. TODO: change search to a visitor pattern """
        if currentNode is None:
            return nodeList
        elif pos in range(currentNode.pos-leniency, currentNode.pos+leniency+1):
            nodeList.append(currentNode)
        if pos <= currentNode.pos:
            return self.getNodes(nodeList, pos, leniency, currentNode.leftChild)
        else:
            return self.getNodes(nodeList, pos, leniency, currentNode.rightChild)


class Variant:
    """ Variant class storing details of a single variant called
    params:
        chrom (str): chromosome/reference name variant was detected in
        pos (int): start position
        end (int): end position
        type (str): variant type e.g. SNP, INDEL, DEL, TSNL
        tool (str): variant calling tool used
        ref (str): reference allele (if SNP or small indel, otherwise '*')
        alt (str): alternate allele (if SNP or small indel, otherwise '*') """
    def __init__(self, chrom, pos, end, SVtype, length, tool, ref, alt):
        self.chrom = chrom
        self.pos = int(pos)
        self.end = int(end)
        self.type = SVtype
        self.length = int(length)
        self.tool = tool
        self.ref = ref
        self.alt = alt


    def __str__(self):
        return("{}:{}".format(self.chrom, self.pos))

    def __repr__(self):
        return("{}:{}".format(self.chrom, self.pos))

    def __eq__(self, other):
        if self.chrom == other.chrom:
            return self.pos == other.pos
        else:
            return False

    def __ne__(self, other):
        if self.chrom == other.chrom:
            return self.pos != other.pos
        else:
            return True

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.pos < other.pos
        else:
            return self.chrom < other.chrom

    def __le__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__lt__(self, other)

    def __gt__(self, other):
        if self.chrom == other.chrom:
            return self.pos > other.pos
        else:
            return self.chrom > other.chrom

    def __ge__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__gt__(self, other)

    def __hash__(self):
        chrom_hash = hash(self.chrom)
        pos_hash = hash(self.pos)
        return int("{}{}".format(chrom_hash, pos_hash))


class ReadVCF:
    """ Read in VCF files from the individual variant calling tools
    params:
        tool (str): variant calling tool used to produce vcf file """
    def __init__(self, tool):
        self.tool = tool
        self.SNP_nodes = {}
        self.INDEL_nodes = {}
        self.END_nodes = {} # for building tree of END positions for Pindel

    def read(self, file):
        if self.tool == "Delly":
            return self.__dellyReader(file)
        elif self.tool == "Freebayes":
            return self.__fbReader(file)
        elif self.tool == "GATK":
            return self.__gatkReader(file)
        elif self.tool == "Pindel":
            return self.__pindelReader(file)

    def appendNode(self, node, var_type):
        """
        Params:
          node: Node
          type: which dictionary to append to (SNP or INDEL)
        """
        chrom = node.var.chrom
        if var_type == "SNP":
            if self.SNP_nodes.get(chrom) != None:
                self.SNP_nodes.get(chrom).append(node)
            else:
                self.SNP_nodes[chrom] = [node]
        elif var_type == "INDEL":
            if self.INDEL_nodes.get(chrom) != None:
                self.INDEL_nodes.get(chrom).append(node)
            else:
                self.INDEL_nodes[chrom] = [node]
        elif var_type == "END":
            if self.END_nodes.get(chrom) != None:
                self.END_nodes.get(chrom).append(node)
            else:
                self.END_nodes[chrom] = [node]

    def __dellyReader(self, file):
        try:
            with open(file) as tsv:
                for line in csv.reader(tsv, delimiter="\t"):
                    if line[0].startswith('#'):
                        pass
                    else:
                        chrom = line[0]
                        pos = int(line[1])
                        for info in line[7].split(";"):
                            if info.startswith('SVTYPE'):
                                SVtype = info.split("=")[1]
                            if info.startswith('END'):
                                end = int(info.split("=")[1])
                                length = pos-end

                        var = Variant(chrom, pos, end, SVtype, length, self.tool, 0, 1)
                        node = Node(pos, var)
                        self.appendNode(node, "INDEL")
            tsv.close()
        except FileNotFoundError:
          print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes)

    def __fbReader(self, file):
        try:
            with open(file) as tsv:
                for line in csv.reader(tsv, delimiter="\t"):
                    if line[0].startswith('#'):
                        pass
                    else:
                        chrom = line[0]
                        pos = int(line[1])
                        for info in line[7].split(";"):
                            if info.startswith("TYPE"):
                                SVtype = info.split("=")[1].upper()
                            if info.startswith("LEN"):
                                length = int(info.split("=")[1].split(",")[0])
                        if SVtype == "SNP":
                            ref = line[3]
                            alt = line[4]
                            var_type = "SNP"
                            length = 0 # freebayes vcf reports SNP length as 1
                            var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                            node = Node(pos, var)
                            self.appendNode(node, var_type)
                        elif SVtype == "COMPLEX":
                            ref = line[3]
                            alt = line[4]
                            length = len(alt)
                            # Put each position in range of leniency in INDEL matrix
                            for n in range(0, length, LENIENCY_SMALL):
                                var = Variant(chrom, pos+n, pos+n+LENIENCY_SMALL, SVtype, length, self.tool, 0, 1)
                                node = Node(pos+n, var)
                                self.appendNode(node, "INDEL")
                            # Also put each position of complex in SNP matrix
                            for n in range(0,length):
                                var = Variant(chrom, pos+n, pos+n, SVtype, length, self.tool, ref, alt)
                                node = Node(pos+n, var)
                                self.appendNode(node, "SNP")
                        elif SVtype == "MNP":
                            ref = line[3]
                            alt = line[4]
                            length = len(alt)
                            for n in range(0, length):
                                if ref[n] != alt[n]: # skip pos in MNP that are not variants
                                    var = Variant(chrom, pos+n, pos+n, SVtype, length, self.tool, ref[n], alt[n])
                                    node = Node(pos+n, var)
                                    self.appendNode(node, "SNP")
                        else:
                            ref = 0
                            alt = 1
                            var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                            node = Node(pos, var)
                            self.appendNode(node, "INDEL")
            tsv.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes)

    def __gatkReader(self, file):
        try:
            with open(file) as tsv:
                for line in csv.reader(tsv, delimiter="\t"):
                    if line[0].startswith('#'):
                        pass
                    else:
                        chrom = line[0]
                        pos = int(line[1])
                        length = len(line[4]) - len(line[3])
                        if length == 0:
                            SVtype = "SNP"
                            ref = line[3]
                            alt = line[4]
                            var_type = "SNP"
                        else:
                            SVtype = "INDEL"
                            ref = 0
                            alt = 1
                            var_type = "INDEL"

                        var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                        node = Node(pos, var)
                        self.appendNode(node, var_type)
            tsv.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes)

    def __pindelReader(self, file):
        try:
            with open(file) as tsv:
                for line in csv.reader(tsv, delimiter="\t"):
                    if line[0].startswith('#'):
                        pass
                    else:
                        # if file.endswith("SI.vcf") or file.endswith("D.vcf"):
                        #   chrom, pos, SVtype, length = self.__pindel_D_SI(line)

                        # elif file.endswith("LI.vcf"):
                        #   chrom, pos, SVtype, length = self.__pindel_LI(line)
                        chrom, pos, end, SVtype, length = self.__pindel_D_SI(line)
                        if SVtype != "RPL": # ignore 'replacements from Pindel'
                            var = Variant(chrom, pos, end, SVtype, length, self.tool, 0, 1)
                            node = Node(pos, var)
                            self.appendNode(node, "INDEL")
                            end_node = Node(var.end, var)
                            self.appendNode(end_node, "END")
            tsv.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.INDEL_nodes, self.END_nodes)

    def __pindel_D_SI(self, line):
        chrom = line[0]
        pos = int(line[1])
        for info in line[7].split(";"):
            if info.startswith("SVTYPE"):
                SVtype=info.split("=")[1]
            if info.startswith("SVLEN"):
                length = info.split("=")[1]
            if info.startswith("END"):
                end = info.split("=")[1]
        return chrom, pos, end, SVtype, length

    def __pindel_LI(self, line):
        ref = line[0]
        chrom = line[0]
        pos = int(line[1])
        AD = line[9].split(':')[1]
        mpileup = subprocess.check_output(['grep','{}\t{}\t'.format(ref,pos),'{}/../1_Mapping/{}.mpileup'.format(DIR,FILE)])
        total_cov = 0
        for result in str(mpileup).lstrip('b\'').split('\\n'):
            try:
                cov=result.split('\\t')[3]
                total_cov+=int(cov)
            except IndexError:
                pass
            if total_cov - int(AD) < 0: # if pindel allele depth > coverage detected by samtools
                return chrom, pos, SVtype, length


class ReadBCF:
    """ Read in BCF files from the individual variant calling tools or the
    Final intersected results
    params:
        tool (str): variant calling tool used to produce vcf file """
    def __init__(self, tool):
        self.tool = tool
        self.SNP_nodes = {}
        self.INDEL_nodes = {}
        self.END_nodes = {} # for building tree of END positions for Pindel

    def read(self, file):
        if self.tool == "Delly":
            return self.__dellyReader(file)
        elif self.tool == "Freebayes":
            return self.__fbReader(file)
        elif self.tool == "GATK":
            return self.__gatkReader(file)
        elif self.tool == "Pindel":
            return self.__pindelReader(file)
        elif self.tool == "Combined":
            return self.__combinedReader(file)

    def appendNode(self, node, var_type):
        """
        Params:
          node: Node
          type: which dictionary to append to (SNP or INDEL)
        """
        chrom = node.var.chrom
        if var_type == "SNP":
            if self.SNP_nodes.get(chrom) != None:
                self.SNP_nodes.get(chrom).append(node)
            else:
                self.SNP_nodes[chrom] = [node]
        elif var_type == "INDEL":
            if self.INDEL_nodes.get(chrom) != None:
                self.INDEL_nodes.get(chrom).append(node)
            else:
                self.INDEL_nodes[chrom] = [node]
        elif var_type == "END":
            if self.END_nodes.get(chrom) != None:
                self.END_nodes.get(chrom).append(node)
            else:
                self.END_nodes[chrom] = [node]

    def __dellyReader(self, file):
        try:
            bcf_in = VariantFile(file)
            for rec in bcf_in.fetch():
                chrom = rec.chrom
                pos = rec.pos
                SVtype = rec.info["SVTYPE"]
                end = rec.stop
                length = end-pos
                var = Variant(chrom, pos, end, SVtype, length, self.tool, 0, 1)
                node = Node(pos, var)
                self.appendNode(node, "INDEL")
            bcf_in.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return(self.SNP_nodes, self.INDEL_nodes)

    def __fbReader(self, file):
        try:
            bcf_in = VariantFile(file)
            for rec in bcf_in.fetch():
                chrom = rec.chrom
                pos = rec.pos
                SVtype = rec.info["TYPE"][0].upper()
                ref = rec.ref
                alt = rec.alts[0]
                length = rec.info["LEN"][0]
                if SVtype == "SNP":
                    for n in range(0, len(alt)): # some SNPs recorded with multiple positions
                        if ref[n] != alt[n]: # skip pos in the SNP that are not variants
                            var = Variant(chrom, pos+n, pos+n, SVtype, 0, self.tool, ref[n], alt[n])
                            node = Node(pos+n, var)
                            self.appendNode(node, "SNP")
                elif SVtype == "COMPLEX":
                    # complex variants can include both snps and indels
                    for n in range(0, length, LENIENCY_SMALL):
                        var = Variant(chrom, pos+n, pos+n+LENIENCY_SMALL, SVtype, length, self.tool, ref[n:], alt[n:]) # repr each complex pos into a small indel
                        node = Node(pos+n, var)
                        self.appendNode(node, "INDEL")
                    for n in range(0, length):
                        var = Variant(chrom, pos+n, pos+n, SVtype, 0, self.tool, ref[n:], alt[n:])
                        node = Node(pos+n, var)
                        self.appendNode(node, "SNP")
                elif SVtype == "MNP":
                    for n in range(0, length):
                        if ref[n] != alt[n]: # skip pos in MNP that are not variants
                            var = Variant(chrom, pos+n, pos+n, SVtype, 0, self.tool, ref[n], alt[n])
                            node = Node(pos+n, var)
                            self.appendNode(node, "SNP")
                else: #INDEL
                    length = rec.info["LEN"][0]
                    var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                    node = Node(pos, var)
                    self.appendNode(node, "INDEL")
            bcf_in.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes)

    def __gatkReader(self, file):
        try:
            bcf_in = VariantFile(file)
            for rec in bcf_in.fetch():
                chrom = rec.chrom
                pos = rec.pos
                ref = rec.ref
                alt = rec.alts[0]
                length = len(alt)-len(ref)
                if length == 0:
                    SVtype = "SNP"
                    var_type = "SNP"
                else:
                    SVtype = "INDEL"
                    var_type = "INDEL"
                var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                node = Node(pos, var)
                self.appendNode(node, var_type)
            bcf_in.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes)

    def __pindelReader(self, file):
        try:
            bcf_in = VariantFile(file)
            for rec in bcf_in.fetch():
                chrom = rec.chrom
                pos = rec.pos
                end = rec.stop
                SVtype = rec.info["SVTYPE"]
                length = rec.info["SVLEN"]
                if SVtype != "RPL": # ignore 'replacements'
                    # end position can be incorrectly reported in some cases by pysam
                    if end + length != pos:
                        end = pos + abs(length)
                    if length not in range(-20,20): # ignore small indels
                        var = Variant(chrom, pos, end, SVtype, length, self.tool, 0, 1)
                        node = Node(pos, var)
                        self.appendNode(node, "INDEL")
                        end_node = Node(var.end, var)
                        self.appendNode(end_node, "END")
            bcf_in.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.INDEL_nodes, self.END_nodes)

    def __combinedReader(self, file):
        translocations = [] # need to return translocations to remove false deletion from BST
        try:
            bcf_in = VariantFile(file)
            for rec in bcf_in.fetch():
                chrom = rec.chrom
                pos = rec.pos
                end = rec.stop
                ref = rec.ref
                alt = rec.alts[0]
                SVtype = rec.info["SVTYPE"]
                length = rec.info["LEN"]
                if SVtype == "SNP":
                    var_type = "SNP"
                else:
                    var_type = "INDEL"
                var = Variant(chrom, pos, pos+length, SVtype, length, self.tool, ref, alt)
                node = Node(pos, var)
                self.appendNode(node, var_type)
                if SVtype == "TSLN":
                    translocations.append(node)

            bcf_in.close()
        except FileNotFoundError:
            print("Error: {} file not found".format(file))
        return (self.SNP_nodes, self.INDEL_nodes, translocations)

class MtxColumns:

    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def __str__():
        return("{}:{}".format(chrom, pos))

    def __repr__():
        return("{}:{}".format(chrom, pos))


class VariantMatrix:
    """ Sample x Variant matrix. Store data as a binary 0 = Variant not detected
    1 = Variant detected. """

    def __init__(self, args, delly, fb_snp, fb_indel, gatk_snp, gatk_indel, \
        pindel, failed_snp, failed_indel, samples, low_cov):
        """
        params:
            delly, freebayes, gatk, pindel = [boolean, {sampleA:{chrom1:BST, chrom2:BST}, sampleB:{chrom1:BST, chrom2:BST}}]
            samples = ["sampleA", "sampleB"]
            low_cov = {"sampleA":[chrom_pos, chrom_pos], "sampleB":[...]}
        """
        self.args = args
        self.d_trees = delly
        self.f_snp_trees = fb_snp
        self.f_indel_trees = fb_indel
        self.g_snp_trees = gatk_snp
        self.g_indel_trees = gatk_indel
        self.p_trees = pindel
        self.fail_snp_trees = failed_snp
        self.fail_indel_trees = failed_indel
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
        if self.args.I:
            self.buildIndelMtx()
        if self.args.S:
            self.cleanSnpMtx()
        if self.args.I:
            self.cleanIndelMtx()

    def buildIndelMtx(self):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building Indel Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building Indel Matrix")
        # takes intersect to find columns (pos) for matrix
        if self.args.F and self.args.G:
            self.__freebayesGatk(self.indel_positions, self.f_indel_trees, self.g_indel_trees)
        elif self.args.F:
            self.__freebayesOnly(self.indel_positions, self.f_indel_trees)
        elif self.args.G:
            self.__gatkOnly(self.indel_positions, self.g_indel_trees)

        if self.args.D and self.args.P:
            self.__dellyPindel(self.indel_positions, self.d_trees, self.p_trees)
        elif self.args.D:
            self.__dellyOnly(self.indel_positions, self.d_trees)
        elif self.args.P:
            self.__pindelOnly(self.indel_positions, self.p_trees)

        # build matrix
        self.indel_mtx = pd.DataFrame(np.zeros((len(self.samples), len(self.indel_positions))), index=self.samples, columns=sorted(self.indel_positions), dtype="int")

        # populate matrix
        if self.args.D:
            self.__populateIndelMtx(self.d_trees)
        if self.args.F:
            self.__populateIndelMtx(self.f_indel_trees)
        if self.args.G:
            self.__populateIndelMtx(self.g_indel_trees)
        if self.args.P:
            self.__populateIndelMtx(self.p_trees)
        self.__populateIndelMtx(self.fail_indel_trees)

    def buildSnpMtx(self):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building SNP Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building SNP Matrix")
        if self.args.F and self.args.G:
            self.__freebayesGatk(self.snp_positions, self.f_snp_trees, self.g_snp_trees)
        elif self.args.F:
            self.__freebayesOnly(self.snp_positions, self.f_snp_trees)
        elif self.args.G:
            self.__gatkOnly(self.snp_positions, self.g_snp_trees)

        self.snp_mtx = pd.DataFrame(np.zeros((len(self.samples), len(self.snp_positions))), index=self.samples, columns=sorted(self.snp_positions), dtype="int")

        if self.args.F:
            self.__populateSnpMtx(self.f_snp_trees)
        if self.args.G:
            self.__populateSnpMtx(self.g_snp_trees)
        self.__populateSnpMtx(self.fail_snp_trees)

    def __freebayesGatk(self, positions, f_trees, g_trees):
        """ get col names for matrix, taking intersect of Freebayes and GATK """
        for sample in self.samples:
            for chrom in g_trees[sample].keys():
                allNodes = g_trees[sample][chrom].getAllNodes()
                for node in allNodes:
                    try:
                        # check whether node position is within range of a node position in gatk tree
                        if f_trees[sample][chrom].searchTree(node.pos, LENIENCY_SMALL, f_trees[sample][chrom].root):
                            positions.add(node.var)
                    except KeyError: # if GATK tree has no SV calls for that chrom
                        pass

    def __freebayesOnly(self, positions, f_trees):
        """ get col names for matrix from Freebayes only """
        for sample in self.samples:
            for chrom in f_trees[sample].keys():
                node_list = f_trees[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)

    def __gatkOnly(self, positions, g_trees):
        """ get col names for matrix from GATK only """
        for sample in self.samples:
            for chrom in g_trees[sample].keys():
                node_list = g_trees[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)

    def __dellyPindel(self, positions, d_trees, p_trees):
        """ get col names for matrix, taking intersect of Delly and Pindel """
        for sample in self.samples:
            for chrom in d_trees[sample].keys():
                allNodes = d_trees[sample][chrom].getAllNodes()
                try:
                    for node in allNodes:
                        if p_trees[sample][chrom].searchTree(node.pos, LENIENCY_LARGE, p_trees[sample][chrom].root):
                            positions.add(node.var)
                except KeyError: # no SV for current chrom in Pindel tree
                    pass

    def __dellyOnly(self, positions, d_trees):
        """ get col names for matrix from Delly only """
        for sample in self.samples:
            for chrom in d_trees[sample].keys():
                node_list = d_trees[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)

    def __pindelOnly(self, positions, p_trees):
        """ get col names for matrix from Pindel only """
        for sample in self.samples:
            for chrom in p_trees[sample].keys():
                node_list = p_trees[sample][chrom].getAllNodes()
                for node in node_list:
                    positions.add(node.var)

    def __populateMtx(self, positions, mtx, tree, leniency):
        for sample in self.samples:
            for var in positions:
                chrom = var.chrom
                pos = var.pos
                try:
                    if tree[sample][chrom].searchTree(pos, LENIENCY_SMALL, tree[sample][chrom].root): # if pos in current tree
                        mtx.loc[sample, var] = var.alt # assign alternate
                    elif "{}:{}".format(chrom, pos) in self.low_cov[sample] and (mtx.loc[sample, var] == 0 or mtx.loc[sample, var] == var.ref):
                        mtx.loc[sample, var] = len(self.samples)*3 # assign low coverage
                    elif mtx.loc[sample, var] == 0:
                        mtx.loc[sample, var] = var.ref # assign reference
                except KeyError: # key error raised if no SVs detected in that sample on that chrom
                    if "{}:{}".format(chrom, pos) in self.low_cov[sample] and (mtx.loc[sample, var] == 0 or mtx.loc[sample, var] == var.ref):
                        mtx.loc[sample, var] = len(self.samples)*3
                    elif mtx.loc[sample, var] == 0:
                        mtx.loc[sample, var] = var.ref # reference

    def __populateSnpMtx(self, tree):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Populating SNP Matrix {} \n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Populating SNP Matrix")
        for sample in self.samples:
            f.write("Populating {} at {} \n".format(sample, datetime.datetime.time(datetime.datetime.now())))
            for var in self.snp_positions:
                try:
                    if tree[sample][var.chrom].searchTree(var.pos, LENIENCY_SNP, tree[sample][var.chrom].root):
                        # self.snp_mtx.loc[sample, var] = var.alt
                        self.snp_mtx.loc[sample, var] = 1
                    elif "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and (self.snp_mtx.loc[sample, var] == 0 or self.snp_mtx.loc[sample, var] == var.ref):
                        self.snp_mtx.loc[sample, var] = len(self.samples)*3
                    # elif self.snp_mtx.loc[sample, var] == 0:
                        # self.snp_mtx.loc[sample, var] = var.ref
                except KeyError:
                    if "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and (self.snp_mtx.loc[sample, var] == 0 or self.snp_mtx.loc[sample, var] == var.ref):
                        self.snp_mtx.loc[sample, var] = len(self.samples)*3
                    # elif self.snp_mtx.loc[sample, var] == 0:
                        # self.snp_mtx.loc[sample, var] = var.ref

    def __populateIndelMtx(self, tree):
        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building Indel Matrix {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Populating Indel Matrix")
        for sample in self.samples:
            f.write("Populating {} at {}\n".format(sample, datetime.datetime.time(datetime.datetime.now())))
            for var in self.indel_positions:
                try:
                    if tree[sample][var.chrom].searchTree(var.pos, LENIENCY_SMALL, tree[sample][var.chrom].root):
                        self.indel_mtx.loc[sample, var] = 1
                    elif "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and self.indel_mtx.loc[sample, var] == 0:
                        self.indel_mtx.loc[sample, var] = len(self.samples)*3
                except KeyError:
                    if "{}:{}".format(var.chrom, var.pos) in self.low_cov[sample] and self.indel_mtx.loc[sample, var] == 0:
                        self.indel_mtx.loc[sample, var] = len(self.samples)*3

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

    def getDifference(self):
        var_dif = {}
        # confirmed variant BST (combine indel and snp as some crossover in calls)
        var_pos_nodes = {}
        for var in self.snp_positions:
            if var_pos_nodes.get(var.chrom) is None:
                var_pos_nodes[var.chrom] = [Node(var.pos, var)]
            else:
                var_pos_nodes[var.chrom].append(Node(var.pos, var))
        for var in self.indel_positions:
            if var_pos_nodes.get(var.chrom) is None:
                var_pos_nodes[var.chrom] = [Node(var.pos, var)]
            else:
                var_pos_nodes[var.chrom].append(Node(var.pos, var))
        var_bst = buildTrees(var_pos_nodes)

        # For each tools' tree results, search confirmed positions
        for sample in self.samples:
            var_dif[sample] = []
            # Delly
            for chrom in self.d_trees[sample].keys():
                try:
                      for node in self.d_trees[sample][chrom].getAllNodes():
                        if not var_bst[chrom].searchTree(node.pos, LENIENCY_LARGE, var_bst[chrom].root):
                            var_dif[sample].append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
                except KeyError:
                    pass
        # freebayes
        fb_snp_dif = []
        fb_indel_dif = []
        for chrom in self.f_snp_trees[sample].keys():
            try:
                for node in self.f_snp_trees[sample][chrom].getAllNodes():
                    if not var_bst[chrom].searchTree(node.pos, LENIENCY_SMALL, var_bst[chrom].root):
                        fb_snp_dif.append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
            except KeyError:
                pass
        for chrom in self.f_indel_trees[sample].keys():
            try:
                for node in self.f_indel_trees[sample][chrom].getAllNodes():
                    if not var_bst[chrom].searchTree(node.pos, LENIENCY_SMALL, var_bst[chrom].root):
                        fb_indel_dif.append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
            except KeyError:
                pass
        # Find only calls not found in either matrix
        var_dif[sample].extend(list(set(fb_snp_dif) & (set(fb_indel_dif))))

        # GATK
        for chrom in self.g_snp_trees[sample].keys():
            try:
                for node in self.g_snp_trees[sample][chrom].getAllNodes():
                    if not var_bst[chrom].searchTree(node.pos, LENIENCY_SMALL, var_bst[chrom].root):
                        var_dif[sample].append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
            except KeyError:
                pass
        for chrom in self.g_indel_trees[sample].keys():
            try:
                for node in self.g_indel_trees[sample][chrom].getAllNodes():
                    if not var_bst[chrom].searchTree(node.pos, LENIENCY_SMALL, var_bst[chrom].root):
                        var_dif[sample].append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
            except KeyError:
                pass

        # PINDEL
        for chrom in self.p_trees[sample].keys():
            try:
                for node in self.p_trees[sample][chrom].getAllNodes():
                    if not var_bst[chrom].searchTree(node.pos, LENIENCY_SMALL, var_bst[chrom].root):
                        var_dif[sample].append("{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.tool))
            except KeyError:
                pass


        return var_dif


class AllTrees:
    """ Build BST for each tool (using filtered bcf output for each tool)
    for each sample """
    def __init__(self, args, sample_names):
        self.args = args
        self.sample_names = sample_names
        self.delly_trees = {}
        self.freebayes_snp_trees = {}
        self.freebayes_indel_trees = {}
        self.gatk_snp_trees = {}
        self.gatk_indel_trees = {}
        self.pindel_trees = {}
        self.del_end_trees = {}
        self.pindel_td_trees = {}
        self.td_end_trees = {}
        self.long_insertions = {}
        self.fail_snp_trees = {}
        self.fail_indel_trees = {}

        f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
        f.write("Building BST {}\n".format(datetime.datetime.time(datetime.datetime.now())))
        print("Building BST")
        self.buildAllTrees()
        self.sort_TD()


    def buildAllTrees(self):
        """ Build BST for all variant caller results """
        for sample in self.sample_names:
            f = open("{}/2_SVs/collate_status.txt".format(self.args.o), 'a')
            f.write("Working on sample {} ({}/{}) {}\n".format(sample, self.sample_names.index(sample)+1, len(self.sample_names), datetime.datetime.time(datetime.datetime.now())))
            print("Working on sample {} ({}/{})".format(sample, self.sample_names.index(sample)+1, len(self.sample_names)))
            # store 'failed' SV calls for each tool for each SV type (snp and indel/other)
            snp_fail_nodes = {}
            indel_fail_nodes = {}

            if self.args.D and self.args.I: # Delly
                snp_nodes, indel_nodes = ReadBCF("Delly").read("{0}/2_SVs/Delly/{1}_Results/{1}.PASS.bcf".format(self.args.o, sample))
                self.delly_trees[sample] = self.buildTrees(indel_nodes)

                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Delly").read("{0}/2_SVs/Delly/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.F: # Freebayes
                snp_nodes, indel_nodes = ReadBCF("Freebayes").read("{0}/2_SVs/Freebayes/{1}_Results/{1}.PASS.bcf".format(self.args.o, sample))
                if self.args.S:
                    self.freebayes_snp_trees[sample] = self.buildTrees(snp_nodes)
                if self.args.I:
                    self.freebayes_indel_trees[sample] = self.buildTrees(indel_nodes)

                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("Freebayes").read("{0}/2_SVs/Freebayes/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.G: # gatk
                snp_nodes, indel_nodes = ReadBCF("GATK").read("{0}/2_SVs/GATK/{1}_Results/{1}.PASS.bcf".format(self.args.o, sample))
                if self.args.S:
                    self.gatk_snp_trees[sample] = self.buildTrees(snp_nodes)
                if self.args.I:
                    self.gatk_indel_trees[sample] = self.buildTrees(indel_nodes)

                tmp_snp_fail = {}
                tmp_indel_fail = {}
                tmp_snp_fail, tmp_indel_fail = ReadBCF("GATK").read("{0}/2_SVs/GATK/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                snp_fail_nodes = self.collateFailedNodes(snp_fail_nodes, tmp_snp_fail)
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

            if self.args.P and self.args.I: # Pindel
                indel_nodes, end_nodes = ReadBCF("Pindel").read("{0}/2_SVs/Pindel/{1}_Results/{1}.PASS.bcf".format(self.args.o, sample))
                self.pindel_trees[sample] = self.buildTrees(indel_nodes)
                self.del_end_trees[sample] = self.buildTrees(end_nodes)

                tmp_indel_fail = {}
                tmp_end_fail = {}
                tmp_indel_fail, tmp_end_fail = ReadBCF("Pindel").read("{0}/2_SVs/Pindel/{1}_Results/{1}.FAIL.bcf".format(self.args.o, sample))
                indel_fail_nodes = self.collateFailedNodes(indel_fail_nodes, tmp_indel_fail)

                # Tandem duplication results stored separately
                # Also need to build tree of end positions to calculate duplications from translocations
                td_nodes, td_end_nodes = ReadBCF("Pindel").read("{0}/2_SVs/Pindel/{1}_Results/{1}.TD.bcf".format(self.args.o, sample))
                self.pindel_td_trees[sample] = self.buildTrees(td_nodes)
                self.td_end_trees[sample] = self.buildTrees(td_end_nodes)

            if self.args.S:
                self.fail_snp_trees[sample] = self.buildTrees(snp_fail_nodes)
            if self.args.I:
                self.fail_indel_trees[sample] = self.buildTrees(indel_fail_nodes)

    def sort_TD(self):
        """ Check for duplications and locations which are reported by
        Pindel as deletion and tandem duplication signals. Then removes the
        mis-reported 'deletion' signal from both Pindel and Delly trees """

        if self.args.P:
            print("Finding translocations...")
            for sample in self.sample_names:
                for chrom in self.pindel_td_trees[sample].keys():
                    for node in self.pindel_td_trees[sample][chrom].getAllNodes():
                        # find TD start == DEL1 start
                        try:
                            del1_nodes = self.pindel_trees[sample][chrom].getNodes([], node.var.pos, 10, self.pindel_trees[sample][chrom].root)
                        except KeyError: # no variants in that chrom for that sample
                            del1_nodes = []
                        for del1 in del1_nodes:
                            # find DEL1 end == DEL2 start
                            del2_nodes = self.pindel_trees[sample][chrom].getNodes([], del1.var.end, 10, self.pindel_trees[sample][chrom].root)
                            for del2 in del2_nodes:
                                # check DEL2 end == TD end
                                if del2.var.end in range(node.var.end-10, node.var.end+10):
                                    print("Translocation found in {}:{} original start: {}, new start: {}".format(sample, chrom, del2.var.pos, del1.var.pos))

                                    self.pindel_trees[sample][chrom].delete(del2) # del2 is original position, remove from tree

                                    # change SV type to translocation in Pindel results
                                    del1.var.type = "TSLN"

                                    # change SV type to translocation in Delly results
                                    try:
                                        tsln_nodes = self.delly_trees[sample][chrom].getNodes([], del1.var.pos, 5, self.delly_trees[sample][chrom].root)
                                        for tsln_node in tsln_nodes:
                                            tsln_node.var.type = "TSLN"

                                    except KeyError:
                                        pass

                                    # remove del2 from delly tree
                                    try:
                                        delly_nodes = self.delly_trees[sample][chrom].getNodes([], del2.var.pos, 5, self.delly_trees[sample][chrom].root)
                                        for delly_node in delly_nodes:
                                            self.delly_trees[sample][chrom].delete(delly_node)
                                    except KeyError: # no SV calls on that chrom in Delly
                                        pass

                                    # remove del2 from fail tree
                                    try:
                                        fail_nodes = self.fail_indel_trees[sample][chrom].getNodes([], del2.var.pos, 5, self.fail_indel_trees[sample][chrom].root)
                                        for fail_node in fail_nodes:
                                            self.fail_indel_trees[sample][chrom].delete(fail_node)
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


class RecombinationPos:
    """ store unique start and end ranges in ascending order of start pos
    params:
        gff_rec = gff_parser.Gff object
        samples = re pattern object of sample names """
    def __init__(self, gff_rec, samples):
        self.start = []
        self.end = []
        self.add(gff_rec, samples)

    def add(self, gff_rec, samples):
        for record in gff_rec.records:
            # performs a pattern search for sample name in gff so can analyse
            # any subset of samples from a full gubbins recombination file
            if samples.search(record.attributes["taxa"]) != None:
                self.__add(record.start, record.end)

    def __add(self, start, end):
        if len(self.start) == 0:
            self.start.append(start)
            self.end.append(end)
        else: # only need to check last value as input gff object is sorted by start
            if start in range(self.start[-1], self.end[-1]):
                if end > self.end[-1]:
                    self.end[-1] = end
            else:
                self.start.append(start)
                self.end.append(end)


"""
Additional Functions
"""

def get_sample_names(args):
    sample_names = []
    with open(args.s, "r") as infile:
        for line in infile:
            sample_names.append(line.strip())
    return sample_names


def getLowCov(args, sample_names):
    low_cov = {}
    for sample in sample_names:
        f = open("{}/2_SVs/collate_status.txt".format(args.o), 'a')
        f.write("Getting low coverage regions for {} ({}/{}) {}\n".format(sample, sample_names.index(sample)+1, len(sample_names), datetime.datetime.time(datetime.datetime.now())))
        print("Getting low coverage regions for {} ({}/{})".format(sample, sample_names.index(sample)+1, len(sample_names)))
        low_cov[sample] = set()
        with open("{}/1_Mapping/{}_low_cov.mpileup".format(args.o, sample)) as tsv:
            for line in csv.reader(tsv, delimiter="\t"):
                    pos = "{}:{}".format(line[0], line[1])
                    low_cov[sample].add(pos)
    return low_cov


def getIndex(alist, target):
    """ Fast search for index of target in a sorted list """
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


def saveSVDetails(file, sample, positions, tree):
    outfile = open(file, "a")
    for position in positions:
        chrom, pos = position.split(":")
        try:
            nodeList = tree[sample][chrom].getNodes([], int(pos), LENIENCY_SMALL, tree[sample][chrom].root)
            for node in nodeList:
                outfile.write("{}\t{}\t{}\t{}\t{}\n".format(node.var.chrom, node.var.pos, node.var.type, node.var.length, node.var.tool))
        except KeyError: # key error thrown if no SVs for that chrom in that sample
            pass
    outfile.close()


def remove_recombination(file, mtx):
    """ Remove regions from matrices based on positions in a gff file """
    gff_records = gff_parser.Gff(file)
    sample_names = mtx.index.values.tolist()
    samples = re.compile('|'.join(r'\b{}\b'.format(name) for name in sample_names))
    rec_pos = RecombinationPos(gff_records, samples)
    print("Recombination start positions",rec_pos.start)
    print("Recombination end positions",rec_pos.end)
    to_drop = []
    sorted_cols = sorted(mtx.columns)
    for var in sorted_cols:
        to_del = [] # delete recombination positions to maintain search speed in higher values
        for idx in range(len(rec_pos.start)): # search start pos for recombination
            if var.pos in range(rec_pos.start[idx], rec_pos.end[idx]+1):
                to_drop.append(var)
                break
            elif var.pos > rec_pos.end[idx]:
                to_del.append(idx)
            elif var.pos < rec_pos.start[idx]:
                break

        to_del.reverse()
        for d in to_del:
            rec_pos.start.pop(d)
            rec_pos.end.pop(d)

    print("Removing {} recombination loci".format(len(to_drop)))

    mtx.drop(to_drop, axis = 1, inplace = True)


def writeSnpNexus(pwd, sample_names, num_snps, snp_cols_str, prefix):
    """
    params:
        pwd: working directory for pipeline results
        sample_names: list of sample names
        num_snps: number of unique snps
        snp_cols_str: string representation of each column from snp matrix
    """
    with open("{}/3_Trees/{}_snp.nex".format(pwd, prefix), "w") as nex:
        nex.write("#NEXUS\n\n")
        nex.write("Begin data;\n")
        nex.write("\tDimensions ntax={} nchar={};\n".format(len(sample_names), num_snps))
        nex.write("\tFormat datatype=DNA interleave=no gap=- missing=?;\n")
        nex.write("\tMatrix\n")
        nex_string = ""
        for n in range(len(sample_names)):
            nex_string += "{}\t{}\n".format(sample_names[n], snp_cols_str[n])
        nex.write(nex_string)
        nex.write("\t;\n")
        nex.write("End;\n\n")
        nex.write("begin mrbayes;\n\n")
        nex.write("\tset autoclose=yes;\n")
        nex.write("\tmcmc ngen=1500000 nchains=4 samplefreq=100 burnin=5000;\n")
        nex.write("\tlset nucmodel=4by4 rates=gamma ngammacat=4;\n")
        nex.write("\tsump burnin=5000;\n")
        nex.write("\tsumt burnin=5000 conformat=simple contype=allcompat;\n\n")
        nex.write("end;")
    nex.close()


def writeIndelNexus(pwd, sample_names, num_indels, indel_cols_str, prefix):
    """
    params:
        pwd: working directory for pipeline results
        sample_names: list of sample names
        num_indels: number of unique snps
        indel_cols_str: string representation of each column from snp matrix
    """
    with open("{}/3_Trees/{}_indel.nex".format(pwd, prefix), "w") as nex:
        nex.write("#NEXUS\n\n")
        nex.write("Begin data;\n")
        nex.write("\tDimensions ntax={} nchar={};\n".format(len(sample_names), num_indels))
        nex.write("\tFormat datatype=restriction interleave=no gap=- missing=?;\n")
        nex.write("\tMatrix\n")
        nex_string = ""
        for n in range(len(indel_cols_str)):
          nex_string += "{}\t{}\n".format(sample_names[n], indel_cols_str[n])
        nex.write(nex_string)
        nex.write("\t;\n")
        nex.write("End;\n\n")
        nex.write("begin mrbayes;\n\n")
        nex.write("\tset autoclose=yes;\n")
        nex.write("\tmcmc ngen=1500000 nchains=4 samplefreq=100 burnin=5000;\n")
        nex.write("\tlset nucmodel=4by4 rates=gamma ngammacat=4;\n")
        nex.write("\tsump burnin=5000;\n")
        nex.write("\tsumt burnin=5000 conformat=simple contype=allcompat;\n\n")
        nex.write("end;")
    nex.close()


def writeCombinedNexus(pwd, sample_names, num_snps, num_indels, snp_cols_str, indel_cols_str, prefix):
    with open("{}/3_Trees/{}.nex".format(pwd, prefix), "w") as nex:
        nex.write("#NEXUS\n\n")
        nex.write("Begin data;\n")
        nex.write("\tDimensions ntax={} nchar={};\n".format(len(sample_names), num_snps+num_indels))
        nex.write("\tFormat datatype=mixed(DNA:1-{},Restriction:{}-{}) interleave=yes gap=- missing=?;\n".format(num_snps, num_snps+1, num_snps+num_indels))
        nex.write("\tMatrix\n")
        nex_string = ""
        for n in range(len(snp_cols_str)):
            nex_string += "{}\t{}\n".format(sample_names[n], snp_cols_str[n])
        nex.write(nex_string)
        nex.write("\n")
        nex_string = ""
        for n in range(len(indel_cols_str)):
            nex_string += "{}\t{}\n".format(sample_names[n], indel_cols_str[n])
        nex.write(nex_string)
        nex.write("\t;\n")
        nex.write("End;\n\n")
        nex.write("begin mrbayes;\n\n")
        nex.write("\tcharset DNA = 1-{};\n".format(num_snps))
        nex.write("\tcharset restriction = {}-{};\n".format(num_snps+1, num_snps+num_indels))
        nex.write("\tpartition favored = 2: DNA, restriction;\n")
        nex.write("\tset autoclose=yes;\n")
        nex.write("\tmcmc ngen=1500000 nchains=4 samplefreq=100 burnin=5000;\n")
        nex.write("\tlset nucmodel=4by4 rates=gamma ngammacat=4;\n")
        nex.write("\tsump burnin=5000;\n")
        nex.write("\tsumt burnin=5000 conformat=simple contype=allcompat;\n\n")
        nex.write("end;")
    nex.close()


def writeSnpPhy(pwd, sample_names, num_snps, snp_cols_str, prefix):
    with open("{}/3_Trees/{}_snp.phy".format(pwd, prefix), "w") as phy:
        phy.write("{} {}\n".format(len(sample_names), num_snps))
        phy_string = ""
        for n in range(len(sample_names)):
            phy_string += "{}\t{}\n".format(sample_names[n], snp_cols_str[n])
        phy.write(phy_string)
    phy.close()


def writeIndelPhy(pwd, sample_names, num_indels, indel_cols_str, prefix):
    with open("{}/3_Trees/{}_indel.phy".format(pwd, prefix), "w") as phy:
        phy.write("{} {}\n".format(len(sample_names), num_indels))
        phy_string = ""
        for n in range(len(sample_names)):
            phy_string += "{}\t{}\n".format(sample_names[n], indel_cols_str[n])
        phy.write(phy_string)
        phy.close()


def writeCombinedPhy(pwd, sample_names, num_snps, num_indels, snp_cols_str, indel_cols_str, prefix):
    with open("{}/3_Trees/{}.phy".format(pwd, prefix), "w") as phy:
        phy.write("{} {}\n".format(len(sample_names), num_snps+num_indels))
        phy_string = ""
        for n in range(len(sample_names)):
            phy_string += "{}\t{} {}\n".format(sample_names[n], snp_cols_str[n], indel_cols_str[n])
        phy.write(phy_string)
        phy.close()
    with open("{}/3_Trees/{}.partition".format(pwd, prefix), "w") as part:
        part.write("DNA, p1 = 1-{}\n".format(num_snps))
        part.write("BIN, p2 = {}-{}\n".format(num_snps+1, num_snps+num_indels))
        part.close()


def write_vcf(args, var_mtx, file_name, sample_names, tmp_file):
    """ Writes variant matrix to a vcf file """

    tool_string = ""
    if args.D:
        tool_string += "Delly;"
    if args.F:
        tool_string += "Freebayes;"
    if args.G:
        tool_string += "GATK-HaplotypeCaller;"
    if args.P:
        tool_string += "Pindel;"

    mtx_pos = sorted(var_mtx.columns.values.tolist())

    # split matrix into blocks
    pos_unsorted = var_mtx.columns.values.tolist()
    block_idx = []
    pos_blocks = []
    pos_blocks.append(mtx_pos[0:2000])
    for n in range(2000,len(pos_unsorted),2000):
        block_idx.append(n)
        pos_blocks.append(mtx_pos[n:n+2000])
    mtx_blocks = np.hsplit(var_mtx.as_matrix(), block_idx)

    vcf = open("{}/2_SVs/{}".format(args.o, file_name), "w")
    vcf.write("##fileformat=VCFv4.2\n")
    vcf.write("##fileDate={}\n".format(datetime.datetime.today().strftime('%Y%m%d')))
    vcf.write("##source=SV_pipeline;Tools={}\n".format(tool_string))
    vcf.write("##reference={}\n".format(args.r))
    vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
    vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant detected">\n')
    vcf.write('##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of variant region">\n')
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    vcf.close()

    n = 0
    for b in range(len(pos_blocks)):
        tmp = open("{}/2_SVs/{}".format(args.o, tmp_file), "w")
        n += 1
        print("Writing variant block {}/{} {}".format(n, len(pos_blocks), datetime.datetime.time(datetime.datetime.now())))
        vcf_strings = []
        for var in pos_blocks[b]:
            if var.ref == 0:
                vcf_strings.append("{}\t{}\t.\t*\t*\t.\t.\t{}\t.".format(var.chrom, var.pos, "END={};SVTYPE={};LEN={};".format(var.end, var.type, var.length)))
            else:
                # sample_field.replace("0", str(var.ref)).replace("1", str(var.alt))
                vcf_strings.append("{}\t{}\t.\t{}\t{}\t.\t.\t{}\t.".format(var.chrom, var.pos, var.ref, var.alt, "END={};SVTYPE={};LEN={};".format(var.end, var.type, var.length)))
        vcf_string = "\n".join(string for string in vcf_strings)
        tmp.write(vcf_string)
        tmp.close()
        os.system("cat {}/2_SVs/{} >> {}/2_SVs/{}".format(args.o, tmp_file, args.o, file_name))
        # cleanup
        sample_field = None
        vcf_strings = None
        vcf_string = None
        tmp = None
    os.system("rm {}/2_SVs/{}".format(args.o, tmp_file))


def write_phylo(args, sample_names, prefix, snp_mtx = None, indel_mtx = None):
    """ Make .phy and .nexus files """

    num_snps = 0
    num_indels = 0
    # string representation of each column to make phylo files
    if snp_mtx is not None and len(snp_mtx) > 0:
        # binary representation of snps
    #     snp_cols = snp_mtx.to_string(header=False, index=False).split('\n') # string representation of each col
    #     snp_cols = [vals.replace("N", "?") for vals in snp_cols] # replace any low cov pos with "?"
    #     snp_cols_str = ["".join(vals.split()) for vals in snp_cols] # list of string reprs for col values in each sample

        # Converting to character representation
        snp_cols_char = []
        for n in range(len(sample_names)):
            snp_cols_char.append("")
        for var in snp_mtx.columns.values.tolist():
            vals = "".join(snp_mtx[var].to_string(header=False, index=False).split("\n"))
            vals = vals.replace("0", var.ref).replace("1", var.alt)
            for v in range(len(vals)):
                snp_cols_char[v] += vals[v]
        num_snps = len(snp_mtx.columns.values)


    if indel_mtx is not None and len(indel_mtx) > 0:
        # string representation of each column to make nexus file
        indel_cols = indel_mtx.to_string(header=False, index=False).split('\n')
        indel_cols = [vals.replace("N", "?") for vals in indel_cols]
        indel_cols_str = ["".join(vals.split()) for vals in indel_cols]
        num_indels = len(indel_mtx.columns.values)

    if num_snps > 0:
        writeSnpNexus(args.o, sample_names, num_snps, snp_cols_char, prefix)
        writeSnpPhy(args.o, sample_names, num_snps, snp_cols_char, prefix)
    if num_indels > 0:
        writeIndelNexus(args.o, sample_names, num_indels, indel_cols_str, prefix)
        writeIndelPhy(args.o, sample_names, num_indels, indel_cols_str, prefix)
    if num_snps > 0 and num_indels > 0:
        writeCombinedNexus(args.o, sample_names, num_snps, num_indels, snp_cols_char, indel_cols_str, prefix)
        writeCombinedPhy(args.o, sample_names, num_snps, num_indels, snp_cols_char, indel_cols_str, prefix)
