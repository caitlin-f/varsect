""" Removes recombination regions of a snp/variant matrix using gubbins gff output
Requires: gubbins gff file sorted by start position """
import argparse
import gff_parser
from Bio import SeqIO
from csv_parser import *
import pandas as pd


class RecombinationPos:
    """ store unique start and end ranges in descending order of start pos"""
    def __init__(self):
        self.start = []
        self.end = []

    def add(self, start, end):
        if len(self.start) == 0:
            self.start.append(start)
            self.end.append(end)
        else:
            for i in range(1,len(self.start)+1):
                idx = -i # search from end of list (for speed)
                if start in range(self.start[idx], self.end[idx]):
                    if end > self.end[idx]:
                        self.end[idx] = end
                else:
                    self.start.append(start)
                    self.end.append(end)
                    break

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-i', type=str, required=True, help="input recombination embl file")
    parser.add_argument('-c', type=str, required=True, help="input csv matrix file")
    parser.add_argument('-o', type=str, required=True, help="Output csv matrix file")

    return parser.parse_args()

def main():
    args = parse_args()

    # read in gff files
    gff_records = gff_parser.Gff(args.i)

    # get start and end recombination positions
    pos = RecombinationPos()
    for record in gff_records.records:
        pos.add(record.start, record.end)


    # read in sample/variant matrix
    cols, data = extract_data(args.c)
    mtx = build_matrix(cols, data)
    samples = list(data.keys())


    # remove cols where pos falls in range of a recombination region
    df = pd.DataFrame(mtx, index=samples, columns=cols)

    to_drop = []
    sorted_cols = sorted(cols)
    for c in sorted_cols:
        to_del = [] # delete recombination positions to maintain search speed in higher values
        for idx in range(len(pos.start)): # search start pos for recombination
            if c in range(pos.start[idx], pos.end[idx]):
                to_drop.append(c)
                break
            elif c > pos.end[idx]:
                to_del.append(idx)
            elif c < pos.start[idx]:
                break

        to_del.reverse()
        for d in to_del:
            pos.start.pop(d)
            pos.end.pop(d)

    print(len(to_drop))
    df.drop(to_drop, axis = 1, inplace = True)

    df.to_csv(args.o)




if __name__ == '__main__':
    main()
