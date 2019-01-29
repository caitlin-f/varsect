import numpy as np
import pandas as pd
import argparse

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-i', type=str, required=True, help="Input .csv")
    parser.add_argument('-o', type=str, required=True, help="Output file")

    return parser.parse_args()

def extract_data(file):
    csv = open(file, "r")
    data = {}
    header = False
    for line in csv:
        if header is False:
            cols = line.strip("\n").split(",")[1:]
            header = True
        else:
            rline = line.strip("\n").split(",")
            row = [2 if x == 'N' else x for x in rline[1:]]
            data[rline[0]] = row
    return cols, data

def build_matrix(cols, data):
    mtx = np.zeros((len(data.keys()), len(cols)))
    for e, sample in enumerate(data.keys()):
        mtx[e] = data[sample]
    return pd.DataFrame(mtx, index=list(data.keys()), columns=cols, dtype=int)


def main():
    args = parse_args()

    cols, data = extract_data(args.i)
    mtx = build_matrix(cols, data)

    mtx.drop([col for col,val in mtx.sum().iteritems() if val <= 1], axis=1, inplace=True)

    mtx[mtx == 2] = 'NA'

    mtx.to_csv(args.o)


if __name__ == '__main__':
    main()
