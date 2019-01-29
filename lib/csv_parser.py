""" Provides functionality for reading in very large csv faster than the read.csv function in python """

import numpy as np

def extract_data(file):
    """ Args:
            file (string) filename of csv to read
        Returns:
            cols [] list of column names (header of csv)
            data {rownames:values} dictionary of matrix values by row
    """
    csv = open(file, "r")
    data = {}
    header = False
    for line in csv:
        if header is False:
            cols = line.strip("\n").split(",")[1:]
            col_vals = [int(c.split(":")[1]) for c in cols] # returns pos where colnames is ref:pos
            header = True
        else:
            rline = line.strip("\n").split(",")
            row = [2 if x == 'N' else x for x in rline[1:]]
            data[rline[0]] = row
    return col_vals, data

def build_matrix(cols, data):
    """ Args:
            cols [] list column names output by extract_data()
            data {row:values} dict of matrix values by row output by extract_data()
        Returns:
            numpy.array
    """
    mtx = np.zeros((len(data.keys()), len(cols)))
    for e, sample in enumerate(data.keys()):
        mtx[e] = data[sample]
    return mtx
