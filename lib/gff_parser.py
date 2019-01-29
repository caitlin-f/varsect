""" Provides functionality for parsing gff files, specifically output from gubbins """

class Gff:
    def __init__(self, file):
        self.records = [] # list of record
        self.parse_lines(file)

    def parse_lines(self, file):
        with open(file) as input:
            for line in input:
                if line.startswith("#"):
                    pass
                else:
                    seqname, source, feature, start, end, score, strand, frame, attrs = line.strip().split("\t")
                    attributes = {}
                    for attr in attrs.strip(";").split(";"):
                        key, value = attr.split("=")
                        attributes[key] = value
                    self.records.append(Record(seqname, source, feature, int(start), int(end), float(score), strand, int(frame), attributes))
        self.records.sort()

class Record:
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attributes):
        self.seqname = seqname # name of chrom or scaffold
        self.source = source # program that generated feature / data source
        self.feature = feature # feature type
        self.start = start # start position
        self.end = end # end
        self.score = score # float
        self.strand = strand # + forward, - reverse
        self.frame = frame # 0, 1 or 2 indicating codon base feature starts
        self.attributes = attributes # dictionary of tag-value pairs containing info about each feature

    def __eq__(self, other):
        return self.start == other.start

    def __ne__(self, other):
        return self.start != other.start

    def __lt__(self, other):
        return self.start < other.start

    def __le__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__lt__(self, other)

    def __gt__(self, other):
        return self.start > other.start

    def __ge__(self, other):
        if self.__eq__(self, other):
            return True
        else:
            return self.__gt__(self, other)
