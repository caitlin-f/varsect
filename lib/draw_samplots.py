import argparse
import subprocess
from pysam import VariantFile

def parse_args():
    """ Parse user supplied arguments """
    parser = argparse.ArgumentParser()
    argparse.ArgumentParser(argument_default=False)
    parser.add_argument('-a', type=int, required=True, help="batch array number")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    parser.add_argument('-v', type=str, required=True, help="Name of input vcf")

    return parser.parse_args()

def get_sample_names(args):
    sample_names = []
    with open("{}/sample_names.txt".format(args.o),"r") as infile:
        for line in infile:
            sample_names.append(line.strip())
    return sample_names

def make_bam_arg(args, sample_names):
    bam_files = []
    start = args.a*30
    end = (args.a*30) + 30
    for sample in sample_names[start:end]:
        bam_files.append("{}/1_Mapping/{}.RG.bam".format(args.o, sample))

    return bam_files

def main():
    args = parse_args()
    sample_names = get_sample_names(args)
    samples = sample_names[args.a*30:(args.a*30)+30]
    bam_files = make_bam_arg(args, sample_names)
    print(samples)
    print(bam_files)

    bcf = VariantFile("{}".format(args.v))
    for rec in bcf.fetch():
        if rec.info["SVTYPE"] == "DEL":
            cmd = ['samplot.py']
            cmd.append('-n')
            for sample in samples:
                cmd.append(sample)
            cmd.append('-b')
            for bam in bam_files:
                cmd.append(bam)
            cmd.append('-o')
            cmd.append('{}/4_Images/DEL_{}.{}.png'.format(args.o, rec.pos, args.a))
            cmd.append('-c')
            cmd.append('{}'.format(rec.chrom))
            cmd.append('-t')
            cmd.append('DEL')
            cmd.append('-s')
            cmd.append('{}'.format(rec.pos))
            cmd.append('-e')
            cmd.append('{}'.format(rec.stop))
            cmd.append('-d')
            cmd.append('100')

            print(cmd)
            status = subprocess.call(cmd)

        elif rec.info["SVTYPE"] == "TSLN":
            cmd = ['samplot.py']
            cmd.append('-n')
            for sample in samples:
                cmd.append(sample)
            cmd.append('-b')
            for bam in bam_files:
                cmd.append(bam)
            cmd.append('-o')
            cmd.append('{}/4_Images/TSLN_inserted_{}.{}.png'.format(args.o, rec.pos, args.a))
            cmd.append('-c')
            cmd.append('{}'.format(rec.chrom))
            cmd.append('-t')
            cmd.append('Translocation - insertion position')
            cmd.append('-s')
            cmd.append('{}'.format(rec.pos))
            cmd.append('-e')
            cmd.append('{}'.format(rec.pos+1))
            cmd.append('-w')
            cmd.append('100')

            status = subprocess.call(cmd)

            cmd = ['samplot.py']
            cmd.append('-n')
            for sample in samples:
                cmd.append(sample)
            cmd.append('-b')
            for bam in bam_files:
                cmd.append(bam)
            cmd.append('-o')
            cmd.append('{}/4_Images/TSLN_origin_{}.{}.png'.format(args.o, rec.pos, args.a))
            cmd.append('-c')
            cmd.append('{}'.format(rec.chrom))
            cmd.append('-t')
            cmd.append('Translocation - original position')
            cmd.append('-s')
            cmd.append('{}'.format(rec.stop))
            cmd.append('-e')
            cmd.append('{}'.format(rec.stop+rec.info["LEN"]))
            cmd.append('-d')
            cmd.append('100')

            status = subprocess.call(cmd)


if __name__ == '__main__':
    main()
