import pandas as pd

import csv
import argparse

def get_args():
    """
    parses arguments from the command line
    :return: list of command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='Input config template',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Output file',
                        required=False)
    args = parser.parse_args()
    return args


def main():
        """
        This script generates several types of Transdecoder friendly gff files of Augustus gtf files without genes (only transcript+exon)
        """
        args = get_args()

        header=["seqid","source","type","start","end","score","strand","phase","attributes"]
        gtf = pd.read_csv(args.input, sep="\t", encoding='utf_8', names=header, header=None)
        #print(gtf['type'].unique())
        # gtf.loc[gtf['type'] == 'exon', 'attributes'] = gtf.loc[gtf['type'] == 'exon', 'attributes'] + gtf.loc[gtf['type'] == 'exon', 'attributes'].str.split(';').str[0].split('"').str[1]
        gtf.loc[gtf['type'] == 'exon', 'attributes'] = gtf.loc[gtf['type'] == 'exon', 'attributes'] + " gene_id \"" + gtf.loc[gtf['type'] == 'exon', 'attributes'].str.split(';').str[
                                                               0].str.split('\"').str[1].str.split('.').str[0] + "\";"
        print(gtf)
        gtf.to_csv(args.output, header=None, sep="\t", index=None, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
        main()

