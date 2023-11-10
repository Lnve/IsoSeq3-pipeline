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


    parser.add_argument('-ov',
                        '--MIN_PERCENT_OVERLAP',
                        type=int,
                        help='',
                        required=False,
                        default=50)
    parser.add_argument('-rp',
                        '--MIN_PERCENT_PROT_CODING',
                        type=int,
                        help='',
                        required=False,
                        default=40)
    parser.add_argument('-pp',
                        '--MIN_PERID_PROT_COMPARE',
                        type=int,
                        help='',
                        required=False,
                        default=70)
    parser.add_argument('-fl',
                        '--MIN_PERCENT_LENGTH_FL_COMPARE',
                        type=int,
                        help='',
                        required=False,
                        default=70)
    parser.add_argument('-nfl',
                        '--MIN_PERCENT_LENGTH_NONFL_COMPARE',
                        type=int,
                        help='',
                        required=False,
                        default=70)
    parser.add_argument('-pa',
                        '--MIN_PERCENT_ALIGN_LENGTH',
                        type=int,
                        help='',
                        required=False,
                        default=70)
    parser.add_argument('-pog',
                        '--MIN_PERCENT_OVERLAP_GENE_REPLACE',
                        type=int,
                        help='',
                        required=False,
                        default=80)
    parser.add_argument('-ue',
                        '--MAX_UTR_EXONS',
                        type=int,
                        help='',
                        required=False,
                        default=2)

    # parser.add_argument('-fo',
    #                     '--MIN_FL_ORF_SIZE',
    #                     type=int,
    #                     help='',
    #                     required=False,
    #                     default=1)
    parser.add_argument('-hp',
                        '--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE',
                        type=int,
                        help='',
                        required=False,
                        default=0)
    # parser.add_argument('-t',
    #                     '--TRUST_FL_STATUS',
    #                     type=int,
    #                     help='',
    #                     required=False,
    #                     default=1)
    # parser.add_argument('-c',
    #                     '--GENETIC_CODE',
    #                     type=str,
    #                     help='',
    #                     required=False,
    #                     default="universal",
    #                     choices=["universal", "Euplotes", "Tetrahymena", "Candida", "Acetabularia"])

    args = parser.parse_args()
    return args


def main():

    args = get_args()

    config = open(args.input, 'r')
    lines = config.readlines()

    linestart = "cDNA_annotation_comparer.dbi:--"

    new_conf=[]

    for line in lines:
        if line.startswith(linestart+"MIN_PERCENT_OVERLAP="):
            new_conf.append(line.replace("<__MIN_PERCENT_OVERLAP__>",str(args.MIN_PERCENT_OVERLAP)))
            continue
        if line.startswith(linestart+"MIN_PERCENT_PROT_CODING="):
            new_conf.append(line.replace("<__MIN_PERCENT_PROT_CODING__>",str(args.MIN_PERCENT_PROT_CODING)))
            continue
        if line.startswith(linestart+"MIN_PERID_PROT_COMPARE="):
            new_conf.append(line.replace("<__MIN_PERID_PROT_COMPARE__>",str(args.MIN_PERID_PROT_COMPARE)))
            continue
        if line.startswith(linestart+"MIN_PERCENT_LENGTH_FL_COMPARE="):
            new_conf.append(line.replace("<__MIN_PERCENT_LENGTH_FL_COMPARE__>",str(args.MIN_PERCENT_LENGTH_FL_COMPARE)))
            continue
        if line.startswith(linestart+"MIN_PERCENT_LENGTH_NONFL_COMPARE="):
            new_conf.append(line.replace("<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>",str(args.MIN_PERCENT_LENGTH_NONFL_COMPARE)))
            continue
        if line.startswith(linestart+"MIN_PERCENT_ALIGN_LENGTH="):
            new_conf.append(line.replace("<__MIN_PERCENT_ALIGN_LENGTH__>",str(args.MIN_PERCENT_ALIGN_LENGTH)))
            continue
        if line.startswith(linestart+"MIN_PERCENT_OVERLAP_GENE_REPLACE="):
            new_conf.append(line.replace("<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>",str(args.MIN_PERCENT_OVERLAP_GENE_REPLACE)))
            continue
        if line.startswith(linestart+"MAX_UTR_EXONS="):
            new_conf.append(line.replace("<__MAX_UTR_EXONS__>",str(args.MAX_UTR_EXONS)))
            continue
        # How do I set this? 1 and 0?
        # if line.startswith(linestart+"TRUST_FL_STATUS="):
        #    new_conf.append(line.replace("<__TRUST_FL_STATUS__>",str(args.TRUST_FL_STATUS)))
        # if line.startswith(linestart+"MIN_FL_ORF_SIZE="):
        #   new_conf.append(line.replace("<__MIN_FL_ORF_SIZE__>",str(args.MIN_FL_ORF_SIZE)))
        if line.startswith(linestart+"STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE="):
            new_conf.append(line.replace("<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>",str(args.STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE)))
            continue
        # if line.startswith(linestart + "GENETIC_CODE="):
        #    new_conf.append(line.replace("<__GENETIC_CODE__>", str(args.GENETIC_CODE)))
        new_conf.append(line)

    with open(f'{args.output}', 'w') as f:
        for line in new_conf:
            f.write(line)

if __name__ == '__main__':
        main()
