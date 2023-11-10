import glob
import pandas as pd
import os.path
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def get_args():
    """
    parses arguments from the command line
    :return: list of command line arguments
    """
    parser = argparse.ArgumentParser()    
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Output file',
                        required=True)
    
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    augustus = sorted(glob.glob("output/v2/collapse_*/*_augustus_sqanti/*-vs-augustus_classification.txt"))
    aug_idx = []
    aug_df = pd.DataFrame()
    for i in augustus:
        aug_idx.append(os.path.basename(i))
        # print(os.path.basename(i))
        df = pd.read_csv(i, sep="\t")
        val = df['structural_category'].value_counts()
        # print(val)
        aug_df = pd.concat([aug_df, val], axis=1)
    
    aug_df.columns=aug_idx
    aug_df=aug_df.sort_index(axis=1)
    print(aug_df.fillna(0).T)
    aug_df.to_csv(args.output + "augustus_comparison")
    
    plt.figure(figsize=(4*len(aug_idx), 10))
    plt.plot(aug_df[:4].T, kind='bar')
    plt.savefig(args.output + "augustus_comparison_plot")


    liftoff = sorted(glob.glob("output/v2/collapse_*/*_liftoff_sqanti/*-vs-liftoff_classification.txt"))
    
    lo_idx = []
    lo_df = pd.DataFrame()
    for i in liftoff:
        lo_idx.append(os.path.basename(i))
        df = pd.read_csv(i, sep="\t")
        val = df['structural_category'].value_counts()
        lo_df = pd.concat([lo_df, val], axis=1)
    
    lo_df.columns=lo_idx
    lo_df=lo_df.sort_index(axis=1)
    print(lo_df.fillna(0).T)
    lo_df.to_csv(args.output + "liftoff_comparison")
    
    plt.figure(figsize=(10, 40))
    lo_p = lo_df[:4].T.plot(kind="bar")
    plt.savefig(args.output + "liftoff_comparison_plot")


if __name__ == '__main__':
    main()
