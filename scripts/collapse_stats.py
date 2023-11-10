import argparse
import subprocess
import pandas as pd
import json
import glob
import os.path

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

    # process input files
    # ccs

    pb = sorted(glob.glob("output/v2/collapse_pb/*_augustus."))

    input_reads = []
    kept_ccs = []
    for file in ccs_in:
        content = open(file, 'r')
        lines = content.readlines()
        input = 0
        for line in lines:
            if line.startswith("ZMWs pass filters"):
                kept_ccs.append(int(line.strip("\n").split(":")[1].split(" ")[1]))
                input += int(line.strip("\n").split(":")[1].split(" ")[1])
            if line.startswith("ZMWs fail filters"):
                input += int(line.strip("\n").split(":")[1].split(" ")[1])
        input_reads.append(input)
    # print(input_reads)
    # print(kept_ccs)

    # lima

    lima_in = sorted(glob.glob("output/lima/pool*/*.demux.hifi.lima.summary"))
    kept_lima = []

    for file in lima_in:
        content = open(file, 'r')
        lines = content.readlines()
        for line in lines:
            if line.startswith("ZMWs above all thresholds"):
                kept_lima.append(line.strip("\n").split(":")[1].split(" ")[1])
    # print(kept_lima)

    # refine

    refine_in = sorted(glob.glob("output/refine/pool*/*.ncfl.filter_summary.report.json"))
    pool = []
    samples = []
    all = []
    kept_fl = []
    kept_polya = []

    for file in refine_in:
        pool.append(os.path.basename(file).split(".")[0])
        samples.append(os.path.basename(file).split(".")[1])
        with open(file, 'r') as f:
            d = json.load(f)
            all.append(d["attributes"][1]["value"])
            kept_fl.append(d["attributes"][2]["value"])
            kept_polya.append(d["attributes"][3]["value"])


    # prepare refine for pool dataframe
    refine_pool_dict = {}
    for idx in range(len(pool)):
        if pool[idx] in refine_pool_dict:
            refine_pool_dict[pool[idx]][0] += int(all[idx])
            refine_pool_dict[pool[idx]][1] += int(kept_fl[idx])
            refine_pool_dict[pool[idx]][2] += int(kept_polya[idx])
        else:
            refine_pool_dict[pool[idx]] = []
            refine_pool_dict[pool[idx]].append(int(all[idx]))
            refine_pool_dict[pool[idx]].append(int(kept_fl[idx]))
            refine_pool_dict[pool[idx]].append(int(kept_polya[idx]))
    print(refine_pool_dict)

    # prepare refine for sample dataframe

    refine_sample_dict = {}
    for idx in range(len(samples)):
        if samples[idx] in refine_sample_dict:
            refine_sample_dict[samples[idx]][0] += int(all[idx])
            refine_sample_dict[samples[idx]][1] += int(kept_fl[idx])
            refine_sample_dict[samples[idx]][2] += int(kept_polya[idx])
        else:
            refine_sample_dict[samples[idx]] = []
            refine_sample_dict[samples[idx]].append(int(all[idx]))
            refine_sample_dict[samples[idx]].append(int(kept_fl[idx]))
            refine_sample_dict[samples[idx]].append(int(kept_polya[idx]))
    print(refine_sample_dict)

    # get mapped reads pb
    mapped_pb = []
    mapped_mini = []

    for file in sorted(glob.glob("output/mapping_pb/*.bam")):
        print(file)
        mapped_pb.append(subprocess.run([f"samtools view -c -F 4 {file}"], shell=True, capture_output=True, text=True).stdout.strip("\n"))
    print(mapped_pb)

    # get mapped reads minimap
    for file in sorted(glob.glob("output/mapping_minimap/*.bam")):
        print(file)
        mapped_mini.append(subprocess.run([f"samtools view -c -F 4 {file}"], shell=True, capture_output=True, text=True).stdout.strip("\n"))
    print(mapped_mini)

    # create dataframes
    pool_df = pd.DataFrame({'CCS Input': input_reads,
                            'CCS kept': kept_ccs,
                            'Lima kept': kept_lima,
                            'Full-Length Reads': [x[0] for x in refine_pool_dict.values()],
                            'Full-Length Non-Chimeric Reads': [x[1] for x in refine_pool_dict.values()],
                            'Full-Length Non-Chimeric Reads with Poly-A Tail': [x[2] for x in refine_pool_dict.values()]})
    pool_df.index = sorted(list(set(pool)))

    print(pool_df)

    pool_df.to_csv(args.output + 'pipeline_pools.stats')

    sample_df = pd.DataFrame({'Full-Length Reads': [x[0] for x in refine_sample_dict.values()],
                               'Full-Length Non-Chimeric Reads': [x[1] for x in refine_sample_dict.values()],
                               'Full-Length Non-Chimeric Reads with Poly-A Tail': [x[2] for x in refine_sample_dict.values()]})
                               # 'Mapped by pbmm2': mapped_pb,
                               # 'Mapped by minimap2': mapped_mini})

    
    unique_list = []
    for s in samples:
        if s not in unique_list:
            unique_list.append(s)
    sample_df.index = unique_list
    
    print(sample_df)
    sample_df.to_csv(args.output + 'pipeline_samples.stats')

if __name__ == '__main__':
    main()
