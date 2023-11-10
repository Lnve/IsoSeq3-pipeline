import argparse
import subprocess

def get_args():
    """
    parses arguments from the command line
    :return: list of command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='Input bam',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Output file',
                        required=True)
    parser.add_argument('-d',
                        '--ID-pattern',
                        type=str,
                        help='Pattern to use in the ID field',
                        required=True)
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        help='Number of threats this script asks for',
                        required=False,
                        default=1)
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    #command to get the input line
    RG = subprocess.run([f"samtools view -H {args.input} | grep \"^@RG\""], shell=True, capture_output=True, text=True)
    
    # make dict from RG header
    d = dict(field.split(":") for field in RG.stdout.strip("\n").split("\t")[1:])

    #create new RG header
    d["ID"] = args.ID_pattern

    # create addreplace command
    cmd = "samtools addreplacerg"
    for k, v in d.items():
        cmd += f" -r \"{k}:{v}\""

    flags = f" -w --threads {args.threads} --output-fmt BAM -o {args.output} {args.input}"
    cmd += flags

    # run subprocess
    print(cmd)
    subprocess.run([f"{cmd}"], shell=True)

if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

