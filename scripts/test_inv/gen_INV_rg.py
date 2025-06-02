import argparse
from pybedtools import BedTool

def parse_args():
    parser = argparse.ArgumentParser(description='Find unique regions in BED file A.')
    parser.add_argument('-a', '--inputA', help='Input BED file A', required=True)
    parser.add_argument('-b', '--inputB', help='Input BED file B', required=True)
    parser.add_argument('-o', '--output', help='Output BED file', required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    # Load the bed files
    bedA = BedTool(args.inputA)
    bedB = BedTool(args.inputB)

    # Find unique regions in A that are not in B
    unique_to_A = bedA.subtract(bedB)

    # Save the unique regions to the output file
    unique_to_A.saveas(args.output)

if __name__ == '__main__':
    main()