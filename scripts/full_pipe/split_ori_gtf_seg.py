import csv
import argparse


def load_intervals(interval_file):
    """
    load intervals from a file, skipping the first two lines and the last line.
    """
    intervals = []
    with open(interval_file, 'r') as f:
        lines = f.readlines()

    # skip the first two lines and the last line
    lines = lines[2:-1]

    for line in lines:
        start, end = line.strip().split()[0].split('-')
        start, end = int(start), int(end)
        intervals.append((start, end))
    return intervals


def split_exon_by_intervals(exon, intervals):
    """
    split a single exon by given intervals.
    """
    chrom, exon_start, exon_end, strand, attributes = exon
    split_exons = []
    current_start = exon_start

    # iterate over intervals and check if splitting is needed
    for interval_start, interval_end in sorted(intervals):
        if interval_start > exon_end:
            break  # subsequent intervals do not overlap with exon
        if interval_end < current_start:
            continue  # current interval is before exon

        # overlap detected, split the exon
        if interval_start <= exon_end and interval_end >= current_start:
            if current_start < interval_start:
                split_exons.append((current_start, interval_start - 1))
            split_exons.append((max(current_start, interval_start), min(exon_end, interval_end)))
            current_start = interval_end + 1

    # add any remaining part of the exon
    if current_start <= exon_end:
        split_exons.append((current_start, exon_end))

    return [(chrom, start, end, strand, attributes) for start, end in split_exons]


def process_gtf(gtf_file, intervals, output_file):
    """
    handle a GTF file, splitting exons based on intervals.
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        for row in reader:
            # skip comment lines
            if row[0].startswith('#'):
                writer.writerow(row)
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = row

            # handle non-exon features
            if feature != 'exon':
                writer.writerow(row)
                continue

            start, end = int(start), int(end)
            exon = (chrom, start, end, strand, attributes)


            # split the exon by chromosome intervals
            split_exons = split_exon_by_intervals(exon, intervals)

            # write the split exons
            for split_exon in split_exons:
                chrom, split_start, split_end, strand, attributes = split_exon
                writer.writerow([chrom, source, feature, split_start, split_end, score, strand, frame, attributes])


def process_gtf2(gtf_file, intervals, output_file):
    """
    Processes a GTF file, splitting exons based on intervals.
    Parameters:
        gtf_file : str
            Path to the input GTF file.
        intervals : list
            List of intervals grouped by chromosome.
        output_file : str
            Path to write the split GTF file.
    """
    with open(gtf_file, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')  # GTF files are tab-delimited
        
        for row in reader:
            # Skip comment lines
            if row[0].startswith('#'):
                outfile.write('\t'.join(row) + '\n')
                continue

            # Parse the GTF fields
            chrom, source, feature, start, end, score, strand, frame, attributes = row

            # Write non-exon features as-is
            if feature != 'exon':
                outfile.write('\t'.join(row) + '\n')
                continue

            # Process exon features
            start, end = int(start), int(end)
            exon = (chrom, start, end, strand, attributes)

            # Split exons based on intervals
            split_exons = split_exon_by_intervals(exon, intervals)

            # Write the split exons
            for split_exon in split_exons:
                chrom, split_start, split_end, strand, attributes = split_exon
                print(split_exon)
                # Write the line using proper GTF formatting
                outfile.write(
                    f"{chrom}\t{source}\t{feature}\t{split_start}\t{split_end}\t{score}\t{strand}\t{frame}\t{attributes}\n"
                )


def main():
    """
    main function to handle command-line arguments and execute the processing.
    """
    parser = argparse.ArgumentParser(description="split exons in a GTF file based on given intervals.")
    parser.add_argument("-g", "--gtf", required=True, help="input GTF file path")
    parser.add_argument("-i", "--intervals", required=True, help="intervals file path")
    parser.add_argument("-o", "--output", required=True, help="output GTF file path")

    args = parser.parse_args()

    # load intervals file
    intervals = load_intervals(args.intervals)

    # process GTF file
    process_gtf2(args.gtf, intervals, args.output)


if __name__ == "__main__":
    main()