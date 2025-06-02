import argparse
import pysam

def count_and_print_reads_in_interval(bam_file, chrom, pos1, pos2):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    count = 0
    total_reads = 0
    for read in samfile.fetch(chrom, pos1-1, pos2):
        # Ensure the read covers the whole interval

        if read.reference_start <= pos1 and read.reference_end >= pos2:
            total_reads += 1
            # Check the positions are not clipped
            # print(read.query_name, read.reference_start, read.preference_end, read.cigarstring)

            cigar_tuples = read.cigartuples
            read_start = read.reference_start+1
            pos1_matched = False
            pos2_matched = False
            for (op, length) in cigar_tuples:
                # print(op, length)
                if op not in [0, 7]:  # 0 and 7 correspond to match operations
                    read_start += length
                    continue
                # print("parse cigar:",read_start, pos1, read_start + length)
                if read_start <= pos1 < read_start + length:
                    pos1_matched = True
                if read_start <= pos2 < read_start + length:
                    pos2_matched = True
                read_start += length
            if pos1_matched and pos2_matched:
                # print(read.query_name, read.cigarstring)
                count += 1
            # else:
                # print("Read {} does not cover the interval [{}, {}]".format(read.query_name, pos1, pos2))
    samfile.close()
    return count, total_reads
# def perfect_match_at_position(read, position):
#     cigar_tuples = read.cigartuples
#     read_start = read.reference_start+1

#     for (op, length) in cigar_tuples:
#         if op not in [0, 7]:  # 0 and 7 correspond to match operations
#             read_start += length
#             continue
#         if read_start <= position < read_start + length:
#             return True
#         read_start += length

#     return False

# def count_and_print_reads_in_interval(bam_file, chrom, pos1, pos2):
#     samfile = pysam.AlignmentFile(bam_file, "rb")
#     count = 0
#     total_reads = 0

#     for read in samfile.fetch(chrom, pos1-1, pos2):
#         # Ensure the read covers the whole interval
#         if read.reference_start <= pos1 and read.reference_end >= pos2:
#             total_reads += 1

#             if perfect_match_at_position(read, pos1):
#                 if not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
#                     mate = next(samfile.fetch(chrom, read.next_reference_start, read.next_reference_start+1), None)
#                     if mate and perfect_match_at_position(mate, pos2):
#                         count += 1

#     samfile.close()
    return count, total_reads
def main():
    parser = argparse.ArgumentParser(description='Count reads in a given interval.')
    parser.add_argument('bam_file', type=str, help='Path to BAM file.')
    parser.add_argument('chrom', type=str, help='Chromosome name.')
    parser.add_argument('pos1', type=int, help='Start position of the interval.')
    parser.add_argument('pos2', type=int, help='End position of the interval.')
    args = parser.parse_args()

    num_reads, total_reads = count_and_print_reads_in_interval(args.bam_file, args.chrom, args.pos1, args.pos2)
    print(f"Number of reads spanning the interval [{args.pos1}, {args.pos2}] on {args.chrom}: {num_reads}")
    print(f"Total number of reads spanning the interval [{args.pos1}, {args.pos2}] on {args.chrom}: {total_reads}")

if __name__ == '__main__':
    main()