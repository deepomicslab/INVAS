import pysam
import sys

# accept command-line arguments
bam_file_path = sys.argv[1]
distance_threshold = int(sys.argv[2])
output_file_path = sys.argv[3]
min_reads_in_cluster = 10  # cluster中最少的reads数目

# define a function to process each cluster of reads
def process_cluster_reads(cluster):
    left_soft_clip_count = 0
    right_soft_clip_count = 0
    start = cluster['start']
    end = cluster['end']
    
    # check reads covering the start of the cluster, count left soft clips
    for read in cluster['reads']:
        if read.reference_start <= start < read.reference_end:
            # check if CIGAR starts with soft clip
            if read.cigartuples[0][0] == 4:
                left_soft_clip_count += 1

    # check reads covering the end of the cluster, count right soft clips
    for read in cluster['reads']:
        if read.reference_start < end <= read.reference_end:
            # check if CIGAR ends with soft clip
            if read.cigartuples[-1][0] == 4:
                right_soft_clip_count += 1

    # return the counts
    return left_soft_clip_count, right_soft_clip_count

# open the BAM file
bam_file = pysam.AlignmentFile(bam_file_path, "rb")

# for storing clusters of reads
clusters = []
prev_chrom = ""
# iterate through each read in the BAM file
for read in bam_file.fetch():
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    
    # get the chromosome of the read
    chrom = bam_file.get_reference_name(read.reference_id)
    if chrom != prev_chrom:
        print(f"Processing {chrom}")
        prev_chrom = chrom

    # try to place the read into an existing cluster
    placed_in_cluster = False
    for cluster in clusters:
        # check if chromosome matches and read's start position is within the cluster's end position plus threshold
        if read.reference_start <= cluster['end'] + distance_threshold and chrom == cluster['chrom']:
            cluster['reads'].append(read)
            cluster['end'] = max(cluster['end'], read.reference_end)
            placed_in_cluster = True
            break

    # if not placed in any cluster, create a new cluster
    if not placed_in_cluster:
        clusters.append({
            'chrom': chrom,
            'start': read.reference_start,
            'end': read.reference_end,
            'reads': [read]
        })

# handle the last cluster and write to output if it meets criteria
out_bed = open(output_file_path, 'w')
for cluster in clusters:
    if len(cluster['reads']) >= min_reads_in_cluster:
        left_clip_count, right_clip_count = process_cluster_reads(cluster)
        if left_clip_count > 0 and right_clip_count > 0:
            out_bed.write(f"{cluster['chrom']}\t{cluster['start']}\t{cluster['end']}\n")
            print(f"Cluster on {cluster['chrom']} from {cluster['start']} to {cluster['end']}, "
                f"contains {len(cluster['reads'])} reads, "
                f"with {left_clip_count} left soft clips, "
                f"and {right_clip_count} right soft clips.")

# close files
bam_file.close()
out_bed.close()