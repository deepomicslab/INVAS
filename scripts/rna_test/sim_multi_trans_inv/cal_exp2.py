import pandas as pd
import pyreadr
import gffutils
import argparse
import sys

def create_db(gtf_file_path):
    """
    Creates a gffutils database from a GTF file.
    """
    db = gffutils.create_db(
        gtf_file_path,
        dbfn=':memory:',
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True
    )
    return db

def get_transcript_length(gtf_db, transcript_id):
    """
    Returns the length of a specified transcript using the gffutils database.
    
    :param gtf_db: A gffutils database object
    :param transcript_id: The ID of the transcript to find the length of
    :return: The length of the transcript or None if transcript not found
    """
    length = 0
    try:
        # Fetching the transcript using its ID
        transcript = gtf_db[transcript_id]
        # Loop through each feature of the transcript
        for exon in gtf_db.children(transcript, featuretype='exon', order_by='start'):
            # Calculate the length of the exon and add it to the total length
            length += exon.end - exon.start + 1
        return length
    except gffutils.exceptions.FeatureNotFoundError:
        return None  # Transcript not found

# 加载数据
def load_data(csv_file):
    return pd.read_csv(csv_file)

# 计算FPKM
def calculate_fpkm(df, total_mapped_reads):
    df['FPKM'] = (1e9 * df['count']) / (df['length'] * total_mapped_reads)
    return df

# 计算TPM
def calculate_tpm(df):
    df['TPM'] = (df['FPKM'] / df['FPKM'].sum()) * 1e6
    return df

def main(rda_files, gtf_files):
    # Create a dictionary to hold the total read counts per transcript
    total_read_counts = {}
    # Create a dictionary to track occurrences of transcript IDs
    transcript_occurrences = {}

    # Create an empty dataframe to hold expression data
    exp_df = pd.DataFrame(columns=['transcript_id', 'length', 'count'])

    # Process each pair of RDA and GTF files
    for rda_file, gtf_file in zip(rda_files, gtf_files):
        # Read the RDA file
        result = pyreadr.read_r(rda_file)
        counts_matrix_df = result[next(iter(result.keys()))]  # Assuming the first key contains the counts matrix
        trans_read_counts = counts_matrix_df.iloc[:, 0].to_dict()  # Assuming the first column is the sample of interest

        # Create a database from the GTF file
        gtf_db = create_db(gtf_file)

        # Accumulate read counts and transcript lengths
        for trans_id, read_count in trans_read_counts.items():
            # If transcript ID is already seen, increment the occurrence count and create a new ID
            if trans_id in transcript_occurrences:
                transcript_occurrences[trans_id] += 1
                new_trans_id = f"{trans_id}_{transcript_occurrences[trans_id]}"
            else:
                transcript_occurrences[trans_id] = 1
                new_trans_id = trans_id

            # Get the length of the transcript
            transcript_length = get_transcript_length(gtf_db, trans_id)

            # Add to the total read counts
            if new_trans_id not in total_read_counts:
                total_read_counts[new_trans_id] = 0
            total_read_counts[new_trans_id] += read_count

            # Add row to the dataframe
            exp_df = exp_df.append({
                'transcript_id': new_trans_id,
                'length': transcript_length,
                'count': total_read_counts[new_trans_id]
            }, ignore_index=True)

    # Calculate total mapped reads
    total_mapped_reads = exp_df['count'].sum()

    # Calculate FPKM and TPM
    exp_df = calculate_fpkm(exp_df, total_mapped_reads)
    exp_df = calculate_tpm(exp_df)

    # Output results
    print(exp_df)

    # Optionally save results to a CSV file
    # exp_df.to_csv('expression_results.csv', index=False)

# Command line argument handling
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate transcript expression levels.')
    parser.add_argument('--rda', dest='rda_files', nargs='+', required=True, help='The .rda files containing read counts.')
    parser.add_argument('--gtf', dest='gtf_files', nargs='+', required=True, help='The .gtf files containing gene annotations.')
    args = parser.parse_args()

    # Ensure the number of RDA and GTF files match
    if len(args.rda_files) != len(args.gtf_files):
        print("The number of RDA files must match the number of GTF files.")
        sys.exit(1)

    main(args.rda_files, args.gtf_files)