import pandas as pd
import pyreadr
import gffutils

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

# 主函数
def main(csv_file):
    # 加载数据
    # expression_data = load_data(csv_file)
    rows = []
    for trans_id, read_count in trans_read_counts.items():
        # Get the length of the transcript
        transcript_length = get_transcript_length(gtf_db, trans_id)
        # 将数据添加到列表中
        rows.append({'transcript_id': trans_id, 'length': transcript_length, 'count': read_count})

    # 一次性将列表转换成DataFrame
    exp_df = pd.DataFrame(rows, columns=['transcript_id', 'length', 'count'])



    # 计算总的mapped reads数
    total_mapped_reads = exp_df['count'].sum()

    # 计算FPKM
    exp_df = calculate_fpkm(exp_df, total_mapped_reads)

    # 计算TPM
    exp_df = calculate_tpm(exp_df)

    # 输出结果
    print(exp_df)
    
    # # 保存结果到新的CSV文件
    # expression_data.to_csv('expression_results.csv', index=False)

# 调用主函数
if __name__ == "__main__":
    import sys


    
    if len(sys.argv) < 3:
        print("Usage: python calculate_expression.py <read_count.rda> <truth.gtf>")
        sys.exit(1)
    # result = pyreadr.read_r('/run/media/wangxuedong/One Touch/simulation/5_trans_2inv/1_trans_inv1/sim_counts_matrix.rda')
    result = pyreadr.read_r(sys.argv[1])
    counts_matrix_df = result['counts_matrix']
    trans_read_counts = counts_matrix_df['sample_01'].to_dict()
    # Create a database from the GTF file
    gtf_db = create_db(sys.argv[2])

    
    
    main(sys.argv[2])

# transcript_id,length,count
# transcript1,1500,300
# transcript2,2000,150