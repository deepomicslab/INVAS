import gffutils
import pyfaidx
import argparse


def get_seq(seq_segments, fasta):
    seq = ''
    for seqs in seq_segments:
        if seqs[0][0].strand == '+':
            seq += ''.join([s[1].seq.upper() for s in seqs])
        else:
            seq += ''.join([s[1].reverse.complement.seq.lower() for s in seqs[::-1]])
    return seq



def get_gene_sequences(gtf_file, genome_file, output_file):
    # 创建数据库
    db = gffutils.create_db(gtf_file, ':memory:', disable_infer_genes=True,
    disable_infer_transcripts=True)
    print(db)

    # 创建FASTA索引
    fasta = pyfaidx.Fasta(genome_file)
    print(db.features_of_type('gene'))
    # 打开输出文件
    with open(output_file, 'w') as out:
        # 遍历每个基因
        for gene in db.features_of_type('gene'):
            transcripts = list(db.children(gene, featuretype='transcript'))
            for transcript in transcripts:
                exons = list(db.children(transcript, featuretype='exon'))

                sequences = []  
                prev_exon_strand = None
                seq_segments = []
                sub_seq_segments = []
                for e in exons:
                    seq = fasta.get_seq(e.seqid, e.start, e.end)
                    # print(e.strand, prev_exon_strand, sub_seq_segments, seq_segments)
                    if e.strand != prev_exon_strand:
                        if prev_exon_strand == None:
                            sub_seq_segments.append((e,seq))
                        elif len(sub_seq_segments) > 0:
                            # print("?", e.start, e.strand)
                            seq_segments.append(sub_seq_segments)
                            sub_seq_segments = [(e,seq)]
                    else:
                        # print(e.start, e.strand)
                        sub_seq_segments.append((e, seq))
                    prev_exon_strand = e.strand
                seq_segments.append(sub_seq_segments)
                    
                sequence =  get_seq(seq_segments, fasta)
                print("seq_segments",seq_segments)
                out.write('>' + transcript.id + '\n') 
                out.write(sequence + '\n')
            
def main():
    parser = argparse.ArgumentParser(description='Generate gene sequences from GTF and genome FASTA file.')
    parser.add_argument('-g','--gtf', help='Input GTF file', required=True)
    parser.add_argument('-f','--fasta', help='Input genome FASTA file', required=True)
    parser.add_argument('-o','--output', help='Output file', required=True)
    args = parser.parse_args()

    get_gene_sequences(args.gtf, args.fasta, args.output)

if __name__ == "__main__":
    main()