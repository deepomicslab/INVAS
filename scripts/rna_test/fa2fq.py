import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def convert_fasta_to_fastq(fasta_file, fastq_file):
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')  # Remove the alphabet argument
    fastq_records = []
    for fasta in fasta_sequences:
        fastq_record = SeqRecord(Seq(str(fasta.seq)), id=fasta.id, description="")
        fastq_record.letter_annotations["phred_quality"] = [40] * len(fasta)
        fastq_records.append(fastq_record)
    SeqIO.write(fastq_records, fastq_file, "fastq")

def main():
    parser = argparse.ArgumentParser(description="Convert a FASTA file to a FASTQ file.")
    parser.add_argument("fasta_file", help="The name of the FASTA file.")
    parser.add_argument("fastq_file", help="The name of the output FASTQ file.")
    args = parser.parse_args()

    convert_fasta_to_fastq(args.fasta_file, args.fastq_file)

if __name__ == "__main__":
    main()
