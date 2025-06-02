import argparse
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def generate_exon_sequences(gtf_file, fasta_file, output_file, extend_length=0):
    # Load the GTF file using gffutils
    db = gffutils.create_db(gtf_file, ':memory:')

    # Load the FASTA file using BioPython
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Prepare a list to hold the new SeqRecord objects
    new_records = []

    # Iterate over each exon in the GTF file
    for gene in db.features_of_type('gene'):
            # print(gene)
        if gene.id == 'ENSG00000225293':
            for exon in list(db.children(gene, featuretype='exon')):
                # Calculate the new start and end positions
                exon_len = exon.end - exon.start
                extend_length = int(exon_len * 0.5)
                # start = max(0, exon.start - extend_length - 1)  # -1 because GFF is 1-based
                # end = exon.end + extend_length

                start = max(0, exon.start )  # -1 because GFF is 1-based
                end = exon.end
                

                # Extract the sequence from the FASTA file
                sequence = fasta_sequences[exon.seqid].seq[start:end]
                print(exon.seqid,exon_len, len(sequence))

                # Create a new SeqRecord object and add it to the list
                seq_id = "{}-{}-{}".format(exon.seqid, start, end)

                new_record = SeqRecord(sequence, id=seq_id, description='{}-{}-{}'.format(exon.seqid, exon.start, exon.end))
                new_records.append(new_record)

            # Write the new sequences to the output FASTA file
            SeqIO.write(new_records, output_file, "fasta")

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate exon sequences from a GTF and a reference FASTA file.')
parser.add_argument('--gtf', required=True, help='Input GTF file')
parser.add_argument('--fasta', required=True, help='Input FASTA file')
parser.add_argument('--output', required=True, help='Output FASTA file')
args = parser.parse_args()

# Call the function with command line arguments
generate_exon_sequences(args.gtf, args.fasta, args.output)