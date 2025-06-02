import bisect
from Bio import SeqIO

def parse_fasta_header(header):
    # Assuming the header is in the format '>chr:start-end'
    chrom, start, end = header.split(" ")[1].split('-')
    return int(start), int(end)

def binary_search(headers, target):
    # Sort headers based on the start position
    headers.sort(key=parse_fasta_header)
    starts = [parse_fasta_header(header)[0] for header in headers]
    ends = [parse_fasta_header(header)[1] for header in headers]

    index = bisect.bisect(starts, target)

    # If the target location is within a sequence
    if index > 0 and target <= ends[index - 1]:
        prev_header = headers[index - 2] if index > 1 else None
        curr_header = headers[index - 1]
        next_header = headers[index] if index < len(headers) else None
    else:
        # If the target is not within a sequence, get the previous and next headers
        prev_header = headers[index - 1] if index > 0 else None
        curr_header = None
        next_header = headers[index] if index < len(headers) else None

    return prev_header, curr_header, next_header

fasta_dict = SeqIO.to_dict(SeqIO.parse("ENSG00000225293.fa", "fasta"))
headers = [fasta_dict[header].description for header in fasta_dict.keys()]
target = 16387693


prev_header, curr_header, next_header = binary_search(headers, target)

print(f'Previous sequence: {prev_header}')
print(f'Current sequence: {curr_header}')
print(f'Next sequence: {next_header}')

fasta_index = SeqIO.index("ENSG00000225293.fa", "fasta")
if prev_header:
    seq_record = fasta_index[prev_header.split(" ")[0]]
    print(f'Previous sequence content: {seq_record.seq}')
if curr_header:
    seq_record = fasta_index[curr_header.split(" ")[0]]
    print(f'Current sequence content: {seq_record.seq}')
if next_header:
    seq_record = fasta_index[next_header.split(" ")[0]]
    print(f'Next sequence content: {seq_record.seq}')