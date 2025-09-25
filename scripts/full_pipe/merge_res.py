import argparse
from collections import defaultdict

def parse_gtf(file_path, skip_header_lines=0):
    """
    parse a GTF file and return a dictionary mapping transcript_id
    """
    transcripts = defaultdict(lambda: {"transcript_line": None, "exon_lines": []})

    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#") or skip_header_lines > 0:  # skip header lines
                skip_header_lines -= 1
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # ignore invalid lines
            
            feature_type = fields[2]  # "transcript" or "exon"
            attributes = fields[8]

            # extract transcript_id
            transcript_id = extract_attribute(attributes, "transcript_id")
            print(transcript_id)
            if not transcript_id:
                continue  # skip if no transcript_id found
            
            # store the line based on feature type
            if feature_type == "transcript":
                transcripts[transcript_id]["transcript_line"] = line.strip().replace("StringTie", "Seqflow")
            elif feature_type == "exon":
                transcripts[transcript_id]["exon_lines"].append(line.strip().replace("StringTie", "Seqflow"))
    
    return transcripts

def extract_attribute(attributes, key):
    """
    extract the specified key (e.g., transcript_id) from the attribute field of a GTF file.
    """
    for attr in attributes.split(";"):
        attr = attr.strip()
        if attr.startswith(key):
            return attr.split('"')[1]
    return None

def merge_gtf(file1, file2, output_file, skip_header_file1=1, skip_header_file2=2):
    """
    Merge two GTF files, ensuring the integrity of each transcript (including transcript and exon).
    """
    # Parse the two GTF files
    transcripts1 = parse_gtf(file1, skip_header_lines=skip_header_file1)
    transcripts2 = parse_gtf(file2, skip_header_lines=skip_header_file2)

    # Merge the data
    merged_transcripts = {}

    # add the first file's content
    merged_transcripts.update(transcripts1)

    
    # add the second file's content, overwriting if transcript_id already exists
    merged_transcripts.update(transcripts2)
    print(merged_transcripts)

    # format exon line , with start < end and filter exon length < 2
    for transcript_id, transcript_data in merged_transcripts.items():
        for i in range(len(transcript_data["exon_lines"])):
            fields = transcript_data["exon_lines"][i].strip().split("\t")
            start = int(fields[3])
            end = int(fields[4])
            if start > end:
                fields[3], fields[4] = fields[4], fields[3]
            if end - start < 2:
                transcript_data["exon_lines"][i] = None
            else:
                transcript_data["exon_lines"][i] = "\t".join(fields)
        transcript_data["exon_lines"] = list(filter(None, transcript_data["exon_lines"]))
        
    
    # write to output file
    with open(output_file, "w") as out_f:
        for transcript_id, transcript_data in merged_transcripts.items():
            # write the transcript line first
            if transcript_data["transcript_line"]:
                out_f.write(transcript_data["transcript_line"] + "\n")
            # write all associated exon lines
            for exon_line in transcript_data["exon_lines"]:
                out_f.write(exon_line + "\n")
    
    print(f"Merged GTF file has been saved to: {output_file}")

# main function to handle command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two GTF files and retain complete transcripts (transcript + exon lines).")
    parser.add_argument("file1", type=str, help="Path to the first GTF file (1-line header).")
    parser.add_argument("file2", type=str, help="Path to the second GTF file (2-line header).")
    parser.add_argument("output_file", type=str, help="Path to save the merged GTF file.")
    args = parser.parse_args()
    
    merge_gtf(args.file1, args.file2, args.output_file)