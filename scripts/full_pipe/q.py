#!/usr/bin/env python3
"""
Simple BAM to Kallisto quantification script
Extract reads (including unmapped) from one BAM file and perform Kallisto quantification
"""

import pysam
import gzip
import subprocess
import argparse
from pathlib import Path

def extract_reads_from_bam(bam_path, fq1_path, fq2_path, include_unmapped=True):
    """Extract paired reads from BAM file to FASTQ files (including unmapped reads)"""
    print(f"Extracting reads from {bam_path}...")
    
    # Open BAM file
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    
    # Store paired reads
    read_pairs = {}
    unmapped_count = 0
    mapped_count = 0
    
    # Read all reads
    for read in bam:
        qname = read.query_name
        
        # Count mapped/unmapped
        if read.is_unmapped:
            unmapped_count += 1
        else:
            mapped_count += 1
        
        # Skip unmapped reads if not including them
        if not include_unmapped and read.is_unmapped:
            continue
            
        if qname not in read_pairs:
            read_pairs[qname] = {}
        
        if read.is_read1:
            read_pairs[qname]['R1'] = read
        elif read.is_read2:
            read_pairs[qname]['R2'] = read
    
    bam.close()
    
    print(f"  Total mapped reads: {mapped_count}")
    print(f"  Total unmapped reads: {unmapped_count}")
    
    # Write to FASTQ files
    count = 0
    single_count = 0
    
    with gzip.open(fq1_path, 'wt') as f1, gzip.open(fq2_path, 'wt') as f2:
        for qname, reads in read_pairs.items():
            if 'R1' in reads and 'R2' in reads:
                # Write read1
                r1 = reads['R1']
                if r1.query_qualities:
                    qual1 = ''.join([chr(q + 33) for q in r1.query_qualities])
                else:
                    qual1 = 'I' * len(r1.query_sequence)  # Default quality if missing
                f1.write(f"@{qname}\n{r1.query_sequence}\n+\n{qual1}\n")
                
                # Write read2
                r2 = reads['R2']
                if r2.query_qualities:
                    qual2 = ''.join([chr(q + 33) for q in r2.query_qualities])
                else:
                    qual2 = 'I' * len(r2.query_sequence)  # Default quality if missing
                f2.write(f"@{qname}\n{r2.query_sequence}\n+\n{qual2}\n")
                
                count += 1
            else:
                single_count += 1
    
    print(f"  Extracted {count} read pairs")
    print(f"  Skipped {single_count} unpaired reads")
    return count

def run_kallisto(ref_fa, fq1, fq2, output_dir, threads=8):
    """Run Kallisto quantification"""
    print("Running Kallisto...")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    kallisto_idx = output_dir / "kallisto.idx"
    kallisto_out = output_dir / "kallisto_quant"
    kallisto_out.mkdir(exist_ok=True)
    
    # Build index
    print("  Building Kallisto index...")
    subprocess.run([
        'kallisto', 'index',
        '-i', str(kallisto_idx),
        str(ref_fa)
    ], check=True, capture_output=True)
    
    # Run quantification
    print("  Running quantification...")
    subprocess.run([
        'kallisto', 'quant',
        '-i', str(kallisto_idx),
        '-o', str(kallisto_out),
        '-t', str(threads),
        '--bootstrap-samples', '100',
        str(fq1), str(fq2)
    ], check=True, capture_output=True)
    
    print(f"  Kallisto quantification complete. Results in: {kallisto_out}")

def main():
    parser = argparse.ArgumentParser(description='Extract reads from BAM file and run Kallisto quantification')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--ref', required=True, help='Reference FASTA file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('--skip-unmapped', action='store_true', 
                       help='Skip unmapped reads (default: include unmapped reads)')
    
    args = parser.parse_args()
    
    # Create path objects
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    # Set file paths (no temp directory needed for single BAM)
    fq1 = output_dir / "reads_R1.fastq.gz"
    fq2 = output_dir / "reads_R2.fastq.gz"
    
    try:
        # Step 1: Extract reads from BAM (including unmapped)
        print("\nStep 1: Extracting reads from BAM")
        include_unmapped = not args.skip_unmapped
        extract_reads_from_bam(args.bam, fq1, fq2, include_unmapped)
        
        # Step 2: Run Kallisto
        print("\nStep 2: Running Kallisto quantification")
        run_kallisto(args.ref, fq1, fq2, output_dir, args.threads)
        
        print("\nPipeline completed successfully!")
        
    except Exception as e:
        print(f"\nError: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())