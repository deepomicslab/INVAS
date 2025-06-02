import os
import subprocess
import argparse


def index_bam(bam_file):
    """
    Ensure that a BAM file is indexed.
    """
    if not os.path.exists(f"{bam_file}.bai"):
        print(f"Indexing BAM file: {bam_file}...")
        subprocess.run(["samtools", "index", bam_file], check=True)


def index_reference(ref_fasta):
    """
    Ensure that the reference FASTA file is indexed.
    """
    if not os.path.exists(f"{ref_fasta}.fai"):
        print(f"Indexing reference FASTA file: {ref_fasta}...")
        subprocess.run(["samtools", "faidx", ref_fasta], check=True)


def call_variants_freebayes(region_bam, ref_fasta, region, output_prefix):
    """
    Call variants using FreeBayes and normalize them with vcfallelicprimitives.
    """
    calls_vcf = f"{output_prefix}.calls.vcf"

    print("Calling variants with FreeBayes and normalizing with vcfallelicprimitives...")
    freebayes_cmd = [
        "freebayes",
        "-f", ref_fasta,
        "-r", region,
        region_bam
    ]
    vcfallelicprimitives_cmd = [
        "vcfallelicprimitives", "-kg"
    ]

    # Use a shell pipeline to connect FreeBayes and vcfallelicprimitives
    with open(calls_vcf, "w") as vcf_out:
        freebayes_proc = subprocess.Popen(freebayes_cmd, stdout=subprocess.PIPE)
        vcfallelic_proc = subprocess.Popen(vcfallelicprimitives_cmd, stdin=freebayes_proc.stdout, stdout=vcf_out)
        vcfallelic_proc.wait()
        freebayes_proc.wait()

    print(f"Variants called and normalized: {calls_vcf}")
    return calls_vcf


def phase_variants(calls_vcf, region_bam, ref_fasta, output_prefix):
    """
    Phase variants using WhatsHap.
    """
    phased_vcf = f"{output_prefix}.phased.vcf.gz"

    print("Phasing variants with WhatsHap...")
    whatshap_cmd = [
        "whatshap", "phase",
        "--reference", ref_fasta,
        "-o", phased_vcf,
        calls_vcf, region_bam
    ]
    subprocess.run(whatshap_cmd, check=True)
    # tabix the phased VCF
    subprocess.run(["tabix","-f", phased_vcf], check=True)

    print("Phasing completed.")
    return phased_vcf


def extract_region(bam_file, ref_fasta, chrom, start, end, extend, output_prefix):
    """
    Extract reads for a specific region, call variants, and phase them.
    """
    # Extend the region
    extended_start = max(1, start - extend)  # Ensure start is not less than 1
    extended_end = end + extend
    region = f"{chrom}:{extended_start}-{extended_end}"

    # Output file paths
    region_bam = f"{output_prefix}.region.bam"

    # Ensure BAM and reference are indexed
    index_bam(bam_file)
    index_reference(ref_fasta)

    # Extract reads for the region using `-o` option
    print(f"Extracting reads for region {region} into {region_bam}...")
    samtools_view_cmd = [
        "samtools", "view", "-b", "-o", region_bam, bam_file, region
    ]
    subprocess.run(samtools_view_cmd, check=True)

    # Index the region-specific BAM
    print("Indexing extracted BAM file...")
    subprocess.run(["samtools", "index", region_bam], check=True)

    # Check if no reads in region_bam use samtools view -c
    if int(subprocess.check_output(["samtools", "view", "-c", region_bam])) < 10:
        print(f"No reads found in region {region}. Exiting...")
        return None


    # Call variants using FreeBayes and normalize them
    calls_vcf = call_variants_freebayes(region_bam, ref_fasta, region, output_prefix)

    # Phase variants
    phased_vcf = phase_variants(calls_vcf, region_bam, ref_fasta, output_prefix)

    print("Pipeline completed successfully.")
    return phased_vcf


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Pipeline to call and phase variants from a WGS BAM file given a genomic region."
    )
    parser.add_argument("--bam", required=True, help="Path to the input WGS BAM file.")
    parser.add_argument("--ref", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g., 'chr1').")
    parser.add_argument("--start", type=int, required=True, help="Start coordinate of the region.")
    parser.add_argument("--end", type=int, required=True, help="End coordinate of the region.")
    parser.add_argument("--extend", type=int, default=0, help="Number of bases to extend the region on both sides (default: 0).")
    parser.add_argument("--output", required=True, help="Prefix for output files.")
    return parser.parse_args()


def main():
    # Parse command-line arguments
    args = parse_args()

    # Run the pipeline
    phased_vcf = extract_region(
        bam_file=args.bam,
        ref_fasta=args.ref,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        extend=args.extend,
        output_prefix=args.output
    )
    print(f"Phased VCF file: {phased_vcf}")


if __name__ == "__main__":
    main()