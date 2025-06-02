import os
import argparse
import subprocess
import re
from utils import *


def run_command(command):
    """
    Run a shell command and handle errors.

    Args:
        command (str): Command to execute.
    """
    try:
        print(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error while running command: {command}")
        print(e)
        exit(1)

def parse_unique_genes_from_file(file_path):
    """
    Parse the file and extract unique gene names from the second column.

    Parameters:
        file_path (str): Path to the input file.

    Returns:
        set: A set of unique gene names.
    """
    unique_genes = set()  # To store unique gene names
    gene_region_dict = {}  # To store gene regions
    chrom=""

    with open(file_path, "r") as file:
        for line in file:
            # Split each line into columns by tab
            columns = line.strip().split("\t")
            if len(columns) > 1:
                # Extract the second column (e.g., 'PTPRF:43990857-44089343')
                gene_info = columns[1]
                # Extract the gene name before the colon
                gene_name = gene_info.split(":")[0]
                unique_genes.add(gene_name)
                # Extract the gene region
                # split gene_info by "," and "-" and ":" in the same time
                gene_part, region_part = gene_info.split(":")
                print(gene_part, region_part)
                region=region_part.split(",")[0]
                gene_start, gene_end = region.split("-")
                # res=re.split(r'[,:-]', gene_info)
                # gene_start, gene_end = res[1], res[2]
                chrom = columns[0].split(":")[0]
                gene_region_dict[gene_name] = (chrom, gene_start, gene_end)

    return unique_genes, gene_region_dict

def check_graph_file(file_path):
    """
    Check if the graph file is too short.

    Args:
        file_path (str): Path to the graph file.

    Returns:
        bool: True if the file is long enough, False otherwise.
    """
    segs=[]
    with open(file_path, "r") as file:
        for line in file:
                segs.append(line.strip())
    if len(segs) <=5:
        return False
    return True

def check_stringtie_gtf(file_path):
    """
    Check if the StringTie GTF file contains too few isoforms.

    Args:
        file_path (str): Path to the StringTie GTF file.

    Returns:
        bool: True if the file contains too few isoforms, False otherwise.
    """
    with open(file_path, "r") as file:
        count = 0
        # skip the first two header
        for line in file:
            if line.startswith("#"):
                continue
            items = line.strip().split("\t")
            if items[2] == "transcript":
                count += 1
    if count < 1:
        return True
    return False


def process_single_gene(
    gene, input_dir, output_dir, ref_genome, rna_map_bam, rna_unmap_bam, gene_dict
):
    """
    Process a single gene through the entire pipeline.

    Args:
        gene (str): Gene name.
        input_dir (str): Input directory.
        output_dir (str): Output directory.
        ref_genome (str): Reference genome file.
        rna_map_bam (str): RNA mapped BAM file.
        rna_unmap_bam (str): RNA unmapped BAM file.
    """
    gene_dir=os.path.join(output_dir, gene)
    os.makedirs(gene_dir, exist_ok=True)
    gene_region = gene_dict[gene]
    gene_region_str = f"{gene_region[0]}_{gene_region[1]}_{gene_region[2]}"
    gene_map_with_remap_bam = os.path.join(gene_dir, f"{gene_region_str}_hisat.map_unmapremap.s.bam")
    gene_unmap_bwa_bam = os.path.join(gene_dir, f"{gene_region_str}_still_unmap_bwa.s.bam")
    single_candidate_gene_file = os.path.join(gene_dir, "gene.txt")
    # Define file paths dynamically based on the gene name
    gene_rna_gtf = os.path.join(gene_dir, f"{gene}_stringtie_rna.gtf")
    inv_bed = os.path.join(gene_dir, f"{gene}_inv.bed")
    inversion_exons = os.path.join(gene_dir, f"{gene}_inversion.exons")
    final_inv_bed = os.path.join(gene_dir, f"{gene}_final_inv.bed")
    final_inverson_exons = os.path.join(gene_dir, f"{gene}_final_inv.exons")
    recover_bam = os.path.join(gene_dir, f"{gene}_recover.bam")
    sorted_recover_bam = os.path.join(gene_dir, f"{gene}_recover.sorted.bam")
    gtf2seggtrf_dir=os.path.join(gene_dir, "gtf2seggtrf")
    os.makedirs(gtf2seggtrf_dir, exist_ok=True)
    seg_file = f"{gene}_seg.o"
    trf_file = f"{gene}_trf.o"
    split_gtf = os.path.join(gene_dir, f"{gene}_split.gtf")
    merged_gtf = os.path.join(gene_dir, f"{gene}_merged.gtf")
    graph_file = os.path.join(gene_dir, f"{gene}_graph.g")
    frg_file=os.path.join(gene_dir, f"{gene}_frg.f")
    inv_junc_file=os.path.join(gene_dir, f"{gene}_inv.junc")
    run_flow_dir=os.path.join(gene_dir, "run_flow")
    os.makedirs(run_flow_dir, exist_ok=True)
    local_stringtie_gtf=os.path.join(gene_dir, f"{gene}_stringtie_local.gtf")
    flow_gtf=os.path.join(gene_dir, f"{gene}_flow.gtf")
    chromosome=gene_dict[gene][0]

    # Run StringTie to process the RNA BAM file
    run_command(f"""stringtie -i "{gene_map_with_remap_bam}" -o "{gene_rna_gtf}" """)

    # Check if the gene_rna_gtf file contains too few isoforms

    if check_stringtie_gtf(gene_rna_gtf):
        print(f"Warning: {gene_rna_gtf} contains too few isoforms.")
        return

    # Run consensus to parse real inversion regions
    run_command(
        f"""python "{script_dir}"/consensus_inv_region.py --bam_file {gene_map_with_remap_bam} --gene_file {single_candidate_gene_file} --output_file {inv_bed} --inv_exon_file {inversion_exons}"""
    )
    
    # get final inversion bed and exons
    run_command(
        f"""python "{script_dir}"/handle_excp.py {inv_bed} {final_inv_bed} {final_inverson_exons}"""
    )

    # Handle BAM files and recover the BAM
    run_command(
        f"""python "{script_dir}"/handle_inv_gene3.py -b1 {gene_unmap_bwa_bam} -b2 {gene_map_with_remap_bam} -inv {final_inv_bed} -outbam {recover_bam} -ref {ref_genome}"""
    )
    # sort and index the recover bam
    run_command(f"""samtools sort -o {sorted_recover_bam} {recover_bam}""")
    run_command(f"""samtools index {sorted_recover_bam}""")

    # Convert GTF to segment and TRF, and add inversion to segment
    run_command(
        f"""python "{script_dir}"/gtf2segtrf.py {gene_rna_gtf} {seg_file} {trf_file} {final_inverson_exons} {gtf2seggtrf_dir}"""
    )

    if not os.path.exists(f"{gtf2seggtrf_dir}/{seg_file}"):
        print(f"Warning: {gtf2seggtrf_dir} not build.")

        return
    elif not check_graph_file(f"{gtf2seggtrf_dir}/{seg_file}"):
        print(f"Warning: {gtf2seggtrf_dir}/{seg_file} too short.")
        return
    # Split GTF
    run_command(
        f"""python "{script_dir}"/split_ori_gtf_seg.py -g {gene_rna_gtf} -i "{gtf2seggtrf_dir}"/{seg_file} -o {split_gtf}"""
    )

    # Link segments
    run_command(
        f"""python "{script_dir}"/link_inv_frag4.py {sorted_recover_bam} {split_gtf} "{gtf2seggtrf_dir}"/{seg_file} {final_inverson_exons} "{gtf2seggtrf_dir}"/{trf_file} {graph_file} {frg_file} {inv_junc_file}"""
    )

    # check if graph file
    print(check_graphg_file(graph_file))
    if not check_graphg_file(graph_file):
        print(f"Warning: {graph_file} not pass check.")
        return
    # Run seqFlow (adjust this command as necessary for your environment)
    try:
        run_command(
            f""" "{script_dir}"/bin/matching -b \
                --model 1 -v 1 \
                -g {graph_file} \
                -r {run_flow_dir}/{gene}.path \
                -c {run_flow_dir}/test.c.path \
                -m {run_flow_dir}/test.m.g \
                --break_c -v 2 -a \
                -t {frg_file} """
        )
    except:
        print(f"Error while running seqFlow for {gene}.")
        return
    # run_command(
    #     f""" "{script_dir}"/bin/matching -b \
    #         --model 1 -v 1 \
    #         -g {graph_file} \
    #         -r {run_flow_dir}/{gene}.path \
    #         -c {run_flow_dir}/test.c.path \
    #         -m {run_flow_dir}/test.m.g \
    #         --break_c -v 2 -a \
    #         -t {frg_file} \
    #         --rna """
    # )

    # Run StringTie locally to get debug BAM status
    run_command(
        f""" "{script_dir}"/bin/stringtie {gene_map_with_remap_bam} \
        -o {local_stringtie_gtf} -v --debug --tmpdir {gene_dir} """
    )

    # Calculate expression
    run_command(
        f"""python "{script_dir}"/generate_exp.py {run_flow_dir}/{gene}.path {gene_dir}/debug_bam_status.txt {sorted_recover_bam} {flow_gtf} {chromosome}"""
    )

    # Merge split GTF with seqFlow GTF
    run_command(f"""python "{script_dir}"/merge_res.py {flow_gtf} {split_gtf} {merged_gtf}""")

    # Call SNVs, indels, and phase VCF
    run_command(f"""python "{script_dir}"/run_vars.py --bam "{args.wgs_bam}" --ref "{ref_genome}" --chrom {chromosome} --start {gene_region[1]} --end {gene_region[2]} --output "{gene_dir}"/wgs""")

    # if no reads in region_bam, continue to next gene, using samtools view -c in region_bam
    g_r=f"{chromosome}:{gene_region[1]}-{gene_region[2]}"
    if int(subprocess.check_output(["samtools", "view", "-c", f"{gene_dir}/wgs.region.bam", g_r])) < 10:
        print(f"No reads found in region {g_r}. Exiting...")
        return

    # Generate the final GTF file
    if args.with_normal_haps: 
        print("Using normal haps.")
        print(
            f"""python "{script_dir}"/generate_final_gtf.py --inversion_exons {final_inverson_exons} --vcf "{gene_dir}"/wgs.phased.vcf.gz --ref {ref_genome} -o "{gene_dir}"/haps -g {merged_gtf} -j {inv_junc_file} --with_normal_haps """
        )
        run_command(
            f"""python "{script_dir}"/generate_final_gtf.py --inversion_exons {final_inverson_exons} --vcf "{gene_dir}"/wgs.phased.vcf.gz --ref {ref_genome} -o "{gene_dir}"/haps -g {merged_gtf} -j {inv_junc_file} --with_normal_haps """
        )
    else:
        run_command(
            f"""python "{script_dir}"/generate_final_gtf.py --inversion_exons {final_inverson_exons} --vcf "{gene_dir}"/wgs.phased.vcf.gz --ref {ref_genome} -o "{gene_dir}"/haps -g {merged_gtf} -j {inv_junc_file} """
        )


def main():
    """
    Main function to orchestrate the pipeline.
    """
    

    # Split BAM by candidate genes

    run_command(
        f"""python "{script_dir}"/get_loci_bam.py -i "{candidate_gene_file}" -b1 "{args.rna_map_with_remap_bam}" -b2 "{args.rna_unmap_bwa_bam}" -o "{args.output_dir}" -e {args.extend_region} """
    )

    # Process each gene
    for gene in candidate_genes:
        # if gene.lower() != "lbhd1":
        #     continue
        # if error, continue to next gene
        try:
            process_single_gene(
                gene,
                args.input_dir,
                args.output_dir,
                args.ref_genome,
                args.rna_map_with_remap_bam,
                args.rna_unmap_bwa_bam,
                gene_dict
            )
        except Exception as e:
            print(f"Error while processing gene {gene}: {e}")
            continue


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Pipeline for processing genes.")
    parser.add_argument("--input_dir", required=True, help="Input directory containing gene data.")
    parser.add_argument("--sample_name", required=True, help="Sample name.")
    parser.add_argument("--extend_region", type=int, default=1000, help="Extension region size.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--ref_genome", required=True, help="Reference genome file.")
    parser.add_argument("--rna_map_with_remap_bam", required=True, help="RNA mapped and remapped BAM file.")
    parser.add_argument("--rna_unmap_bwa_bam", required=True, help="RNA unmapped BWA BAM file.")
    parser.add_argument("--wgs_bam", required=True, help="wgs BAM file.")
    parser.add_argument("--with_normal_haps", required=False, default=False, action="store_true", help="Whether to use normal haps.")
    args = parser.parse_args()


    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Get support tools for candidate genes
    candidate_gene_file = os.path.join(args.output_dir, "candidate_gene.txt")
    run_command(f"""python "{script_dir}"/get_candidate_vote.py -i "{args.input_dir}" -m 1 -o "{candidate_gene_file}" """)

    # Parse candidate gene files
    candidate_genes, gene_dict = parse_unique_genes_from_file(candidate_gene_file)
    main()