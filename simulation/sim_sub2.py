#!/usr/bin/env python3

import argparse
import os
import random
import subprocess
import numpy as np
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import gffutils
import pybedtools
import shutil

class InversionSimulator:
    def __init__(self, ref_fasta, gtf_file, exclude_bed, output_dir, wgs_depth, rna_depth):
        self.ref_fasta = ref_fasta
        self.gtf_file = gtf_file
        self.exclude_bed = exclude_bed
        self.output_dir = output_dir
        self.wgs_depth = wgs_depth
        self.rna_depth = rna_depth
        
        # Create output directories
        self.wgs_dir = os.path.join(output_dir, "wgs")
        self.rna_dir = os.path.join(output_dir, "rna")
        self.truth_dir = os.path.join(output_dir, "truth")
        
        for directory in [self.wgs_dir, self.rna_dir, self.truth_dir]:
            os.makedirs(directory, exist_ok=True)
        
        # Load reference genome
        print("Loading reference genome...")
        self.ref_sequences = {record.id: record for record in SeqIO.parse(ref_fasta, "fasta")}
        
        # Create GTF database
        print("Creating GTF database...")
        if not os.path.exists(gtf_file + ".db"):
            self.db = gffutils.create_db(gtf_file, gtf_file + ".db", force=True, 
                                         merge_strategy='merge', disable_infer_genes=True, verbose=True)
        else:
            self.db = gffutils.FeatureDB(gtf_file + ".db")
        
        # Load exclude regions
        print("Loading exclude regions...")
        self.exclude_regions = pybedtools.BedTool(exclude_bed)
        
        # Initialize truth file
        self.truth_file = os.path.join(self.truth_dir, "inversions.tsv")
        with open(self.truth_file, 'w') as f:
            f.write("gene_id\tgene_name\tchromosome\tgene_start\tgene_end\tinversion_start\t"
                    "inversion_end\tinversion_type\tzygosity\ttranscript_id\ttranscript_fraction\t"
                    "is_inverted\taffected_exons\texon_positions\n")
        
        self.modified_genome = None
        self.modified_transcripts = {}
        self.gene_transcript_map = defaultdict(list)
        self.processed_genes = set()
        self.original_transcripts = {}
        
    def get_genes(self):
        """Get forward strand protein-coding genes that don't overlap with excluded regions"""
        genes = []
        for gene in self.db.features_of_type('gene'):
            if ('gene_type' in gene.attributes and 'protein_coding' in gene.attributes['gene_type'] 
                and gene.strand == '+'):  
                if gene.seqid not in self.ref_sequences:
                    continue
                    
                gene_bed = pybedtools.BedTool(f"{gene.seqid} {gene.start-1} {gene.end}", from_string=True)
                if not gene_bed.intersect(self.exclude_regions).count():
                    genes.append(gene)
        
        return genes
    
    def get_transcripts(self, gene):
        """Get all transcripts for a gene"""
        return list(self.db.children(gene, featuretype='transcript'))
    
    def get_exons(self, transcript):
        """Get all exons for a transcript in order"""
        exons = list(self.db.children(transcript, featuretype='exon'))
        exons.sort(key=lambda x: x.start)  
        return exons
    
    def determine_inversion_type(self):
        """随机选择反转类型"""
        inversion_types = [
            "internal_exons", 
            "single_internal_exon" 
        ]
        weights = [0.7, 0.3]  
        return random.choices(inversion_types, weights=weights)[0]
    
    def determine_zygosity(self):
        """Randomly determine zygosity"""
        return random.choice(["homozygous", "heterozygous"])
    
    def select_internal_region_for_inversion(self, gene, transcript):
        """选择基因内部区域进行反转，避免第一个和最后一个外显子"""
        exons = self.get_exons(transcript)
        if len(exons) < 3:  
            return None, None, None, None
        
        internal_exons = exons[1:-1]
        if not internal_exons:
            return None, None, None, None
        
        if len(internal_exons) == 1:
            selected_exons = internal_exons
            inversion_type = "single_internal_exon"
        else:
            inversion_type = self.determine_inversion_type()
            
            if inversion_type == "single_internal_exon":
                selected_exon_idx = random.randint(0, len(internal_exons) - 1)
                selected_exons = [internal_exons[selected_exon_idx]]
            else:  # internal_exons mode
                max_exons = min(3, len(internal_exons))
                num_exons = random.randint(2, max_exons)
                start_idx = random.randint(0, len(internal_exons) - num_exons)
                selected_exons = internal_exons[start_idx:start_idx + num_exons]
        
        inversion_start = selected_exons[0].start
        inversion_end = selected_exons[-1].end
        
        exon_ids = [exon.id for exon in selected_exons]
        print(f"  Selected {len(selected_exons)} internal exons for inversion: {exon_ids}")
        print(f"  Inversion region: {inversion_start}-{inversion_end}")
        
        return selected_exons, inversion_start, inversion_end, inversion_type
    
    def build_transcript_sequence(self, transcript):
        exons = self.get_exons(transcript)
        transcript_seq = ""
        chromosome = transcript.seqid
        ref_seq = str(self.ref_sequences[chromosome].seq)
        
        for exon in exons:
            transcript_seq += ref_seq[exon.start-1:exon.end]
        
        return transcript_seq
    
    def build_inverted_transcript_sequence(self, transcript, inversion_start, inversion_end):
        exons = self.get_exons(transcript)
        inverted_seq = ""
        chromosome = transcript.seqid
        ref_seq = str(self.ref_sequences[chromosome].seq)
        
        for exon in exons:
            if exon.end < inversion_start or exon.start > inversion_end:
                inverted_seq += ref_seq[exon.start-1:exon.end]
            elif exon.start >= inversion_start and exon.end <= inversion_end:
                exon_seq = ref_seq[exon.start-1:exon.end]
                inverted_exon_seq = str(Seq(exon_seq).reverse_complement())
                inverted_seq += inverted_exon_seq
            else:
                if exon.start < inversion_start:
                    inverted_seq += ref_seq[exon.start-1:inversion_start-1]
                    
                overlap_start = max(exon.start, inversion_start)
                overlap_end = min(exon.end, inversion_end)
                overlap_seq = ref_seq[overlap_start-1:overlap_end]
                inverted_overlap = str(Seq(overlap_seq).reverse_complement())
                inverted_seq += inverted_overlap
                
                if exon.end > inversion_end:
                    inverted_seq += ref_seq[inversion_end:exon.end]
        
        return inverted_seq
    
    def apply_inversion(self, sequence, start, end):
        """Apply inversion to a sequence"""
        start_idx = start - 1
        end_idx = end - 1
        
        inversion_seq = sequence[start_idx:end_idx+1]
        inverted_seq = str(Seq(inversion_seq).reverse_complement())
        
        new_sequence = sequence[:start_idx] + inverted_seq + sequence[end_idx+1:]
        return new_sequence
    
    def create_modified_genome_region(self, gene, inversion_start, inversion_end, zygosity):
        """Create modified genome region with inversion (gene ± 1000bp)"""
        chromosome = gene.seqid
        ref_seq = str(self.ref_sequences[chromosome].seq)
        
        chrom_length = len(ref_seq)
        region_start = max(1, min(gene.start, inversion_start) - 1000)
        region_end = min(chrom_length, max(gene.end, inversion_end) + 1000)
        
        original_region_seq = ref_seq[region_start-1:region_end]
        
        relative_inversion_start = inversion_start - region_start + 1
        relative_inversion_end = inversion_end - region_start + 1
        
        if relative_inversion_start < 1:
            relative_inversion_start = 1
        if relative_inversion_end > len(original_region_seq):
            relative_inversion_end = len(original_region_seq)
        
        modified_region_seq = self.apply_inversion(
            original_region_seq, 
            relative_inversion_start, 
            relative_inversion_end
        )
        
        if zygosity == "homozygous":
            self.modified_genome = {
                f"{chromosome}:{region_start}-{region_end}_hom_inv": {
                    "sequence": modified_region_seq,
                    "fraction": 1.0
                }
            }
        else:  # heterozygous
            self.modified_genome = {
                f"{chromosome}:{region_start}-{region_end}_het_ref": {
                    "sequence": original_region_seq,
                    "fraction": 0.5
                },
                f"{chromosome}:{region_start}-{region_end}_het_inv": {
                    "sequence": modified_region_seq,
                    "fraction": 0.5
                }
            }
        
        return {
            "chromosome": chromosome,
            "region_start": region_start,
            "region_end": region_end,
            "relative_inversion_start": relative_inversion_start,
            "relative_inversion_end": relative_inversion_end
        }
    
    def create_modified_transcripts(self, gene, transcripts, affected_exons, inversion_start, inversion_end, inversion_type, zygosity):
        """Create modified transcripts with inversion"""
        gene_id = gene.id
        self.modified_transcripts = {}
        self.original_transcripts = {} 
        
        total_transcripts = len(transcripts)
        if total_transcripts == 0:
            return
        
        transcript_weights = {}
        for transcript in transcripts:
            transcript_id = transcript.id
            transcript_weights[transcript_id] = np.random.lognormal(0, 1)
        
        total_weight = sum(transcript_weights.values())
        for tid in transcript_weights:
            transcript_weights[tid] /= total_weight
        
        for transcript in transcripts:
            transcript_id = transcript.id
            
            original_transcript_seq = self.build_transcript_sequence(transcript)
            
            self.original_transcripts[transcript_id] = original_transcript_seq
            
            exons = self.get_exons(transcript)
            transcript_contains_affected_exons = any(exon.id in [ae.id for ae in affected_exons] for exon in exons)
            
            
            if transcript_contains_affected_exons:
                if zygosity == "homozygous":
                    normal_fraction = transcript_weights[transcript_id] * 0.5
                    inv_fraction = transcript_weights[transcript_id] * 0.5
                    
                    self.modified_transcripts[transcript_id] = {
                        "original": original_transcript_seq,
                        "fraction": normal_fraction,
                        "gene_id": gene_id,
                        "is_inverted": False
                    }
                    
                    inverted_transcript_seq = self.build_inverted_transcript_sequence(
                        transcript, inversion_start, inversion_end)
                    
                    self.modified_transcripts[f"{transcript_id}_inv"] = {
                        "original": inverted_transcript_seq,
                        "fraction": inv_fraction,
                        "gene_id": gene_id,
                        "is_inverted": True
                    }
                else:  # heterozygous
                    normal_fraction = transcript_weights[transcript_id] * 0.7
                    inv_fraction = transcript_weights[transcript_id] * 0.3
                    
                    self.modified_transcripts[transcript_id] = {
                        "original": original_transcript_seq,
                        "fraction": normal_fraction,
                        "gene_id": gene_id,
                        "is_inverted": False
                    }
                    
                    inverted_transcript_seq = self.build_inverted_transcript_sequence(
                        transcript, inversion_start, inversion_end)
                    
                    self.modified_transcripts[f"{transcript_id}_inv"] = {
                        "original": inverted_transcript_seq,
                        "fraction": inv_fraction,
                        "gene_id": gene_id,
                        "is_inverted": True
                    }
            else:
                self.modified_transcripts[transcript_id] = {
                    "original": original_transcript_seq,
                    "fraction": transcript_weights[transcript_id],
                    "gene_id": gene_id,
                    "is_inverted": False
                }
        
        total_fraction = sum(info["fraction"] for info in self.modified_transcripts.values())
        
        if abs(total_fraction - 1.0) > 0.001:
            adjustment_factor = 1.0 / total_fraction
            for tid in self.modified_transcripts:
                self.modified_transcripts[tid]["fraction"] *= adjustment_factor
        
        for tid, info in self.modified_transcripts.items():
            self.gene_transcript_map[info["gene_id"]].append(tid)
    
    def write_fasta(self, sequences, filename):
        """Write sequences to a FASTA file"""
        with open(filename, 'w') as f:
            for seq_id, seq_info in sequences.items():
                if isinstance(seq_info, dict) and "sequence" in seq_info:
                    f.write(f">{seq_id}\n{seq_info['sequence']}\n")
                elif isinstance(seq_info, dict) and "original" in seq_info:
                    f.write(f">{seq_id}\n{seq_info['original']}\n")
                else:
                    f.write(f">{seq_id}\n{seq_info}\n")
    
    def write_truth_info(self, gene, affected_exons, inversion_start, inversion_end, inversion_type, zygosity):
        """Write truth information to file including both normal and inverted transcripts"""
        with open(self.truth_file, 'a') as f:
            gene_id = gene.id
            gene_name = gene.attributes.get('gene_name', [gene_id])[0]
            chromosome = gene.seqid
            gene_start = gene.start
            gene_end = gene.end
            
            exon_ids = []
            exon_positions = []
            
            for exon in affected_exons:
                exon_ids.append(exon.id)
                
                overlap_start = max(exon.start, inversion_start)
                overlap_end = min(exon.end, inversion_end)
                exon_positions.append(f"{exon.id}:{overlap_start}-{overlap_end}")
            
            for transcript_id, info in self.modified_transcripts.items():
                if info['gene_id'] != gene_id:
                    continue
                    
                expression_fraction = info['fraction']
                
                is_inverted = info['is_inverted']
                is_inverted_str = "Yes" if is_inverted else "No"
                
                if is_inverted:
                    affected_exons_field = ','.join(exon_ids)
                    exon_positions_field = ','.join(exon_positions)
                else:
                    affected_exons_field = ""
                    exon_positions_field = ""
                
                f.write(f"{gene_id}\t{gene_name}\t{chromosome}\t{gene_start}\t{gene_end}\t"
                        f"{inversion_start}\t{inversion_end}\t{inversion_type}\t{zygosity}\t"
                        f"{transcript_id}\t{expression_fraction:.6f}\t{is_inverted_str}\t"
                        f"{affected_exons_field}\t{exon_positions_field}\n")
    
    def simulate_wgs_reads_for_gene(self, gene_id):
        """Simulate WGS reads for a specific gene region using wgsim，"""
        if not self.modified_genome:
            return
        
        gene_wgs_dir = os.path.join(self.wgs_dir, gene_id)
        os.makedirs(gene_wgs_dir, exist_ok=True)
        
        reads_r1_files = []
        reads_r2_files = []
        
        for seq_id, seq_info in self.modified_genome.items():
            sequence = seq_info["sequence"]
            fraction = seq_info["fraction"]
            
            allele_fasta = os.path.join(gene_wgs_dir, f"{seq_id.replace(':', '_')}.fa")
            with open(allele_fasta, 'w') as f:
                f.write(f">{seq_id}\n{sequence}\n")
            
            region_size = len(sequence)
            allele_depth = self.wgs_depth * fraction
            read_pairs = int(allele_depth * region_size / 300)  # 300 = 2*150 (paired reads)
            
            r1_out = os.path.join(gene_wgs_dir, f"{seq_id.replace(':', '_')}_1.fq")
            r2_out = os.path.join(gene_wgs_dir, f"{seq_id.replace(':', '_')}_2.fq")
            
            reads_r1_files.append(r1_out)
            reads_r2_files.append(r2_out)
            
            wgsim_cmd = [
                "wgsim",
                "-e", "0",      # base error rate
                "-d", "500",        # outer distance between the two reads
                "-s", "50",         # standard deviation of fragment length
                "-N", str(read_pairs),  # number of read pairs
                "-1", "150",        # length of the first read
                "-2", "150",        # length of the second read
                "-r", "0",      # rate of mutations
                "-R", "0",      # rate of indels
                allele_fasta,
                r1_out,
                r2_out
            ]
            
            print(f"  Simulating WGS reads for gene {gene_id}, allele {seq_id} (size: {region_size}bp, depth: {allele_depth:.1f}x, reads: {read_pairs})")
            subprocess.run(wgsim_cmd, check=True)
        
        combined_r1 = os.path.join(gene_wgs_dir, "reads_1.fq")
        combined_r2 = os.path.join(gene_wgs_dir, "reads_2.fq")
        
        with open(combined_r1, 'w') as out_r1, open(combined_r2, 'w') as out_r2:
            for r1_file, r2_file in zip(reads_r1_files, reads_r2_files):
                if os.path.exists(r1_file) and os.path.exists(r2_file):
                    with open(r1_file, 'r') as in_r1, open(r2_file, 'r') as in_r2:
                        out_r1.write(in_r1.read())
                        out_r2.write(in_r2.read())
        
        return (combined_r1, combined_r2)
    
    def simulate_rna_reads_for_gene(self, gene_id):
        """Simulate RNA-seq reads for a specific gene using wgsim"""
        if gene_id not in self.gene_transcript_map:
            return
        
        transcript_ids = self.gene_transcript_map[gene_id]
        if not transcript_ids:
            return
            
        gene_rna_dir = os.path.join(self.rna_dir, "genes", gene_id)
        os.makedirs(gene_rna_dir, exist_ok=True)
        
        for transcript_id in transcript_ids:
            if transcript_id not in self.modified_transcripts:
                continue
                
            transcript_info = self.modified_transcripts[transcript_id]
            transcript_seq = transcript_info["original"]
            expression_fraction = transcript_info["fraction"]
            
            transcript_fasta = os.path.join(gene_rna_dir, f"{transcript_id}.fa")
            with open(transcript_fasta, 'w') as f:
                f.write(f">{transcript_id}\n{transcript_seq}\n")
            
            transcript_length = len(transcript_seq)
            read_pairs = int(self.rna_depth * transcript_length * expression_fraction / 300)
            
            if expression_fraction > 0 and read_pairs < 10:
                read_pairs = 10
            
            if expression_fraction == 0 or read_pairs == 0:
                continue
            
            r1_out = os.path.join(gene_rna_dir, f"{transcript_id}_1.fq")
            r2_out = os.path.join(gene_rna_dir, f"{transcript_id}_2.fq")
            
            wgsim_cmd = [
                "wgsim",
                "-e", "0",      # base error rate
                "-d", "350",        # outer distance between the two reads
                "-s", "30",         # standard deviation of fragment length
                "-N", str(read_pairs),  # number of read pairs
                "-1", "150",        # length of the first read
                "-2", "150",        # length of the second read
                "-r", "0",       # rate of mutations
                "-R", "0",       # rate of indels
                transcript_fasta,
                r1_out,
                r2_out
            ]
            
            print(f"  Simulating RNA reads for transcript {transcript_id} (length: {transcript_length}bp, reads: {read_pairs})")
            subprocess.run(wgsim_cmd, check=True)
            
        gene_r1 = os.path.join(gene_rna_dir, f"{gene_id}_1.fq")
        gene_r2 = os.path.join(gene_rna_dir, f"{gene_id}_2.fq")
        
        with open(gene_r1, 'w') as out_r1, open(gene_r2, 'w') as out_r2:
            for transcript_id in transcript_ids:
                r1_file = os.path.join(gene_rna_dir, f"{transcript_id}_1.fq")
                r2_file = os.path.join(gene_rna_dir, f"{transcript_id}_2.fq")
                
                if os.path.exists(r1_file) and os.path.exists(r2_file):
                    with open(r1_file, 'r') as in_r1, open(r2_file, 'r') as in_r2:
                        out_r1.write(in_r1.read())
                        out_r2.write(in_r2.read())
        
        return (gene_r1, gene_r2)
    
    def merge_all_reads(self):
        """Merge all gene-level reads into master files"""
        master_rna_r1 = os.path.join(self.rna_dir, "rna_reads_1.fq")
        master_rna_r2 = os.path.join(self.rna_dir, "rna_reads_2.fq")
        
        master_wgs_r1 = os.path.join(self.wgs_dir, "wgs_reads_1.fq")
        master_wgs_r2 = os.path.join(self.wgs_dir, "wgs_reads_2.fq")
        
        print("Merging RNA-seq reads...")
        with open(master_rna_r1, 'w') as out_r1, open(master_rna_r2, 'w') as out_r2:
            for gene_id in self.gene_transcript_map:
                gene_rna_dir = os.path.join(self.rna_dir, "genes", gene_id)
                r1_file = os.path.join(gene_rna_dir, f"{gene_id}_1.fq")
                r2_file = os.path.join(gene_rna_dir, f"{gene_id}_2.fq")
                
                if os.path.exists(r1_file) and os.path.exists(r2_file):
                    with open(r1_file, 'r') as in_r1, open(r2_file, 'r') as in_r2:
                        out_r1.write(in_r1.read())
                        out_r2.write(in_r2.read())
        
        print("Merging WGS reads...")
        with open(master_wgs_r1, 'w') as out_r1, open(master_wgs_r2, 'w') as out_r2:
            for gene_id in self.processed_genes:
                gene_wgs_dir = os.path.join(self.wgs_dir, gene_id)
                r1_file = os.path.join(gene_wgs_dir, "reads_1.fq")
                r2_file = os.path.join(gene_wgs_dir, "reads_2.fq")
                
                if os.path.exists(r1_file) and os.path.exists(r2_file):
                    with open(r1_file, 'r') as in_r1, open(r2_file, 'r') as in_r2:
                        out_r1.write(in_r1.read())
                        out_r2.write(in_r2.read())
        
        print("All reads merged successfully.")
    
    def is_gene_qualified(self, gene):
        if gene.id in self.processed_genes:
            return False
            
        if gene.strand != '+':
            return False
            
        gene_length = gene.end - gene.start + 1
        if gene_length > 10000:
            return False
            
        if 'gene_type' not in gene.attributes or 'protein_coding' not in gene.attributes['gene_type']:
            return False
            
        transcripts = self.get_transcripts(gene)
        if not transcripts:
            return False
            
        has_sufficient_exons = False
        for transcript in transcripts:
            if len(self.get_exons(transcript)) >= 3:
                has_sufficient_exons = True
                break
                
        return has_sufficient_exons
    
    def get_qualified_genes(self, all_genes, min_required):
        qualified_genes = []
        random.shuffle(all_genes) 
        
        for gene in all_genes:
            if self.is_gene_qualified(gene):
                qualified_genes.append(gene)
                if len(qualified_genes) >= min_required:
                    break
                    
        return qualified_genes
    
    def simulate_inversions(self, num_genes=50):
        """Simulate gene inversions for a number of randomly selected forward strand protein-coding genes"""
        print(f"Selecting {num_genes} forward strand protein-coding genes for inversion simulation...")
        
        all_genes = self.get_genes()
        print(f"Found {len(all_genes)} forward strand protein-coding genes.")
        
        qualified_genes = self.get_qualified_genes(all_genes, num_genes * 2)
        
        print(f"Initial screening found {len(qualified_genes)} qualified forward strand protein-coding genes with ≥3 exons and length ≤ 10kb.")
        
        successful_genes = 0
        remaining_attempts = len(qualified_genes)
        
        while successful_genes < num_genes and remaining_attempts > 0:
            gene = qualified_genes.pop(0)
            remaining_attempts -= 1
            
            gene_id = gene.id
            gene_length = gene.end - gene.start + 1
            print(f"Processing gene {successful_genes+1}/{num_genes}: {gene_id} (length: {gene_length}bp, remaining attempts: {remaining_attempts})")
            
            transcripts = self.get_transcripts(gene)
            
            if not transcripts:
                print(f"  Skipping gene {gene_id} - no transcripts found")
                continue
                
            reference_transcript = random.choice(transcripts)
            
            zygosity = self.determine_zygosity()
            zygosity = "heterozygous"
            
            affected_exons, inversion_start, inversion_end, inversion_type = self.select_internal_region_for_inversion(gene, reference_transcript)
            
            if not affected_exons or inversion_start is None or inversion_end is None:
                print(f"  Skipping gene {gene_id} - couldn't select internal region for inversion")
                
                if len(qualified_genes) < (num_genes - successful_genes) * 2:
                    more_genes = self.get_qualified_genes(all_genes, (num_genes - successful_genes) * 2)
                    for g in more_genes:
                        if g.id not in [x.id for x in qualified_genes]:
                            qualified_genes.append(g)
                    print(f"  Added {len(more_genes)} more qualified genes, now have {len(qualified_genes)} candidates.")
                    
                continue  
            
            region_info = self.create_modified_genome_region(gene, inversion_start, inversion_end, zygosity)
            
            self.create_modified_transcripts(gene, transcripts, affected_exons, 
                                        inversion_start, inversion_end, inversion_type, zygosity)
            
            self.write_truth_info(gene, affected_exons, inversion_start, inversion_end, 
                                inversion_type, zygosity)
            
            self.simulate_wgs_reads_for_gene(gene_id)
            
            self.simulate_rna_reads_for_gene(gene_id)
            
            self.processed_genes.add(gene_id)
            successful_genes += 1
            
            print(f"  Successfully processed gene {gene_id}. Progress: {successful_genes}/{num_genes}")
            
            if len(qualified_genes) < (num_genes - successful_genes) * 2 and (num_genes - successful_genes) > 0:
                more_genes = self.get_qualified_genes(all_genes, (num_genes - successful_genes) * 2)
                for g in more_genes:
                    if g.id not in [x.id for x in qualified_genes]:
                        qualified_genes.append(g)
                print(f"  Added {len(more_genes)} more qualified genes, now have {len(qualified_genes)} candidates.")
        
        if successful_genes < num_genes:
            print(f"Warning: Only managed to successfully simulate inversions for {successful_genes}/{num_genes} genes.")
        
        self.merge_all_reads()
        
        print(f"Simulation completed successfully with {successful_genes} genes!")
        
def main():
    parser = argparse.ArgumentParser(description='Simulate gene inversions in forward strand protein-coding genes')
    parser.add_argument('--ref', required=True, help='Reference genome FASTA file (hg38.fa)')
    parser.add_argument('--gtf', required=True, help='GTF annotation file (hg38.gtf)')
    parser.add_argument('--exclude', required=True, help='BED file with regions to exclude (hg38.exclude.bed)')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--wgs_depth', type=float, default=30.0, help='WGS sequencing depth')
    parser.add_argument('--rna_depth', type=float, default=50.0, help='RNA-seq sequencing depth')
    parser.add_argument('--num_genes', type=int, default=50, help='Number of genes to simulate inversions for')
    
    args = parser.parse_args()
    
    simulator = InversionSimulator(
        args.ref, 
        args.gtf, 
        args.exclude, 
        args.outdir,
        args.wgs_depth,
        args.rna_depth
    )
    
    simulator.simulate_inversions(args.num_genes)

if __name__ == "__main__":
    main()