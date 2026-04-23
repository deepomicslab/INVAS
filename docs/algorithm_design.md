# INVAS Algorithm Design Choices

This document explains the rationale behind key default parameter and tool selections in INVAS.
It also describes the exact parameters used in the SSW-based back-splicing junction detection algorithm.

---

## Table of Contents

1. [Why HISAT2 as the Default RNA-seq Aligner?](#1-why-hisat2-as-the-default-rna-seq-aligner)
2. [Why NCBI RefSeq as the Default Annotation Database?](#2-why-ncbi-refseq-as-the-default-annotation-database)
3. [SSW-Based Back-Splicing Junction Detection: Parameters and Rationale](#3-ssw-based-back-splicing-junction-detection-parameters-and-rationale)

---

## 1. Why HISAT2 as the Default RNA-seq Aligner?

INVAS uses [HISAT2](https://daehwankimlab.github.io/hisat2/) as the default RNA-seq aligner for the following reasons:

### Memory Efficiency
Alternative splice-aware aligners such as [STAR](https://github.com/alexdobin/STAR) can require **>30 GB of RAM** just for genome index loading, which may be prohibitive in clinical settings or on standard desktop workstations. HISAT2 maintains a **low memory footprint** (~4–8 GB for a human genome index), making INVAS accessible to a broader range of researchers and institutional environments.

### Splice-Junction Awareness
HISAT2 uses a **hierarchical graph FM index (HGFM)**, combining a global FM index of the genome with thousands of local FM indexes of genomic regions. This architecture makes it particularly robust for aligning reads **across annotated and novel splice sites** — a critical requirement for detecting structural variant breakpoints such as intragenic inversions where reads bridge non-canonical junctions.

### Compatibility with Downstream Processing
The HISAT2 output BAM format integrates seamlessly with the rest of the INVAS preprocessing pipeline (Picard duplicate marking, BWA-based remapping of unmapped reads, SSW-based soft-clip analysis).

### Alternative Aligners
Users are **not restricted to HISAT2**. The `-x` parameter accepts any HISAT2-format index path. If you have already aligned your RNA-seq data with another aligner (e.g., STAR, STAR-Fusion), you can provide the resulting BAM file directly to the downstream preprocessing scripts. INVAS only requires a position-sorted, indexed BAM as input to `process_rna_common.sh`.

---

## 2. Why NCBI RefSeq as the Default Annotation Database?

INVAS is pre-configured to use [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) annotation files for the following reasons:

### Curated, Non-Redundant Transcript Set
RefSeq focuses on a **curated, non-redundant set of transcripts** with experimental or functional validation. In contrast, databases like Ensembl or GENCODE contain a much larger volume of **computationally predicted isoforms** that lack direct experimental support.

### Clean Splice Graph Construction
INVAS constructs a **conjugate alternative splice graph (CASG)** from reference annotations as a baseline for *de novo* inversion detection. Using a highly curated database:
- Simplifies the splice graph structure.
- Reduces false-positive inversion calls that might arise from aligning reads to poorly defined or predicted transcripts.
- Establishes a clean reference against which truly novel inverted junctions can be identified.

### Ready-to-Use Files Provided
To facilitate a **quick start**, pre-formatted RefSeq annotation BED files for **hg19** and **hg38** are directly provided in the repository:

```
scripts/full_pipe/annotation/hg19_genes.gencode.bed
scripts/full_pipe/annotation/hg38_genes.gencode.bed
```

This allows users to run the full pipeline immediately **without manual database downloading or format conversion**.

### Using Custom Annotations (GENCODE / Ensembl)
Users who prefer GENCODE or Ensembl annotations can supply a custom GTF/BED file.  
Simply format it to match the expected BED schema and pass it via the `$gene_annotation_bed` parameter (see [README.md](../README.md) for usage).

---

## 3. SSW-Based Back-Splicing Junction Detection: Parameters and Rationale

INVAS uses the **SIMD Striped Smith–Waterman (SSW) library** to resolve soft-clipped read tails at candidate inversion breakpoints. This section documents the exact parameters and explains why SSW outperforms seed-based aligners for this task.

### Alignment Scoring Parameters

| Parameter | Value |
|-----------|-------|
| Match score | +2 |
| Mismatch penalty | −2 |
| Gap open penalty | 3 |
| Gap extension penalty | 1 |
| Masking length | ⌊read\_len / 2⌋ |
| Orientation | Both forward and reverse-complement; best is kept |

The **masking length** (⌊read\_len / 2⌋) suppresses spurious secondary alignment peaks in repetitive regions.

### Filtering Thresholds

Alignments are retained only if they satisfy **all** of the following:

| Criterion | Threshold |
|-----------|-----------|
| Minimum aligned length | ≥ 12 nt |
| Match proportion (identity) | ≥ 0.85 (i.e., ≥ 85% identity after gap penalties) |
| Maximum gaps (aligned length < 25 nt) | ≤ 1 gap |
| Maximum gaps (aligned length ≥ 25 nt) | ≤ 2 gaps |
| Ambiguity filter (best vs. suboptimal score) | (nScore − nScore2) / nScore ≥ 0.10 |

When a candidate breakpoint falls within **blacklisted or segmental-duplication regions**, filters are tightened:

| Criterion | Tightened Threshold |
|-----------|---------------------|
| Minimum aligned length | ≥ 15 nt |
| Match proportion | ≥ 0.90 (≥ 90% identity) |

Additionally, whenever both 5′ and 3′ clipped evidence are present, they must be **strand- and position-consistent**.

### Why SSW Instead of Seed-Based Aligners?

Soft-clipped segments near inversion junctions are frequently **very short (5–20 nt)**. We evaluated BWA-MEM (including supplementary alignments), BWA ultrashort/short presets, minimap2, STAR, Bowtie2, BLAST, and HISAT2 on simulated intragenic-inversion transcripts (*n* = 100) across overall, 5′-only, and 3′-only clip scenarios.

| Aligner | Performance on clips ≥ 25 nt | Performance on clips 5–20 nt |
|---------|------------------------------|-------------------------------|
| BWA-MEM (supplementary) | Good | Degrades; inconsistent within ±5 bp of breakpoint |
| minimap2 / STAR / HISAT2 | Adequate | Seed-length limited; lower sensitivity |
| Bowtie2 | Adequate | Underperforms on ultrashort seeds |
| BLAST | Sensitive | Less practical; harder to constrain junction placement |
| **SSW (windowed)** | **Excellent** | **Highest hit rate; most accurate breakpoint localisation** |

SSW is applied in a **constrained window** around each WGS-supported breakpoint (typically ±200 bp, refined to ±50 bp), so it does not rely on seeds and **optimises true local alignment**. This yields the highest hit rate and most accurate breakpoint localisation (within ±5 bp) for 5–20 nt clips — precisely the length range that defines inverted splice junctions.

### Why the SSW Library Specifically?

Our implementation wraps the [SSW C library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library), which extends Farrar's SIMD algorithm to return:
- The **optimal local alignment score**
- The **full CIGAR/path and precise alignment endpoints**
- A **heuristic suboptimal (second-best) score and location** in linear space with negligible extra cost

These outputs are critical:
- Precise CIGAR and endpoints enable junction pinpointing within ±5 bp.
- The suboptimal score enables the ambiguity filter `(nScore − nScore2) / nScore ≥ 0.10`, analogous to mapping-quality estimation.

Compared with other SIMD SW implementations (e.g., SSEARCH), SSW provides these features at substantially better performance in short-read genomic contexts.

---

## References

- Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. *Nat Methods.* 2015;12(4):357–360.
- Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics.* 2013;29(1):15–21.
- O'Leary NA, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res.* 2016;44(D1):D733–D745.
- Zhao M, et al. SSW Library: An SIMD Smith-Waterman C/C++ Library for Use in Genomic Applications. *PLOS ONE.* 2013;8(12):e82138.
