filename="$1"

while read -r sample_name wgs_sample1 rna_sample; do
    echo "Processing $sample_name with WGS sample $wgs_sample1 and RNA sample $rna_sample"
    gene_file=../$sample_name/truth/inversions.tsv
    # 去掉第一行得到第一列的unique gene list
    awk 'NR>1 {print $1}' $gene_file | sort | uniq > $gene_file.unique
    # 逐行读取gene list
    while read -r gene; do
        echo "Processing gene $gene"
        # 逐行读取gene list中的每个gene
        # 如果../$sample_name/rna/genes/$gene/analysis/stringtie/和../$sample_name/rna/genes/$gene/analysis/stringtie/$gene.stringtie.fasta都存在 则运行
        # stringite
        if [ -d ../$sample_name/rna/genes/$gene/analysis/stringtie/ ] && [ -f ../$sample_name/rna/genes/$gene/analysis/stringtie/$gene.stringtie.fasta ]; then
            python eva_seqs_single2.py \
            --ref ../$sample_name/truth_seq/$gene/transcripts.fa \
            --asm ../$sample_name/rna/genes/$gene/analysis/stringtie/$gene.stringtie.fasta \
            --out ../$sample_name/eva/$gene/analysis/stringtie/seq_eva2/ \
            --threads 30
        else
            echo "Skipping gene $gene in stringtie"
        fi
        # cufflinks
        if [ -d ../$sample_name/rna/genes/$gene/analysis/cufflinks/ ] && [ -f ../$sample_name/rna/genes/$gene/analysis/cufflinks/transcripts.gtf ]; then
            python eva_seqs_single2.py \
            --ref ../$sample_name/truth_seq/$gene/transcripts.fa \
            --asm ../$sample_name/rna/genes/$gene/analysis/cufflinks/$gene.cufflinks.fasta \
            --out ../$sample_name/eva/$gene/analysis/cufflinks/seq_eva2/ \
            --threads 30
        else
            echo "Skipping gene $gene in cufflinks"
        fi
        # trinity
        if [ -d ../$sample_name/rna/genes/$gene/analysis/trinity/ ] && [ -f ../$sample_name/rna/genes/$gene/analysis/trinity/Trinity.fasta ]; then
            python eva_seqs_single2.py \
            --ref ../$sample_name/truth_seq/$gene/transcripts.fa \
            --asm ../$sample_name/rna/genes/$gene/analysis/trinity/Trinity.fasta \
            --out ../$sample_name/eva/$gene/analysis/trinity/seq_eva2/ \
            --threads 30
        else
            echo "Skipping gene $gene in trinity"
        fi
        
        
    done < "$gene_file.unique"
done < "$filename"