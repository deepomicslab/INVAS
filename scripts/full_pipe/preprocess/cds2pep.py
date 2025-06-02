import argparse

codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

# 函数将CDS序列转换为氨基酸序列
def translate_cds(cds):
    protein = ""
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, '')
            protein += amino_acid
    return protein

# 设置命令行参数
parser = argparse.ArgumentParser(description="Translate a CDS sequence into an amino acid sequence.")
parser.add_argument("cds", type=str, help="The CDS sequence to translate.")

# 解析命令行参数
args = parser.parse_args()

# 执行翻译
protein_sequence = translate_cds(args.cds.upper())

# 输出结果
print(f"Protein Sequence: {protein_sequence}")