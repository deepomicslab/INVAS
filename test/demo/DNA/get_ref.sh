#!/usr/bin/env bash
set -euo pipefail

REF=../../../../wgs_ref/hg19.fa
CHR=chr12
GENE_START=71533490
GENE_END=71533627
FLANK=300

# 计算三段坐标
WIN_START=$((GENE_START - FLANK))
WIN_END=$((GENE_END + FLANK))
L_START=${WIN_START}
L_END=$((GENE_START - 1))
M_START=${GENE_START}
M_END=${GENE_END}
R_START=$((GENE_END + 1))
R_END=${WIN_END}

# 建索引
[ -f "${REF}.fai" ] || samtools faidx "${REF}"

# 提取三段
samtools faidx "${REF}" "${CHR}:${L_START}-${L_END}" > L.fa
samtools faidx "${REF}" "${CHR}:${M_START}-${M_END}" > M.fa
samtools faidx "${REF}" "${CHR}:${R_START}-${R_END}" > R.fa

# 辅助：取纯序列一行
seq1() { grep -v '^>' "$1" | tr -d '\n' ; }

# 反向互补（awk）
revcomp() {
  awk 'BEGIN{FS=OFS=""}
    function rc(s,  i,c,b,out){
      for(i=length(s); i>=1; i--){
        c=substr(s,i,1)
        if(c=="A") b="T"; else if(c=="C") b="G";
        else if(c=="G") b="C"; else if(c=="T") b="A";
        else if(c=="a") b="t"; else if(c=="c") b="g";
        else if(c=="g") b="c"; else if(c=="t") b="a";
        else b=c;
        out=out b
      } return out
    }
    {print rc($0)}'
}

# 拼接 L + rc(M) + R
LSEQ=$(seq1 L.fa)
MSEQ=$(seq1 M.fa | revcomp)
RSEQ=$(seq1 R.fa)

# 输出：原始窗口（不换行）
echo ">${CHR}_${WIN_START}_${WIN_END}_REF" > TSPAN8.hg19.flank300.ref.fa
samtools faidx "${REF}" "${CHR}:${WIN_START}-${WIN_END}" | grep -v '^>' | tr -d '\n' >> TSPAN8.hg19.flank300.ref.fa
echo >> TSPAN8.hg19.flank300.ref.fa

# 输出：带 inversion 的窗口（不换行）
OUT=TSPAN8.hg19.flank300.inv.fa
echo ">${CHR}_${WIN_START}_${WIN_END}_INV_${M_START}_${M_END}" > "${OUT}"
echo -n "${LSEQ}${MSEQ}${RSEQ}" >> "${OUT}"
echo >> "${OUT}"

# 简要校验：两序列长度应一致
echo -n "ref length: "
grep -v '^>' TSPAN8.hg19.flank300.ref.fa | wc -c
echo -n "inv length: "
grep -v '^>' TSPAN8.hg19.flank300.inv.fa | wc -c

# 清理
rm -f L.fa M.fa R.fa
