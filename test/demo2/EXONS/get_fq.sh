# ===== 配置参数 =====
FA1="exon_1+2+3+4+5+6+9+.fa"; D1=80
FA2="exon_1+2+3+4-5+6+9+.fa";  D2=40
FA3="exon_1+2+5+6+9+.fa";  D3=40

READ_LEN=150
ERR=0
SEED=2025
R_SNV=0
R_INDEL=0
FRAC_INDEL_EXT=0

# 统一质量分数（Phred），推荐 40；如需 60，则改为 60
FIX_Q=60

# ===== 工具函数 =====
fa_len() { awk '/^>/ {next} {l+=length($0)} END{print l+0}' "$1"; }
calc_N() {
  L="$1"; DEPTH="$2"; RL="$3"
  awk -v L="$L" -v D="$DEPTH" -v RL="$RL" 'BEGIN{N=int((D*L)/(2.0*RL)+0.5); if(N<1)N=1; print N}'
}

# 将 FASTQ 的质量行替换为固定质量（Phred+33）
fix_qual() {
  in="$1"; out="$2"; q="$3"
  # 计算对应字符
  qchar=$(awk -v q="$q" 'BEGIN{printf "%c", q+33}')
  awk -v qchar="$qchar" '
    NR%4==2{seq=$0; print; next}
    NR%4==0{print gensub(/./, qchar, "g", substr(seq,1,length(seq))); next}
    {print}
  ' "$in" > "$out"
}

run_wgsim() {
  FA="$1"; DEPTH="$2"
  STEM="${FA##*/}"; STEM="${STEM%.*}"
  L=$(fa_len "$FA")
  N=$(calc_N "$L" "$DEPTH" "$READ_LEN")
  OUT1="${STEM}__${DEPTH}x_R1.fastq"
  OUT2="${STEM}__${DEPTH}x_R2.fastq"
  TMP1="${OUT1%.fastq}.tmp.fastq"
  TMP2="${OUT2%.fastq}.tmp.fastq"
  echo "[INFO] $FA len=$L depth=${DEPTH}x -> N_pairs=$N -> $OUT1 / $OUT2"

  # 先用 wgsim 生成临时文件
  wgsim -N "$N" -1 "$READ_LEN" -2 "$READ_LEN" \
        -e "$ERR" -r "$R_SNV" -R "$R_INDEL" -X "$FRAC_INDEL_EXT" \
        -S "$SEED" \
        "$FA" "$TMP1" "$TMP2"

  # 固定质量分数
  fix_qual "$TMP1" "$OUT1" "$FIX_Q"
  fix_qual "$TMP2" "$OUT2" "$FIX_Q"
  rm -f "$TMP1" "$TMP2"

  # 压缩
  gzip -f "$OUT1" "$OUT2"
  echo "[OK] $(basename "$OUT1").gz and $(basename "$OUT2").gz written with fixed Q=$FIX_Q"
}

# ===== 执行 =====
run_wgsim "$FA1" "$D1"
run_wgsim "$FA2" "$D2"
run_wgsim "$FA3" "$D3"