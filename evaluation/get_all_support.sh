filename=$1
# 检查文件是否存在
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi
res_dir="../res_plus/"


# 使用循环读取文件的每一行，每行包含 WGS 和 RNA 的样本名
ALL_SUPPORT_FILE="$res_dir"/all_support.txt
while read -r sample_name wgs_sample1 rna_sample; do
    echo "precessing $sample_name,$wgs_sample1"
    sample_res1_dir="$res_dir"/$sample_name/
    python get_support.py -i $sample_res1_dir -m 1 -o $sample_res1_dir/support.txt
    if [ -s $ALL_SUPPORT_FILE ];then
	    echo "" >> $ALL_SUPPORT_FILE
	fi
    cat $sample_res1_dir/support.txt >> $ALL_SUPPORT_FILE

	
    
done < "$filename"



