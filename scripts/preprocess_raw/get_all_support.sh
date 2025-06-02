filename=$1
# 检查文件是否存在
if [ ! -f "$filename" ]; then
    echo "File not found: $filename"
    exit 1
fi
res_dir="res/"
if [ ! -d "$res_dir" ]; then
    mkdir -p "$res_dir"
fi
# 使用循环读取文件的每一行, cat all the support files
ALL_SUPPORT_FILE="$res_dir"/all_support.txt
while read -r line; do
    echo "precessing $line"
    sample_res_dir="$res_dir"/$line
    python get_support.py -i $sample_res_dir -m 1 -o $sample_res_dir/support.txt
    cat $sample_res_dir/support.txt >> $ALL_SUPPORT_FILE
    
done < "$filename"
echo "All done"