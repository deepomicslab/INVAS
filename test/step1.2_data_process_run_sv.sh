conda activate /data6/wangxuedong/invas_nc/invas_data_prep
reference_genome="/data6/wangxuedong/invas_nc/wgs_ref/hg19.fa"
sample_name="test"
# Create output directory
mkdir -p sv_out/$sample_name
threads=2

input_wgs_bam="/data6/wangxuedong/invas_nc/INVAS/test/demo/DNA/inv_wgs.hg19.bam"
# Run SV detection
bash ../scripts/full_pipe/preprocess/process_sv_common.sh -s $sample_name \
    -i $input_wgs_bam \
    -o sv_out/$sample_name \
    -r $reference_genome \
    -t $threads \
    -c delly,manta,svaba
