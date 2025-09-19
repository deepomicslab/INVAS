conda activate invas_data_prep
sample_name="test"
# Create output directory
mkdir -p rna_out/$sample_name
threads=2

# Run RNA processing
bash ../scripts/full_pipe/preprocess/process_rna_common.sh -n $sample_name \
    -b $input_bam \
    -o rna_out/$sample_name \
    -r $reference_genome \
    -x $hisat_index \
    -t $threads \
    -p scripts/full_pipe/bin/picard.jar