#!/bin/bash

# Path to the reference genome
REF="path-to/hg38.fa"
GTF="hg38.gtf"
EXCLUDE="custom_sv_exclusion_list.bed"
NUM_GENES=10000

# Set WGS depth range (30-100, step 20)
for wgs_depth in $(seq 30 20 110); do
    # Set RNA depth range (100-1000, step 150)
    for rna_depth in $(seq 150 150 1050); do
        # Create output directory name
        outdir="sim_${wgs_depth}_${rna_depth}"
        
        echo "Running simulation with WGS depth: ${wgs_depth}, RNA depth: ${rna_depth}"
        
        # Run the sim_sub2.py command
        python ../sim_sub2.py \
            --ref "${REF}" \
            --gtf "${GTF}" \
            --exclude "${EXCLUDE}" \
            --outdir "${outdir}" \
            --wgs_depth "${wgs_depth}" \
            --rna_depth "${rna_depth}" \
            --num_genes "${NUM_GENES}"
        
        # Check if the command executed successfully
        if [ $? -eq 0 ]; then
            echo "Simulation completed successfully for ${outdir}"
        else
            echo "Error: Simulation failed for ${outdir}"
        fi
    done
done

echo "All simulations completed!"