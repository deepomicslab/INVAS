#!/bin/bash

set -e

case "$1" in
    "step1"|"rna")
        shift
        export PATH="/opt/conda/envs/invas_data_prep/bin:$PATH"
        bash scripts/full_pipe/preprocess/process_rna_common.sh "$@"
        ;;
    "step2"|"sv")
        shift
        export PATH="/opt/conda/envs/invas_assembly/bin:$PATH"
        bash scripts/full_pipe/preprocess/process_sv_common.sh "$@"
        ;;
    "step3"|"candidate")
        shift
        export PATH="/opt/conda/envs/invas_assembly/bin:$PATH"
        bash scripts/full_pipe/preprocess/combine_sv_rna.sh "$@"
        ;;
    "step4"|"assembly")
        shift
        export PATH="/opt/conda/envs/invas_assembly/bin:$PATH"
        python scripts/full_pipe/main.py "$@"
        ;;
    "bash"|"shell")
        exec /bin/bash
        ;;
    *)
        echo "INVAS Docker Container"
        echo ""
        echo "Usage: docker run invas <command> [options]"
        echo ""
        echo "Commands:"
        echo "  step1, rna        - RNA-Seq data processing (invas_data_prep env)"
        echo "  step2, sv         - WGS SV detection (invas_assembly env)"
        echo "  step3, candidate  - Candidate inversion detection"
        echo "  step4, assembly   - Core transcriptome assembly"
        echo "  bash, shell       - Interactive shell"
        ;;
esac