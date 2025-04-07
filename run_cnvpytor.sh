#!/bin/bash

# 顯示使用說明
if [[ "$1" == "-h" || "$1" == "--help" || "$1" == "" ]]; then
    echo ""
    echo "CNVpytor analysis pipeline usage:"
    echo ""
    echo "  ./run_cnvpytor.sh [Options]"
    echo ""
    echo "Options:"
    echo "  -b,  --bin <int>           Bin size (default: 10000)"
    echo "  -t,  --thread <int>        Number of threads (default: 8)"
    echo "  -r,  --root <file>         Output path of .pytor file (default: file.pytor)"
    echo "  -bam <file>                Input BAM file"
    echo "  -vcf <file>                Input VCF file (required if -baf true)"
    echo "  -s,  --sample <name>       Sample name from VCF (required if -baf true, default: test_sample)"
    echo "  -o,  --output <file>       Output VCF filename (default: output.vcf)"
    echo "  -baf <true|false>          Whether to run BAF + combined call (default: false)"
    echo "  -h,  --help                Show this help message"
    echo ""
    exit 0
fi

# parameter
BIN_SIZE=10000
THREAD=8
ROOT_FILE="file.pytor"
BAM_FILE=""
VCF_FILE=""
SAMPLE_NAME="test_sample"
OUTPUT_VCF="output.vcf"
USE_BAF=false

# ===== parser =====
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bin) BIN_SIZE="$2"; shift ;;
        -t|--thread) THREAD="$2"; shift ;;
        -r|--root) ROOT_FILE="$2"; shift ;;
        -bam) BAM_FILE="$2"; shift ;;
        -vcf) VCF_FILE="$2"; shift ;;
        -s|--sample) SAMPLE_NAME="$2"; shift ;;
        -o|--output) OUTPUT_VCF="$2"; shift ;;
        -baf) USE_BAF="$2"; shift ;;
        *) echo "unknown parameter: $1"; exit 1 ;;
    esac
    shift
done
OUTPUT_BAF_VCF="${OUTPUT_VCF%.vcf}_with_baf.vcf"

# check argment
if [[ -z "$BAM_FILE" ]]; then
    echo "must provide -bam BAM file"
    exit 1
fi

if [[ "$USE_BAF" == true ]] && ([[ -z "$VCF_FILE" ]]); then
    echo "BAF nead -vcf snp file"
    exit 1
fi

# 1. 載入 BAM 建立 RD
cnvpytor -root "$ROOT_FILE" -rd "$BAM_FILE" -j "$THREAD"

# 2. 建立 histogram（RD binning）
cnvpytor -root "$ROOT_FILE" -his "$BIN_SIZE" -j "$THREAD"

# 3. CNV 區段切分
cnvpytor -root "$ROOT_FILE" -partition "$BIN_SIZE" -j "$THREAD"

# 4. 初步 CNV calling（不含 BAF）
cnvpytor -root "$ROOT_FILE" -call "$BIN_SIZE" -j "$THREAD"

    cnvpytor -root "$ROOT_FILE" -view "$BIN_SIZE" <<EOF
set print_filename $OUTPUT_VCF
print calls
EOF

# ===== 是否執行 BAF 分析與 combined call =====
if [[ "$USE_BAF" == true ]]; then
    cnvpytor -root "$ROOT_FILE" -snp "$VCF_FILE" -sample "$SAMPLE_NAME" -j "$THREAD"

    cnvpytor -root "$ROOT_FILE" -mask_snps -j "$THREAD"

    cnvpytor -root "$ROOT_FILE" -baf "$BIN_SIZE" -j "$THREAD"
    
    cnvpytor -root "$ROOT_FILE" -call combined "$BIN_SIZE" -j "$THREAD"

    cnvpytor -root "$ROOT_FILE" -view "$BIN_SIZE" <<EOF
set callers combined_mosaic
set print_filename $OUTPUT_BAF_VCF
print calls
EOF
fi
