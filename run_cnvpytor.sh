#!/bin/bash
SCRIPT_START_TIME=$(date +%s)

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
    echo "  -o,  --output <file>       Output VCF filename (default: output.vcf)"
    echo "  -baf <true|false>          Whether to run BAF + combined call (default: false)"
    echo "  -mask_snp <fasta>          Optional. Input strict mask FASTA to generate CNVpytor-compatible mask file"
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
MASK_SNP=""
MASK_SNP_FLAG=""
FASTA_FILE=""
BED_FILE=""


# ===== parser =====
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--bin) BIN_SIZE="$2"; shift ;;
        -t|--thread) THREAD="$2"; shift ;;
        -r|--root) ROOT_FILE="$2"; shift ;;
        -bam) BAM_FILE="$2"; shift ;;
        -vcf) VCF_FILE="$2"; shift ;;
        -o|--output) OUTPUT_VCF="$2"; shift ;;
        -baf) USE_BAF="$2"; shift ;;
        -mask_snp) MASK_SNP_FLAG=1 MASK_SNP="$2"; shift ;;
        *) echo "unknown parameter: $1"; exit 1 ;;
    esac
    shift
done
OUTPUT_BAF_VCF="${OUTPUT_VCF%.vcf}_with_baf.vcf"

if [[ -n "$MASK_SNP_FLAG" ]]; then
    if [[ -z "$MASK_SNP" ]]; then
        echo " ERROR: -mask_snp option given but no file specified!"
        exit 1
    fi

    echo "▶ Detected mask SNP FASTA: $MASK_SNP"
    echo "▶ Generating CNVpytor-compatible mask file..."
    cnvpytor -root "$ROOT_FILE" -mask "$MASK_SNP" -make_mask_file
fi

# check argment
if [[ -z "$BAM_FILE" ]]; then
    echo "▶ CNV call must provide -bam BAM file"
    exit 1
fi

if [[ "$USE_BAF" == true ]] && [[ -z "$VCF_FILE" ]]; then
    echo "▶ BAF needs -vcf SNP file, generating via bcftools (it may take a long time)..."

    # 若沒有提供 reference fasta，就請使用者輸入
    if [[ -z "$FASTA_FILE" ]]; then
        echo -n " No reference FASTA provided. Please enter the path to the reference FASTA file: "
        read FASTA_FILE

        # 若使用者還是沒輸入，就退出
        if [[ -z "$FASTA_FILE" ]]; then
            echo "ERROR: You must provide a reference FASTA file to continue."
            exit 1
        fi
    fi

    # 詢問是否提供 BED 檔案
    echo -n "Optional: Enter BED file to limit regions (leave blank to use whole genome): "
    read BED_FILE

    # 組裝 mpileup 指令
    MPILEUP_CMD=(bcftools mpileup --threads "$THREAD" -f "$FASTA_FILE" -a AD,DP)

    if [[ -n "$BED_FILE" ]]; then
        MPILEUP_CMD+=(-R "$BED_FILE")
    fi

    MPILEUP_CMD+=(-Ou "$BAM_FILE")

    # 執行 SNP calling
    "${MPILEUP_CMD[@]}" | \
      bcftools call --threads "$THREAD" -mv -Oz -o snp.vcf.gz

    # 更新 VCF_FILE 變數
    VCF_FILE="snp.vcf.gz"

    # 建立 index
    bcftools index "$VCF_FILE"
fi


if [[ "$USE_BAF" == true ]]; then
    SAMPLE_NAME=$(bcftools query -l "$VCF_FILE" | head -n 1)
    echo "▶ Auto-detected sample: $SAMPLE_NAME"
fi


# 1. 載入 BAM 建立 RD
SECONDS=0
cnvpytor -root "$ROOT_FILE" -rd "$BAM_FILE" -j "$THREAD"
STEP1_TIME=$SECONDS

# 2. 建立 histogram（RD binning）
SECONDS=0
cnvpytor -root "$ROOT_FILE" -his "$BIN_SIZE" -j "$THREAD"
STEP2_TIME=$SECONDS

# 3. CNV 區段切分
SECONDS=0
cnvpytor -root "$ROOT_FILE" -partition "$BIN_SIZE" -j "$THREAD"
STEP3_TIME=$SECONDS

# 4. 初步 CNV calling（不含 BAF）
SECONDS=0
cnvpytor -root "$ROOT_FILE" -call "$BIN_SIZE" -j "$THREAD"
STEP4_TIME=$SECONDS

    cnvpytor -root "$ROOT_FILE" -view "$BIN_SIZE" <<EOF
set print_filename $OUTPUT_VCF
print calls
EOF

# ===== 是否執行 BAF 分析與 combined call =====
if [[ "$USE_BAF" == true ]]; then
    SECONDS=0
    cnvpytor -root "$ROOT_FILE" -snp "$VCF_FILE" -sample "$SAMPLE_NAME" -j "$THREAD" -nofilter
    STEP5_TIME=$SECONDS

    SECONDS=0
    cnvpytor -root "$ROOT_FILE" -mask_snps -j "$THREAD"
    STEP6_TIME=$SECONDS

    SECONDS=0
    cnvpytor -root "$ROOT_FILE" -baf "$BIN_SIZE" -j "$THREAD"
    STEP7_TIME=$SECONDS
    
    SECONDS=0
    cnvpytor -root "$ROOT_FILE" -call combined "$BIN_SIZE" -j "$THREAD"
    STEP8_TIME=$SECONDS

    cnvpytor -root "$ROOT_FILE" -view "$BIN_SIZE" <<EOF
set callers combined_mosaic
set print_filename $OUTPUT_BAF_VCF
print calls
EOF
fi

# ===== 時間統整輸出 =====
echo ""
echo "========== CNVpytor 執行時間總結 =========="
printf "%-35s %6s sec\n" "Step 1: Load BAM and extract RD"         "$STEP1_TIME"
printf "%-35s %6s sec\n" "Step 2: Create RD histogram (binning)"   "$STEP2_TIME"
printf "%-35s %6s sec\n" "Step 3: Partition RD segments"           "$STEP3_TIME"
printf "%-35s %6s sec\n" "Step 4: CNV Calling (RD only)"           "$STEP4_TIME"


TOTAL_TIME=$((STEP1_TIME + STEP2_TIME + STEP3_TIME + STEP4_TIME))

if [[ "$USE_BAF" == true ]]; then
    printf "%-35s %6s sec\n" "Step 5: Load SNPs from VCF"           "$STEP5_TIME"
    printf "%-35s %6s sec\n" "Step 6: Mask unreliable SNPs"        "$STEP6_TIME"
    printf "%-35s %6s sec\n" "Step 7: Calculate BAF"               "$STEP7_TIME"
    printf "%-35s %6s sec\n" "Step 8: Combined CNV calling"        "$STEP8_TIME"
    TOTAL_TIME=$((TOTAL_TIME + STEP5_TIME + STEP6_TIME + STEP7_TIME + STEP8_TIME))
fi

SCRIPT_END_TIME=$(date +%s)
ELAPSED_TIME=$((SCRIPT_END_TIME - SCRIPT_START_TIME))
ELAPSED_MIN=$((ELAPSED_TIME / 60))
ELAPSED_SEC=$((ELAPSED_TIME % 60))
TOTAL_MIN=$((TOTAL_TIME / 60))
TOTAL_SEC=$((TOTAL_TIME % 60))

echo "-------------------------------------------"
echo "Total measured steps time: ${TOTAL_MIN} min ${TOTAL_SEC} sec"
echo "Total wall-clock time: ${ELAPSED_MIN} min ${ELAPSED_SEC} sec"
echo "==========================================="

# you can add some filter command after "cnvpytor -root "$ROOT_FILE" -view "$BIN_SIZE" <<EOF"
# set Q0_range -1 0.5        # filter calls with more than half not uniquely mapped reads
# set p_range 0 0.0001       # filter non-confident calls 
# set p_N 0 0.5              # filter calls with more than 50% Ns in reference genome 
# set size_range 50000 inf   # filter calls smaller than 50kbp 
# set dG_range 100000 inf    # filter calls close to gaps in reference genome (<100kbp)
