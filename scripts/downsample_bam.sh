INPUT_FILE=$1
N_TARGET=$(($2 * 1000000))
SAMPLE_NAME=$3
OUTPUT_DIR=cell_line_subs

echo "[INFO] Computing number of reads..."
N_READS=$(samtools view -c $INPUT_FILE)
echo "Total number of reads: $N_READS"
echo "[INFO] Subsampling bam-file..."
SAMPLE_FRAC=$(echo "$N_TARGET / $N_READS" | bc -l)
samtools view -s $SAMPLE_FRAC -b $INPUT_FILE > $OUTPUT_DIR/$SAMPLE_NAME-$2.bam
