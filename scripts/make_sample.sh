#!/bin/bash

for arg in "$@"; do
   case "$arg" in
      pid=*) PID="${arg#*=}" ;;
      bed=*) BED="${arg#*=}" ;;
   esac
done

if [ -z "$PID" ]; then
        PID=$$
fi

echo $PID

if [ -z $BED ]; then
	echo "Generating random bed file for simulation"
	Rscript ../simulate_bed.R $PID # the file will be saved as simulationPID.bed                                   │
        BED_FILE=bed_files2/simulation$PID.bed
else
	BED_FILE=$BED
	echo "Using $4 as bed file for simulation"	
fi

REFERENCE=$1

#echo "Running simulation from $BED_FILE" # here $1 is the bed file to simulate from
echo "Using reference file $REFERENCE"

OUTPUT_DIR="output/simulation_$PID"
OUTPUT_BAM="simulation_$PID"

docker run --rm -v /home/lottaer/cnv-detection/simulated_data/:/my_data nabavilab/cnv-sim \
	./cnv-sim.py -o /my_data/$OUTPUT_DIR --read_length 100 --n_reads 0 --cnv_list /my_data/$BED_FILE genome /my_data/$1 

CANCER=bam_files/$OUTPUT_BAM
CONTROL=bam_files/{$OUTPUT_BAM}_control.bam

N_BPs=$((249 * 1000000)) # number of basepairs in the first chromosome
READ_LEN=100
N_TARGET=$(echo "$N_BPs * $2 / $READ_LEN" | bc -l)

echo "N_BPs = $N_BPs"
echo "Read length = $READ_LEN"
echo "N_TARGET = $N_TARGET"

# do the alignment ( only do cnv_1 and control 1 to start with)
../bwa-0.7.15/bwa mem $1 $OUTPUT_DIR/cnv_1.fastq | samtools sort -o $CANCER
../bwa-0.7.15/bwa mem $1 $OUTPUT_DIR/control_1.fastq | samtools sort -o $CONTROL

PURITY_INT=$3
PURITY=$( echo $PURITY_INT / 100 | bc -l )
# N_TARGET=$(($2 * 1000000)) # million of reads
echo "[INFO] Mixing sample of purity $PURTIY and $N_TARGET reads (coverage $2)"

LOWER_LIMIT=0.5
PURITY_CNV=$( Rscript ../purity_script.R -i $CANCER -b $BED_FILE)
echo "THE ESTIMATED PURITY IS $PURITY_CNV"

if awk "BEGIN {exit !($PURITY_CNV < $LOWER_LIMIT)}"; then
	rm $CANCER
	rm $CONTROL
	echo "Exiting: Too low purity estimated"
	exit 1
fi

echo "[INFO] Computing size of $CANCER"
# subsample the file with alterations
N_CANCER=$(samtools view -c $CANCER)
FACTOR_CANCER=$( echo "$N_TARGET / $N_CANCER" | bc -l )
WEIGHT_CANCER=$( echo "$FACTOR_CANCER * $PURITY / $PURITY_CNV" | bc -l )

echo "[INFO] Subsampling $CANCER"
samtools view -s $WEIGHT_CANCER -b $CANCER > bam_files/sub_$PID.bam

METADATA_FILE=metadata_eval/metadata_$PID.txt

touch $METADATA_FILE
N_CANCER_SUB=$(samtools view -c bam_files/sub_$PID.bam)
echo "$N_CANCER_SUB\n" > $METADATA_FILE

# repeat for the control sample
echo "[INFO] Computing size of $CONTROL"
N_CONTROL=$(samtools view -c $CONTROL)
N_SAMPLE_CONTROL=$( echo "$N_TARGET - $N_TARGET * $PURITY / $PURITY_CNV" | bc -l )
FACTOR_CONTROL=$( echo "$N_SAMPLE_CONTROL / $N_CONTROL" | bc -l )

echo "Factor control = $FACTOR_CONTROL"

echo "[INFO] Subsampling control"
samtools view -s $FACTOR_CONTROL -b $CONTROL > bam_files/control_sub_$PID.bam
N_CONTROL_SUB=$(samtools view -c bam_files/control_sub_$PID.bam)

echo "$N_CONTROL_SUB\n" >> $METADATA_FILE

CNV_WITH_EXTENSION=$(basename $BED_FILE)
CNV_FILENAME="${CNV_WITH_EXTENSION%.*}"

# Naming convention: name_nreads_purity
echo "[INFO] Merging files.."
samtools merge -f mixed_samples_final/${CNV_FILENAME}_$2_$PURITY_INT.bam bam_files/sub_$PID.bam bam_files/control_sub_$PID.bam

TARGET_CONTROL=$(samtools view -c mixed_samples_final/${CNV_FILENAME}_$2_$PURITY_INT.bam) # number of files in the subsampled cancer
TARGET_PROP=$( echo "$TARGET_CONTROL / $N_CONTROL" | bc -l )                                     

samtools view -s $TARGET_PROP -b $CONTROL > bam_files/control_sub_est_$PID.bam  
N_TARGET_SUB=$(samtools view -c bam_files/control_sub_est_$PID.bam)
echo "$N_TARGET_SUB" >> $METADATA_FILE

# here we normalize the samples using the control downsampled to have the same number of reads
Rscript ../process_normal_sample.R -n bam_files/control_sub_est_$PID.bam -c mixed_samples_final/${CNV_FILENAME}_$2_$PURITY_INT.bam

#echo "[INFO] Cleaning files.."
#rm bam_files/sub_$$.bam
#rm bam_files/control_sub_$$.bam
#rm bam_files/control_sub_est_$$.bam

echo "Finishing..."
