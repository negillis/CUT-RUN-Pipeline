######## Annotated Slurm script for Processing CUT&RUN Data from .fastq to .narrowPeak #########
#### Author: Noelle Gillis, PhD
#### Date Modified: 2/22/23

#!/bin/bash
#SBATCH --job-name=cutnrun_pipeline
#SBATCH --partition= PartitionName
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-user= youremail.here@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=cutnrun_pipeline.log
#SBATCH --error=cutnrun_pipeline.log

cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo "$SLURM_ARRAY_TASK_ID"

# Load necessary modules - check what you have installed on your server vs. what is in your bash
# module load fastqc
# module load cutadapt
# module load samtools
# module load deeptools
# module load macs2
# module load multiqc

# Create directory for FastQC reports
mkdir fastqc_reports

# Create directory for trimmed fastq files
mkdir trimmed_fastq

# Run FastQC on all input FASTQ files and perform adapter trimming with TrimGalore
for f in *_R1_001.fastq.gz; do  #check the base name of your files and modify as needed
  base=$(basename ${f} _R1_001.fastq.gz)
  fastqc -o fastqc_reports ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz; # Generate FastQC report
  trim_galore --paired --illumina --fastqc --output_dir trimmed_fastq ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz; # Trim adapters using TrimGalore
done

# Create directory for BAM files
mkdir bam_files

# Map trimmed fastq files to reference genome using Bowtie2 and convert SAM to BAM format
for f in trimmed_fastq/*_R1_001_val_1.fq.gz; do
  base=$(basename $f _R1_001_val_1.fq.gz); # Extract basename of input file
  r2="$base""_R2_001_val_2.fq.gz"
  bowtie2 --very-sensitive-local --dovetail --no-unal --no-mixed --no-discordant -x /home/langeca/gilli431/software/GRCh38_noalt_as/GRCh38_noalt_as -1 $f -2 trimmed_fastq/$r2 | samtools view -Sb - > bam_files/$base.bam; # Map and convert to BAM
done

# Create directory for filtered BAM files
mkdir filtered_bam_files

# Filter BAM files for fragments that are 120 bp or less using samtools and awk
for f in bam_files/*.bam; do
  base=$(basename $f .bam); # Extract basename of input file
  samtools view -h $f | awk '{if ($9 <= 120) {print $0}}' | samtools view -b - > filtered_bam_files/$base.filtered.bam; # Filter BAM and save as new file
done

# Sort and index the filtered BAM files with samtools
for f in filtered_bam_files/*.bam; do
  base=$(basename $f .bam)
  samtools sort $f > filtered_bam_files/$base.sorted.bam
  samtools index filtered_bam_files/$base.sorted.bam
done

# Create directory for bigWig files
mkdir bigwig_files

# Generate bigWig files from filtered BAM files using bamCoverage
for f in filtered_bam_files/*sorted.bam; do
  base=$(basename $f .sorted.bam); # Extract basename of input file
  bamCoverage -b $f -o bigwig_files/$base.bw --binSize 10 --normalizeUsing RPKM --smoothLength 30 --extendReads 200; # Generate bigWig file
done

# Create directory for peak calls
mkdir peak_calls

# Call peaks using MACS2 on filtered BAM files
for f in filtered_bam_files/*.sorted.bam; do
  base=$(basename $f .sorted.bam); # Extract basename of input file
  macs2 callpeak --keep-dup=all -f BAMPE -t $f -g 2.9e9 -n $base --outdir peak_calls; # Call peaks using MACS2 and save output in new directory
done

# Run MultiQC to generate a quality control report of all the steps
multiqc . --filename multiqc_report.html --interactive

echo "Pipeline completed!"
