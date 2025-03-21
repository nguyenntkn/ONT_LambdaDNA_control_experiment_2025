

# Before proceeding, create a appropriate working directory:
mk /home/user/Work_stuff
cd /home/user/Work_stuff

# In the ONT run, we have raw sequencing files.
# If enabled basecalling, the resulting files will be written into a folder called "fastq_pass"
# Basecalled files are stored as .fastq.gz (compressed?)
# fastq_pass folder stores reads with high quality score, as opposed to fastq_fail.
# I like to create a copy of that folder (fastq_pass_backup) and use this for analysis.

# Data processing:
# To use cocalc or other cloud-based services, we need to compress the folder containing 
# sequencing data and upload to cloud environment.
# Unzip that folder:
unzip fastq_pass_backup.zip

# STEP 1: Running FastQC for quality control
# Need to open conda environment to access fastqc tool.
anaconda2023

# Activate fastqc environment
conda activate fastqc_env

# Create a directory to store QC reports.
mkdir fastqc_report

# Run FastQC on all files with .fastq end. (-o: specifies output directory)
/home/user/FastQC/fastqc fastq_pass/*.fastq -o fastqc_reports/

# Exit fastqc environment.
# NOTE: this will also exit conda and the current terminal.
exit

# Need to manually go back to anaconda2023 environment.
anaconda2023

# STEP 2: Use fastp to trim or filter out problematic sequences.
# Activate fastp environment
conda activate fastp_env

# Trim or filter out problematic sequences.
# Use a loop to loop through all fastq files. 
mkdir trimmed_fastq
for file in fastq_pass/*.fastq; do
  fastp -i "fastq_pass/$file" -o "trimmed_fastq/trimmed_${file}"
done

# Exit fastqc environment.
# NOTE: this will also exit conda and the current terminal.
exit

# Need to manually go back to anaconda2023 environment.
anaconda2023

# Optional: fastqc again to confirm the quality improvements.============== HERE!!
mkdir trimmed_qc_fastq
/home/user/FastQC/fastqc trimmed_fastq/*.fastq -o trimmed_qc_fastq/

# Get reference genome
# This one doesn't really work. Need some fixing
# ref_genome_lambda=https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1?report=fasta
# wget $ref_genome_lambda -O lambda_genome.fasta

# Index genome
# minimap2 -d lambda_genome.mmi /home/user/Work_stuff/lambda_genome.fasta 


# mkdir /home/user/Work_stuff/Sam_files
# Use minimap2 for alignment




for file in /home/user/Work_stuff/trimmed_qc_fastq/*.fastq; do
    filename=$(basename "$file")  # Extract the filename from the full path
    minimap2 -a /home/user/lambda_genome.mmi "$file" > "/home/user/Work_stuff/Sam_files/${filename%.fastq}_aligned.sam"
done





for file in /home/user/Work_stuff/Sam_files/*.sam; do
    # Extract the base filename without extension
    filename=$(basename "$file" .sam)
    
    # Convert SAM to BAM
    samtools view -Sb "$file" > "/home/user/Work_stuff/Bam_files/${filename}.bam"
    
    # Sort the BAM file
    samtools sort "/home/user/Work_stuff/Bam_files/${filename}.bam" -o "/home/user/Work_stuff/Bam_files/${filename}_sorted.bam"
    
    # Optional: Remove unsorted BAM file to save space
    rm "/home/user/Work_stuff/Bam_files/${filename}.bam"
done


for bam_file in /home/user/Work_stuff/Bam_files/*.bam; do
    samtools index "$bam_file" "/home/user/Work_stuff/Bam_files/$(basename "$bam_file").bai"
done













