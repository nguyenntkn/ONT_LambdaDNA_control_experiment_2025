# Author: Tran Khanh Nguyen Nguyen
# Date: 21 March 2025

# Install FastQC for quality check, fast to trim reads, minimap2 for alignment, and
# samtools for converting into .bam alignment files.


# First need to open conda environment and activate it.
anaconda2023
conda create --prefix ./bioenv
conda activate bioenv

# Install FastQC for quality check of .fastq files optained after basecalling.
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
# Go to FastQC file and make it executable. Do this only once on the computer.
cd FastQC
chmod +x fastqc

# Install fastp to trim reads
conda create --name fastp_env   # Download..
conda activate fastp_env        # Activate fastq environment.
conda install -c bioconda fastp # Install fastq into conda environment.
fastp --version                 # Verify installation.
exit # NOTE: this will completely exit out of the program and the current terminal. 
     # To go back to conda environment and continue installing other tools, run anaconda2023.

# Install minimap2, a fast and versatile sequence alignment tool, useful for long read like Nanopore.
conda create --name minimap2_env   # Download.
conda activate minimap2_env        # Activate minimap2 environment.
conda install -c bioconda minimap2 # Install minimap2 into conda environment.
minimap2 --version                 # Verify installation.
exit # NOTE: this will completely exit out of the program and the current terminal. 
     # To go back to conda environment and continue installing other tools, run anaconda2023.
     
# Install samtools for converting into bam files.
conda create --name samtools_env   # Download.
conda activate samtools_env        # Activate minimap2 environment.
conda install -c bioconda samtools # Install minimap2 into conda environment.
samtools --version                 # Verify installation.
exit # NOTE: this will completely exit out of the program and the current terminal. 
     # To go back to conda environment and continue installing other tools, run anaconda2023.
