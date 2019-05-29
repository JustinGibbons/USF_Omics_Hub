#!/bin/bash
#SBATCH --ntasks=3
#SBATCH --workdir=/work/j/jgibbons1/Sample_RNA-Seq_on05.28.2019at13.36_6/
#SBATCH --mail-type=ALL
#SBATCH --time=1:55:00
#SBATCH --mem=10000
#SBATCH --qos=rra
#SBATCH --nodes=1
#SBATCH --partition=rra
#SBATCH --mail-user=JGibbons1@mail.usf.edu
#SBATCH --job-name=hisat_pipeline
#SBATCH --array=0-13%14
#SBATCH --output=rna-seq_pipeline.out
module purge
module load apps/hisat2/2.1.0
module load apps/samtools/1.3.1
module load apps/cufflinks/2.2.1
module load apps/subread/1.6.3

declare -a unpaired_array=("/work/j/jgibbons1/Richards/Samples/P2.fastq.gz" "/work/j/jgibbons1/Richards/Samples/P3.fastq.gz" "/work/j/jgibbons1/Richards/Samples/P4.fastq.gz" "/work/j/jgibbons1/Richards/Samples/AC2.fastq.gz" "/work/j/jgibbons1/Richards/Samples/AC4.fastq.gz" "/work/j/jgibbons1/Richards/Samples/AC1.fastq.gz" "/work/j/jgibbons1/Richards/Samples/AC3.fastq.gz" "/work/j/jgibbons1/Richards/Samples/M4.fastq.gz" "/work/j/jgibbons1/Richards/Samples/C_3.fastq.gz" "/work/j/jgibbons1/Richards/Samples/C_5.fastq.gz" "/work/j/jgibbons1/Richards/Samples/C_4.fastq.gz" "/work/j/jgibbons1/Richards/Samples/M3.fastq.gz" "/work/j/jgibbons1/Richards/Samples/M82-2.fastq.gz" "/work/j/jgibbons1/Richards/Samples/M1.fastq.gz")


declare -a hisat_outfile_array=("Alignments/P2.sam" "Alignments/P3.sam" "Alignments/P4.sam" "Alignments/AC2.sam" "Alignments/AC4.sam" "Alignments/AC1.sam" "Alignments/AC3.sam" "Alignments/M4.sam" "Alignments/C_3.sam" "Alignments/C_5.sam" "Alignments/C_4.sam" "Alignments/M3.sam" "Alignments/M82-2.sam" "Alignments/M1.sam")


declare -a sorted_bam_files_array=("Alignments/P2_sorted.bam" "Alignments/P3_sorted.bam" "Alignments/P4_sorted.bam" "Alignments/AC2_sorted.bam" "Alignments/AC4_sorted.bam" "Alignments/AC1_sorted.bam" "Alignments/AC3_sorted.bam" "Alignments/M4_sorted.bam" "Alignments/C_3_sorted.bam" "Alignments/C_5_sorted.bam" "Alignments/C_4_sorted.bam" "Alignments/M3_sorted.bam" "Alignments/M82-2_sorted.bam" "Alignments/M1_sorted.bam")


declare -a cufflinks_output_array=("P2_cufflinks_output" "P3_cufflinks_output" "P4_cufflinks_output" "AC2_cufflinks_output" "AC4_cufflinks_output" "AC1_cufflinks_output" "AC3_cufflinks_output" "M4_cufflinks_output" "C_3_cufflinks_output" "C_5_cufflinks_output" "C_4_cufflinks_output" "M3_cufflinks_output" "M82-2_cufflinks_output" "M1_cufflinks_output")


hisat2 --rna-strandness F -x /work/j/jgibbons1/Richards/Index/Solanum_lycopersicum.SL3.0.dna.toplevel.fa -U ${unpaired_array[$SLURM_ARRAY_TASK_ID]} -S ${hisat_outfile_array[$SLURM_ARRAY_TASK_ID]}


samtools sort -@ 3 -o ${sorted_bam_files_array[$SLURM_ARRAY_TASK_ID]} ${hisat_outfile_array[$SLURM_ARRAY_TASK_ID]}

echo ${sorted_bam_files_array[$SLURM_ARRAY_TASK_ID]}
cufflinks -q -p 3 -o ${cufflinks_output_array[$SLURM_ARRAY_TASK_ID]} ${sorted_bam_files_array[$SLURM_ARRAY_TASK_ID]}


featureCounts -T 3 -t gene -g ID -a /work/j/jgibbons1/Richards/Index/Solanum_lycopersicum.SL3.0.43.gff3 -o feature_counts_output.txt -p Alignments/AC1_sorted.bam Alignments/AC2_sorted.bam Alignments/AC3_sorted.bam Alignments/AC4_sorted.bam Alignments/C_3_sorted.bam Alignments/C_4_sorted.bam Alignments/C_5_sorted.bam Alignments/M1_sorted.bam Alignments/M3_sorted.bam Alignments/M4_sorted.bam Alignments/M82-2_sorted.bam Alignments/P2_sorted.bam Alignments/P3_sorted.bam Alignments/P4_sorted.bam


cuffnorm -q -p 4 -o Cuffnorm_Output --library-norm-method classic-fpkm /work/j/jgibbons1/Richards/Index/Solanum_lycopersicum.SL3.0.43.gff3 -L AC,C_,M,P Alignments/AC1_sorted.bam,Alignments/AC2_sorted.bam,Alignments/AC3_sorted.bam,Alignments/AC4_sorted.bam Alignments/C_3_sorted.bam,Alignments/C_4_sorted.bam,Alignments/C_5_sorted.bam Alignments/M1_sorted.bam,Alignments/M3_sorted.bam,Alignments/M4_sorted.bam,Alignments/M82-2_sorted.bam Alignments/P2_sorted.bam,Alignments/P3_sorted.bam,Alignments/P4_sorted.bam
