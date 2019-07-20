#!/bin/bash
#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=ps6      ### Job Name
#SBATCH --output=ps6.out         ### File in which to store job output
#SBATCH --error=ps6.err          ### File in which to store job error messages
#SBATCH --time=0-01:00:00	### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=7	### Number of cpus (cores) per task
#SBATCH --mail-user=abbiefayeolson@gmail.com
#SBATCH --mail-type=ALL

conda deactivate
conda deactivate
conda deactivate

conda activate bgmp_py3


fq_1=/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1
fq_2=/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2
fq_unmatched=/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched


velveth KMER31 31 \
-separate -fastq -shortPaired $fq_1 $fq_2 -short $fq_unmatched

velveth KMER41 41 \
-separate -fastq -shortPaired $fq_1 $fq_2 -short $fq_unmatched

velveth KMER49 49 \
-separate -fastq -shortPaired $fq_1 $fq_2 -short $fq_unmatched
