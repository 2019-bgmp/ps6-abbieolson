#!/bin/bash
#SBATCH --account=bgmp                  ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp                ### Partition (like a queue in PBS)
#SBATCH --job-name=ps6                  ### Job Name
#SBATCH --output=slurm-%j-%x.out        ### File in which to store job output
#SBATCH --time=0-00:15:00	              ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                       ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1             ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1	              ### Number of cpus (cores) per task
#SBATCH --mail-user=abbiefayeolson@gmail.com
#SBATCH --mail-type=ALL

conda deactivate
conda deactivate
conda deactivate

conda activate bgmp_py3

velvetg /home/afo/bgmp/Bi621/PS6/KMER31/ -cov_cutoff 20
mv /home/afo/bgmp/Bi621/PS6/KMER31/contigs.fa /home/afo/bgmp/Bi621/PS6/KMER31/contigs_31_cov_20.fa

velvetg /home/afo/bgmp/Bi621/PS6/KMER31/ -cov_cutoff 60
mv /home/afo/bgmp/Bi621/PS6/KMER31/contigs.fa /home/afo/bgmp/Bi621/PS6/KMER31/contigs_31_cov_60.fa

velvetg /home/afo/bgmp/Bi621/PS6/KMER31/ -cov_cutoff auto
mv /home/afo/bgmp/Bi621/PS6/KMER31/contigs.fa /home/afo/bgmp/Bi621/PS6/KMER31/contigs_31_cov_auto.fa
