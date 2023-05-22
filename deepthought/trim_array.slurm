#!/bin/bash
#SBATCH --job-name=Trimming
#SBATCH --time=5-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -o trimming_slurm/trim_adap-%A_%a.out
#SBATCH -e trimming_slurm/trim-adap-%A_%a.err

# submit this with sbatch --array=1-127:1 ./trim_array.slurm

R1=$(head -n $SLURM_ARRAY_TASK_ID R1_reads.txt | tail -n 1)
R2=${R1/R1/R2}
R1trimmed=${R1/.fastq.gz/.trimmed.fastq.gz}
R2trimmed=${R2/.fastq.gz/.trimmed.fastq.gz}
R1hits=${R1/.fastq.gz/.hits.tsv}
R2hits=${R2/.fastq.gz/.hits.tsv}
R2adjust=${R2/.fastq.gz/.adjustments}

~/GitHubs/primer-trimming/bin/search-paired-snps -l fastq/$R1 -r fastq/$R2 -p fastq_trimmed/$R1trimmed -q fastq_trimmed/$R2trimmed -j fastq_adapter_matches/$R1hits -k fastq_adapter_matches/$R2hits -a fastq_adapter_matches/$R2adjust
