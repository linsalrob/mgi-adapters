# Using trim-mgi-adapters on Deepthought.

_Note:_ Deepthought is Flinders' HPC, and is pretty much a standard HPC running `slurm`. These instructions may work for you, but they may not.

This assumes that you have a directory called `fastq` with your R1 and R2 reads, and that you want the trimmed fastq files in a directory called `fastq\_trimmed`. We will also put the reports of which sequences find which files into a directory called fastq\_adapter\_matches. Finally, we use a directory called trimming\_slurm for the slurm output files.

# Step 1. Clone the repo and build the code

In your account, clone the repository and build all the code. No modules are needed for this:

```
git clone https://github.com/linsalrob/mgi-adapters.git
cd mgi-adapters
make all
```

# Step 2. Create a file with all the R1 reads, and just their names

```
cd path/where/sequences/are
find fastq -type f -name \*R1* -printf "%f\n" > R1_reads.txt
```

# Step 3. Create a directory for the slurm output files

```
mkdir trimming_slurm
```

**Note:** If you don't do this, the `sbatch` command below will run normally and you will think everything is fine, although it finishes very quickly and has not done anything!

# Step 4. Find out how many R1 reads we have

```
READS=$(wc -l R1_reads.txt | awk '{print $1}');
```

# Step 5. Submit the array job to process all those reads. 

Note that the slurm script makes the directories for fastq\_trimmed and fastq\_adapter\_matches. It also handles both the R1 and R2 files.

```
sbatch --array=1-$READS:1 ~/GitHubs/mgi-adapters/deepthought/trim_array.slurm
```

# Wait for the results!

Here is the whole command in one line, so you can just copy and paste it!

```
mkdir trimming_slurm; find fastq -type f -name \*R1* -printf "%f\n" > R1_reads.txt; READS=$(wc -l R1_reads.txt | awk '{print $1}'); sbatch --array=1-$READS:1 ~/GitHubs/mgi-adapters/deepthought/trim_array.slurm
```



