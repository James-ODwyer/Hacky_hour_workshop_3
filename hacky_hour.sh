# Hacky Hour Workshop 3
# Hello! Welcome! The aim of the workshop, and this script, is to build a phylogenetic tree from nucleotide sequences on Petrichor, CSIRO's High Performance Computing Cluster (HPC). To do this we will be looking at this job submission script - usually job submission scripts will have a .sh (shell) suffix or .pbs. This script can be submitted to the job scheduler on the HPC, SLURM, to be run. Within this script we will provide SLURM with:
#  1) what computing parameters we would like to use for the job, 
#  2) what programs to load (allow us to use them in the job instance on the HPC),
#  3) some variables to be used in later code
#  4) run the tool MAFFT to align our sequences
#  5) run trimal to trim the resulting multiple sequence alignment (MSA)
#  6) run IQTree2 to build a phylogenetic tree from the trimmed MSA

#SBATCH --job-name=hacky3
#SBATCH --time=05:00:00
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --account=OD-220926

# Module load
module load mafft/7.490
module load trimal/1.4.1r22
module load iqtree/2.2.0.5

# BASH variable used when running tools. It is best to specify file and directory loactions in variables
# like this so that you can easily switch up the target and parameters in your commands.
threads=12 # You can use the $SLURM_THREADS_PER_CORE if you want instead 
name='human_adeno_5'
refs='sequences/human_adeno_5_genomes.fa'
target='sequences/mystery_virus.fa'
outgroup='M73260.1_Mastadenovirus_h5_gene'

# Put reference sequences and our mystery virus in the same fasta file, ready for alignment.
cat $refs $target > $name.fa

# Perform multiple sequence alignment using MAFFT:
# --nuc: Indicates that the input data is nucleotide sequences.
# --auto: Automatically selects an appropriate algorithm.
# --reorder: Reorders input sequences to optimize alignment.
# --thread $threads: Specifies the number of threads for parallel processing.
mafft --nuc --auto --reorder --thread $threads $name.fa > alignment/$name.aln

# Trim the alignment using trimAl:
# -in $name.aln: Input alignment file.
# -out $name.trimaln: Output trimmed alignment file.
# -gt 0.8: Remove columns with more than 80% gaps.
# -st 0.001: Remove sequences with less than 0.1% similarity to others.
trimal -in alignment/$name.aln -out alignment/$name.trimaln -gt 0.8 -st 0.001

# Build a phylogenetic tree using IQ-TREE 2:
# -s $name.trimaln: Input alignment file.
# --prefix $name.tree: Prefix for output tree file.
# -m TEST: Model selection (IQ-TREE will choose the best model).
# --alrt 1000: Perform ultrafast bootstrap approximation with 1000 replicates.
# -B 1000: Perform standard bootstrap analysis with 1000 replicates.
# -T $threads: Number of threads for parallel computation.
# -o $outgroup: Designate the specified taxon as the outgroup for tree rooting.
iqtree2 -s alignment/$name.trimaln  --prefix tree/$name.tree -m TEST --alrt 1000 -B 1000 -T $threads -o $outgroup

## Last, we can run this script on the HPC.
# Log in by:
#    iden@petrichor-login.hpc.csiro.au
# Navigate to where you would like to run your script.
# Run this script suing:
#    sbatch hacky_hour.sh
