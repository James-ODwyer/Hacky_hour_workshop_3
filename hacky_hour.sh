#!/bin/bash
# Hacky Hour Workshop 4 - Script to align sequences, trim the alignment, and build a tree.

#SBATCH --job-name=hackhour4
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

indir='References_mastadenovirusC5_strain5'
input='Human_mastadenovirus_C5_all.fasta' 

outgroup='M73260.1'


# Perform multiple sequence alignment using MAFFT:
# --nuc: Indicates that the input data is nucleotide sequences.
# --auto: Automatically selects an appropriate algorithm.
# --reorder: Reorders input sequences to optimize alignment.
# --thread $threads: Specifies the number of threads for parallel processing.
mafft --nuc --auto --reorder --thread $threads $indir/$input > alignment/${input}.aln

# Trim the alignment using trimAl:
# -in ${input}.aln: Input alignment file.
# -out ${input}.aln.trim: Output trimmed alignment file.
# -gt 0.8: Remove columns with more than 80% gaps.
# -st 0.001: Remove sequences with less than 0.1% similarity to others.
trimal -in alignment/${input}.aln -out alignment/${input}.aln.trim -gt 0.8 -st 0.001

# Build a phylogenetic tree using IQ-TREE 2:
# -s ${input}.aln.trim: Input alignment file.
# --prefix ${input}.aln.trim.tree: Prefix for output tree file.
# -m TEST: Model selection (IQ-TREE will choose the best model).
# --alrt 1000: Perform ultrafast bootstrap approximation with 1000 replicates.
# -B 1000: Perform standard bootstrap analysis with 1000 replicates.
# -T $threads: Number of threads for parallel computation.
# -o $outgroup: Designate the specified taxon as the outgroup for tree rooting.
iqtree2 -s alignment/${input}.aln.trim  --prefix tree/${input}.aln.trim.tree -m TEST --alrt 1000 -B 1000 -T $threads -o $outgroup
