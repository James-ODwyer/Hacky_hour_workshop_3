#!/bin/bash
# Hacky Hour Workshop 3
###########################
# Introduction
###########################
# Hello! Welcome! The aim of the workshop, and this script, is to build a phylogenetic tree from 
# nucleotide sequences on Petrichor, CSIRO's High Performance Computing Cluster (HPC). To do this we will
# be looking at this job submission script - usually job submission scripts will have a .sh (shell) 
# suffix or .pbs. This script can be submitted to the job scheduler on the HPC, SLURM, to be run. Within 
# this script we will provide SLURM with:
#  1) what computing parameters we would like to use for the job, 
#  2) what programs to load (allow us to use them in the job instance on the HPC),
#  3) some variables to be used in later code
#  4) run the tool MAFFT to align our sequences
#  5) run trimal to trim the resulting multiple sequence alignment (MSA)
#  6) run IQTree2 to build a phylogenetic tree from the trimmed MSA
# Then outside of this script we will take the trees files we create and make a figure using ITol.

###########################
# 1) Running on HPC
###########################
# First things first, if you have not already, log in to your hpc account with:
#    > ssh iden@petrichor-login.hpc.csiro.au
# Navigate to where you would like to run your script:
#    > cd ~/script/location
# Get access to this workshop repository through git with:
#    > git clone https://github.com/James-ODwyer/Hacky_hour_workshop_3
#    > cd Hacky_hour_workshop_3
# Then you can run this script with:
#    > sbatch hacky_hour_workshop_3.sh
#
# The following commented lines the arguments that SBatch takes when creating a job instance. I believe 
# that these have been covered in the previous workshop, so I won't go into detail on their meaning here.
# What is important is that we do not need much memory or cpu to run all of the following tools with only
# ~12 sequences. 1G RAM and 5 CPU cores is enough.

#SBATCH --job-name=hacky3
#SBATCH --time=05:00:00
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --account=OD-220926

###########################
# 2) Loading in software
###########################
# On a computing cluster you make the programs you want to use available as needed, you load in the 
# programs you need to use in a job at the start. Again, I believe this has been described in the 
# previous workshop. If you ever need to see if Petrichor has a program you need, run 'module avail'.
#
# For this job, we will need the aligning program MAFFT, the MSA trimming tool trimAl, and the 
# maximum-likelihood inferred tree building program IQTree2 (which is excellent). Note that I specify
# what versions of these programs I use in the commands, because knowing exactly what version you used
# when running something like this is a good practice. It will save you hassle with fixing bugs or when
# you are writing about your method down the track.

module load mafft/7.490
module load trimal/1.4.1r22
module load iqtree/2.2.0.5

###########################
# 3) Setting our variables
###########################
# If you have some text that is used in repeated commands, put it as a variable. If you have some parts 
# of your code that can be changed to adapt your script to a new task, put it as a variable. Generally, 
# if you have any constants that might be clearer with a name (like 'pval=0.005' rather than just giving
# '0.005', put it as a variable.
#
# The following is the bits that we would like adjustable for this script: the amount of threads we are 
# going to tell the programs we are using, the name of input files, and the sequence we think would be a 
# good outgroup.

threads=12 # You can use the $SLURM_THREADS_PER_CORE if you want instead 
name='human_adeno_5'
refs='sequences/human_adeno_5_genomes.fa'
target='sequences/mystery_virus.fa'
outgroup='M73260.1_Mastadenovirus_h5_gene'


###############################
# Description of task and data
###############################
# So we want to build a phylogenetic tree, but why? We are building a tree to identify what our mystery 
# sequence is (mystery_virus.fa). We know that it is a DNA virus genome and that it is likely in the 
# mastadenovirus family, probably a human adenovirus. So we have collected some prospective relatives of
# our mystery virus: some Human adenoviruss in human_adeno_all_genomes.fa.
#
# In the following we will take our target virus sequence and the reference sequences and align them
# all together in a MSA. We will trim this alignment so that only informative alignment sites remain,
# and then we will build a tree from this alignment using an algorithm that maximises the likelihood of
# seeing our data (our sequences and how they align) given a model (a tree structure). After that we
# can view the tree with ITol and see if we can come to a conclusion of what our mystery virus is, based
# on how it is related to the other tree members.
# Also, we will be designating an outgroup sequence. As trees don't inherently have directionality - you
# can find how much evolution occured between three members, but you can't say which have a most recent 
# common ancestor with out more clues - we can specify a member that we believe to be the most distant 
# relative, which "roots" the tree, giving directionality. If this is confusing I recommend giving it a
# google, there are good resources that can explain this (with pictures) better than I can.
#
# So, before we get into it we have to solve a problem. We will need to provide a fasta file as input to
# MAFFT as the first step, but we have two fasta files, a reference set of sequences and our mystery 
# virus. We will use the 'cat' command (for conCATenate) to print out the content of both our fasta 
# files, but instead of sending it the terminal we will send the content to a new file using the '>' 
# pipe. Try runnning the following command without the '>' and subsequent code to get a feel for it if
# you are confused. Use 'ls' and 'less' to see what we have created.

cat $refs $target > sequences/$name.fa

#################################
# 4) Multiple sequence alignment
#################################
# Perform multiple sequence alignment using MAFFT:
# --nuc: Indicates that the input data is nucleotide sequences.
# --auto: Automatically selects an appropriate algorithm for alignment (gap and mismatch thresholds,
# appropriate substituion matrix).
# --reorder: Reorders input sequences so that they ordered by similarity.
# --thread $threads: Specifies the number of threads for parallel processing.

mafft --nuc --auto --reorder --thread $threads sequences/$name.fa > alignment/$name.aln

###########################
# 5) Trimming the MSA
###########################
# Trim the alignment using trimAl:
# -in $name.aln: Input alignment file.
# -out $name.trimaln: Output trimmed alignment file.
# -gt 0.8: Remove columns with more than 80% gaps.
# -st 0.001: Remove sequences with less than 0.1% similarity to others.
#
# Alignment trimming is a very important part of this process, the idea is to cut out bits of the MSA 
# that would be a waste of compute time, or misleading for analysis. It is not always necessary, sometime
# you align sequences and they have minimal gaps and are at that sweet spot of not identical, but not 
# completely different, that is most informative. It is, however, VITAL that you look at your alignments
# when making phylogenetic trees - make sure they look good! 
# 
# I recommend taking the $name.aln and $name.trimaln files off to an alignment viewer (likely Geneious)
# and having a look at the MSA. Grab someone from the workshop to discuss it with you!

trimal -in alignment/$name.aln -out alignment/$name.trimaln -gt 0.8 -st 0.001

#################################
# 6) Phylogenetic tree building
#################################
# Build a phylogenetic tree using IQTree 2:
# -s $name.trimaln: Input alignment file.
# --prefix $name.tree: Prefix for output tree file.
# -m TEST: Evolution model selection (IQ-TREE will choose the best model).
# --alrt 1000: Run SH-aLRT test - it similar to bootstraps but more efficient, basically a likelihood of
# internal nodes being where they are. From the paper: "This approach extends the recently proposed 
# approximate likelihood-ratio test and relies on a nonparametric, Shimodaira–Hasegawa–like procedure." 
# Do not ask me questions about this.
# -B 1000: Perform standard bootstrap analysis with 1000 replicates - that is make 1000 tree strucutures 
# that fit the data well and count how many time branches appear. The more often they appear, the more 
# 'stable' that grouping of members is.
# -T $threads: Number of threads for parallel computation.
# -o $outgroup: Designate the specified taxon as the outgroup for tree rooting.
#
# So, I don't exactly want to explain how phylogenetic tree building works because I will just do a worse
# job than the googlable resources out there, but I'll give it a quick go so that you understand some 
# terms: IQtree take your MSA and estimates (using an evolutionary model) the amount of evolution 
# undergone between all sequences. It then generates trees, fits the evolutionary distances it calculates
# to it, and calculates# how likely it is that the sequences have evolved this much, with respect to each
# other, given the tree structure.
# I did my best, it's a very complicated topic. Best read about general tree building and evolutionary 
# models, then about maximum-likelihood algorithmns specifically, then maybe the IQTree paper.
#
# Note that the TEST argument makes things a lot slower, so if you know ahead of time what evolutionary
# model that this tree will use, best to specify it.

iqtree2 -s alignment/$name.trimaln  --prefix tree/$name -m TEST --alrt 1000 -B 1000 -T $threads -o $outgroup

#################################
# Making figures from tree files
#################################
# We get quite a few files from IQTree, the important one being the .treefile. This contains your tree
# structure in newick format (pretty sure), and we can use this file to generate our tree figures. It 
# describes many nodes, the distances between the nodes (in substitutions per site), and any labels
# attached to them - the terminal nodes will have sequence names attached and the internal nodes will
# have boostrap and aLRT scores attached.
#
# Download the tree file. If you are using a mac terminal or Windows Subsytem for Linux, you can use:
#     > scp iden@petrichor.hpc.csiro.au:~/some/path/to/Hacky_hour_workshop_3/tree/human_adeno_all.treefile ./
# Once you have the file on you PC, navigate to the ITol webserver at:
#    https://itol.embl.de/
# Drag and drop your tree into the upload page, and you should be able to view your tree. Have a go with
# an unrooted circular tree view, then try using a specific root. What is the mystery virus?
