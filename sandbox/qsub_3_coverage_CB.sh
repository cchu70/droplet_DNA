#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Memory request for 4G - 64G total
#$ -l h_vmem=24G

# Cores
#$ -pe smp 16
#$ -binding linear:16

# I like single output files
#$ -j y

# Runtime request
#$ -l h_rt=24:00:00

# Email
#$ -M cchu@broadinstitute.org
#$ -m bea  # Sends emails at the beginning (b), end (e), and if the job is aborted (a)

#$ -cwd
#$ -V
#$ -o /broad/thechenlab/ClaudiaC/droplet_DNA/sandbox/results/03/logs/cytoband_coverage.o
#$ -e /broad/thechenlab/ClaudiaC/droplet_DNA/sandbox/results/03/logs/cytoband_coverage.e


######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
reuse -q Anaconda3

##################
### Run script ###
##################

conda activate droplet_dna_py39

# run from sandbox
python /broad/thechenlab/ClaudiaC/droplet_DNA/sandbox/src/binned_coverage.py --bam /broad/thechenlab/Benno/experiments/xBO153/bams/xBO153_G.markdup.bam \
    --bin_size 1000000 \
    --num_workers 16 \
    --output_fn /broad/thechenlab/ClaudiaC/droplet_DNA/sandbox/results/03/binned_coverage.by_CB.txt
