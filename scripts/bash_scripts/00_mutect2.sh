#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Memory request for 4G - 64G total
#$ -l h_vmem=16G

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
#$ -o /broad/thechenlab/ClaudiaC/droplet_DNA/results/logs/00_mutect2.o
#$ -e /broad/thechenlab/ClaudiaC/droplet_DNA/results/logs/00_mutect2.e


######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
reuse -q .gatk-4.2.3.0

##################
### Run script ###
##################

# tumor only
gatk Mutect2 \
  -R reference.fa \
  -I sample.bam \
  --germline-resource af-only-gnomad.vcf.gz \
  --panel-of-normals pon.vcf.gz \
  -O single_sample.vcf.gz
