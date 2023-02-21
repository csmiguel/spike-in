#!/bin/bash
#trim adaptors
# Cutadapt cutadapt 2.10 with Python 3.6.11 was used to trim PCR primers and
# write to ouput with only paired sequences for which both primers were found.
#Run loop to (1) trim forward primer anchored to 5' end in R1 reads and
# reverse primer anchored to 5' end in R2 reads
# cutadapt defaults
# - empty reads are kept and will appear in the output. else use -m
# - maximum error rate (default: 0.1) -e 0.1
# - minimum overlap (default: 3)  -O 3
#
### 16S ###
# primers 515F-Y (5’ GTGYCAGCMGCCGCGGTAA 3’) (Parada, Needham, and Fuhrman 2016)
# and 806R (5’ GGACTACNVGGGTWTCTAAT 3’) (Apprill et al. 2015).
find data/raw/16S*R1.fastq.gz | xargs basename | sed 's/_R1.fastq.gz//' | while read infasta
do
  cutadapt -g ^GTGYCAGCMGCCGCGGTAA -G ^GGACTACNVGGGTWTCTAAT --trimmed-only \
  -o data/intermediate/"$infasta"_R1.cutadapt.fastq.gz \
  -p data/intermediate/"$infasta"_R2.cutadapt.fastq.gz \
  -m 50 \
  data/raw/"$infasta"_R1.fastq.gz data/raw/"$infasta"_R2.fastq.gz
done > output/16S_cutadatp_filtering.txt

### ITS ###
# For fungal library preparation, a fragment of the ITS2 genomic region of around
# 350 bp was amplified using the primers
# ITS3 (5’ GCATCGATGAAGAACGCAGC 3’) and ITS4R (5’ TCCTCCGCTTATTGATATGC 3’)
# (White et al. 1990). However, ITS has introns and its size is very variable.
# Thus the target region can vary in such a way that some reads can even extend over
# the reverser complementary primer in the 3' end.
# I detected a large proportion of the reads carrying the RC primer. Therefore, and
# as recommended by Callahan in the ITS tutorial, I trimmed the primer from 3' ends. I Run
# cutadapt in 2 steps because I am interested in the "--trimmed-only" option to be applied
# only to the anchored 5' primers.
find data/raw/ITS*R1.fastq.gz | xargs basename | sed 's/_R1.fastq.gz//' | while read infasta
do
  # remove 5' primer
  cutadapt -g ^GCATCGATGAAGAACGCAGC -G ^TCCTCCGCTTATTGATATGC --trimmed-only \
  -o data/intermediate/"$infasta"_R1.temp.fastq.gz \
  -p data/intermediate/"$infasta"_R2.temp.fastq.gz \
  data/raw/"$infasta"_R1.fastq.gz data/raw/"$infasta"_R2.fastq.gz
  # remove 3' primer
  cutadapt -a GCATATCAATAAGCGGAGGA \
    -A GCTGCGTTCTTCATCGATGC \
    -o data/intermediate/"$infasta"_R1.cutadapt.fastq.gz \
    -p data/intermediate/"$infasta"_R2.cutadapt.fastq.gz \
    -m 50 \
    data/intermediate/"$infasta"_R1.temp.fastq.gz data/intermediate/"$infasta"_R2.temp.fastq.gz
  # remove temp files
  rm data/intermediate/"$infasta"_R1.temp.fastq.gz
  rm data/intermediate/"$infasta"_R2.temp.fastq.gz
done > output/its_cutadatp_filtering.txt
