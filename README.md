# SOMatic-Network-Analysis

This package is for performing network analysis on Linked Metaclusters from SOMatic

1.) Build reference genome background file from fasta reference genome:
fasta-get-markov hg38.fa > hg38_bg

2.) Run MotifSearch program (choose the one that fits your server architecture - if one doesn't exist, modify the script to work on your server)
./MotifSearchLM-hpc.sh -LinkFolder <Linked Metacluster Folder> -Analysis <Name of this analysis> -Metacluster1 <Number of DNA metaclusters> -Metacluster2 <Number of RNA metaclusters> -MotifDatabase <MEME file> -ReferenceGenome <Ref genome fasta> -Thresh <q value threshold> -BGFile <location of the ref genome background>

3.) Run Motif enrichment
./MotifEnrichment -Metacluster1 <Number of DNA metaclusters> -Metacluster2 <Number of RNA metaclusters> -LinkFolder <Location of the fusion breakup folder> -Analysis <Name of Motif Analysis> -pval <Desired zscore pvalue>
