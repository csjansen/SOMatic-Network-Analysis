
This package is for performing network analysis on Linked Metaclusters from SOMatic

Download the repo and make it:
`git clone http://github.com/csjansen/SOMatic-Network-Analysis`
`cd SOMatic-Network-Analysis`
`make`

1.) Build reference genome background file from fasta reference genome:
fasta-get-markov hg38.fa > hg38_bg

2.) Run MotifSearch program (choose the one that fits your server architecture - if one doesn't exist, modify the script to work on your server)
./MotifSearchLM-hpc.sh -LinkFolder <Linked Metacluster Folder> -Analysis <Name of this analysis> -Metacluster1 <Number of DNA metaclusters> -Metacluster2 <Number of RNA metaclusters> -MotifDatabase <MEME file> -ReferenceGenome <Ref genome fasta> -Thresh <q value threshold> -BGFile <location of the ref genome background>

3.) Run Motif enrichment
./MotifEnrichment.sh -Metacluster1 <Number of DNA metaclusters> -Metacluster2 <Number of RNA metaclusters> -LinkFolder <Location of Linked Metacluster folder> -Analysis <Name of Motif Analysis> -pval <Desired zscore pvalue>

4.) Run Network generation
./MakeNetwork.sh -Analysis <Name of Motif Analysis> -LinkFolder <Location of Linked Metacluster folder> -Metacluster1 <Number of DNA metaclusters> -Metacluster2 <Number of RNA metaclusters> -Output <Output File Location>

Add a MotifFile option if your motif database includes a TSV to convert a motif ID to a gene name (ala HOCOMOCO)

