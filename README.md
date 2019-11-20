# TF-KOMPAS

-------
**This software is under development for publication, currently in the validation stage. It is not yet ready for general use.**
-------

Transcription Factor Kmers Oriented to Models Perform Aligned Searches is a program that calls binding sites using high-throughput kmer data from protein binding microarrays and SELEX-seq. It takes in the following input:

Files:
1. PWM model 
2. kmer data from PBM for SELEX-seq
3. A bed file to search for occurences

Parameters:
1. The minimum positions in the model required to fully describe the binding site. Called the core positions.
2. Binding signal cuttoff to use in the PBM or SELEX-seq data. 
3. A relative position in the model to center the calls on (optional, default is the middle of the PWM rounded down)

Output:
1. Centered sites in 0-base bed format with orientation (+/-)
2. Log file (optional)

Currently, the TF_KOMPAS_kmer_alignment produces aligned kmers that are then used in TF_KOMPAS_Site_Caller to call the sites.

See the ipython notebooks for specific instructions. Additional documentation available shortly. 



