# TF K.O.M.P.A.S.
Kmers Oriented to Models Perform Aligned Searches

-------
**This software is under development for publication, currently in the validation stage. It is not yet ready for general use.**
-------

Transcription Factor Kmers Oriented to Models Perform Aligned Searches is a suite of tools that calls binding sites using high-throughput kmer data from protein binding microarrays and SELEX-seq. It uses a novel algorithm the builds on previous kmer overlap methods. Currently, the following programs are includes:

1. KOMPAS: Kmer Alignment
2. KOMPAS: Site Caller
3. KOMPAS: Gapped Dimeric Calls

Kmer Alignment produces an aligned kmer file that can be used with the Site Caller or Gapped Dimeric Calls scripts. For most transcription factors, even dimeric palindromes, Site Caller is appropriate for use. When the transcription factor is a gapped dimeric palindrome with a full model size that is larger than k to the point where aligned kmers of high enrichment mostly (or only) describe half sites, the Gapped Dimeric Calls script can be used. 

Pipeline Output:
1. Centered sites in 0-base bed format with orientation (+/-)
2. log file (optional)

Currently, the TF_KOMPAS_kmer_alignment produces aligned kmers that are then used in TF_KOMPAS_Site_Caller to call the sites.

See the ipython notebooks for specific instructions. Additional documentation available shortly. 



