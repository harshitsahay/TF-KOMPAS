'''
### TF-KOMPAS: Site caller

Updated 1/14/20

Author: Zachery Mielko

Contributors and Collaborators:
* Tiffany Ho: optimization review, feedback, and validation
* Farica Zhuang: review, feedback, and validation

-----

The script requires the following dependencies:

* pandas
* numpy
* biopython
* pybedtools

The script takes the following as input:
- Sequence information in the form of:
    - A bed file of genomic coordinates with a genome file in FASTA format
    - A FASTA file of sequences
- An aligned kmer file

With the following parameters
- Core definition
- Escore threshold
- If the transcription factor is palindromic or not
- A position to center the calls on (kPosition)

The script gives the following as **output**:
- Bed file of centered sites (Centered_KOMPAS.bed)
    - The Bed file is 0 based and of the format start inclusive, end exclusive. So half-open.
    - Columns: Chromosome, Start, End, Orientation, Threshold Score
        - If given a FASTA file, KOMPAS will use the ID as the first column instead
- Log file (optional)
- Scoring Table (optional)
    - Useful for troubleshooting or validation
    - Gives the score for all positions in all sequences as a tsv file. ScoreTable.tsv
    - Gives the calls for all sequences in a tsv format. CallTable.tsv

By default, the center position will be the midpoint of the core (rounded down), but it can be specifically chosen as an optional parameter.

**Core**

The core is the range of positions in the model that must be described in order to call a binding site. This is a prior given by the user based on biological context. The longer the core, the more strict the site calling is. It does not affect the threshold score given to a site. 

**Threshold Score**

This is the threshold for when the site would be called. It is a score used for the calling of a site, but not for the scoring of how well a site is bound. 

If k > core length (if multiple kmers exist that can fully describe the core), then the second score from **maximum** for those kmers is used. This is due to the requirnment of 2 overlapping kmers.

If k <= core length (if multiple (or 1) kmers are needed to describe the core fully), then the **minimum** score for those kmers is used because only when threshold is set to that minimum or below would the site be called

**Handles both palindromes and non-palindromes.**

If the query is for a palindrome, there is an option isPalindome to set to True. All it does is remove - strand calls, since any match on one strand is a match on another. If this is set to false for a palindrome, it should give + and - matches representing both centered positions, or 2 calls per binding site.

-------
Latest update (1/14/20) - ZM
* New management of .py and .ipynb using utils from https://github.com/kopptr/notebook-tools.
* Minor bug fix

Latest update (1/7/20) - ZM
* Input can be either:
    - Bed coordinates and genome file
    - Fasta file   
* Added additional parameter for scoring tables (optional)
* Fixed a bug where a site would be reported multiple times with different scores

Update (12/5/19) - ZM
* Minor bug fix for when k == core length

Update (11/27/19) - ZM
* Completely redone site calling implementation using numpy arrays.
    - Fixes an issue with overlapping cores, especially at low Escore thresholds. 
    - Now KOMPAS can call overlapping cores in the same manner at any threshold. 
    - Removes a couple dependencies, though they are packeged with python by default
* Fixed a bug with the threshold scoring
* Optimized palindrome mode, should run 2x as fast as non-palindrome mode

Update (11/25/19) - ZM
* Diagnostic tsv file shows all searched peaks, even if no calls were made
* Threshold score is now set up so you can use it as a substitute for running at different thresholds. Represents practical threshold given 2 overlapping, not theoretical. 


Update (11/20/19) - ZM
* Fixed bug in coordinates for 0 vs 1 based. KOMPAS is now fully in-sync with half-open, 0 based formats used by bedtools, MACS, and UCSC (backend)

'''

#nb>
#
################### Parameters ###############
#wkdir = "/Users/ZMielko/Desktop/KOMPAS_Paper/gcPBMAnalysis/gcPBM/"
#kmerFile = '/Users/ZMielko/Desktop/KOMPAS_Paper/gcPBMAnalysis/TFs_Processed/ETS1/ETS1_8mersAligned.txt'
#isFasta = True
## For using coordinate input, can leave as empty string if using fasta
#peakFile = ''
#genomeFile = ''
## For using fasta input, can leave as empty string if using bed files
#fastaFile = '/Users/ZMielko/Desktop/KOMPAS_Paper/gcPBMAnalysis/gcPBM/Ets1.fasta'
## kPosition of the core and where to center the call
## core is right exclusionary, left inclusive [)
#core = (10,14) 
#threshold = 0.4
#isPalindrome = False
## Optional settings
#centerPos = 12 # if set to 'default', will be half position of the core rounded down
#logFile = True
#scoreTable = False # for diagnostic/troubleshooting purposes
## Calculated parameters and error checking
#if centerPos == 'default': # Calculate center position if default
#    centerPos = int((core[0] + core[1])/2)
#    

#py>


import argparse

programDescription = '{K.O.M.P.A.S} Site Caller: calls transcription factor binding sites using PBM data'


parser = argparse.ArgumentParser(description=programDescription)
parser.add_argument('wkdir', type=str, help='Save directory')
parser.add_argument('aKmer', type=str, help='aligned kmer file to use')
parser.add_argument('cStart', type=int, help='Start kPosition of the core, inclusive (int)')
parser.add_argument('cEnd', type=int, help='End kPosition of the core, exclusive (int)')
parser.add_argument('threshold', type=float, help='Escore threshold for calling (float)')
parser.add_argument('-bed', type=str, help='Coordinates to search (.bed)')
parser.add_argument('-genome', type=str, help='Genome file (.fasta/.fa)')
parser.add_argument('-fasta', type=str, help='Fasta file (.fasta/.fa)')
parser.add_argument('-centerPos', type=int, help='Choose a center position based on kPosition, default = midpoint of core rounded down. (int)')
parser.add_argument('-palindrome', action = 'store_true', help='Indicate if the TF is palindromic (True/False)')
parser.add_argument('-logFile',action = 'store_true', help='Generate a log file')
parser.add_argument('-diagnostic',action = 'store_true', help='Generate a score table file for diagnostic purposes')


args = parser.parse_args()

if args.bed is None and args.genome is None and args.fasta is None:
    parser.error("KOMPAS requires either -bed and -genome or -fasta")
if (args.bed and args.genome is None) or (args.genome and args.bed is None):
    parser.error("-bed and -genome need to be specified together")
if (args.fasta and args.genome) or (args.fasta and args.bed):
    parser.error("-fasta argument cannot be used with -bed or -genome")


wkdir = args.wkdir
kmerFile = args.aKmer
isFasta = args.fasta
if args.bed:
    peakFile = args.bed
    genomeFile = args.genome
elif args.fasta:
    isFasta = True
    fastaFile = args.fasta    

core = (args.cStart,args.cEnd) 
threshold = args.threshold
isPalindrome = args.palindrome

logFile = args.logFile
scoreTable = args.diagnostic

if args.centerPos:
    centerPos = args.centerPos 
else:
    centerPos = int((core[0] + core[1])/2)




#>

##### Imports ####

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import reverse_complement
import os
from pybedtools import BedTool

#### Parsing and Conversion Functions ####
def readFasta(file):
    '''
    Reads a fasta file and turns it into a dataframe with forward and rc sequences
    Input: fasta file path
    Output: DataFrame with seq_id, fwd sequence, sequence length, reverse complement
    '''
    entry = []
    with open(file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            entry.append([record.id, str(record.seq).upper(), len(record.seq), str(reverse_complement(record.seq)).upper()])
    arr = np.array(entry)
    df = pd.DataFrame({'seq_id':arr[:, 0], 'fwd':arr[:,1], 'seq_len':arr[:,2], 'rev_comp':arr[:,3]})
    return(df)


def parseID(df):
    """
    Takes the concatinated names given in fasta outputs from bedtool's getfasta 
    and turns them into bed compatible columns
    Input: DataFrame from readFasta
    Output: DataFrame with Chromosome, Start, End columns
    """
    chrom, start, end = [], [], []
    for i in df.seq_id:
        cr = i.split(':')
        pos = cr[1].split('-')
        chrom.append(cr[0])
        start.append(int(pos[0]))
        end.append(int(pos[1]))
    df['Chromosome'] = chrom
    df['Start'] = start
    df['End'] = end
    return(df)

def parseIDFASTA(df):
    """parseIDFASTA takes the fasta ID and converts it into a format comaprable to parseID
    Chromosome = ID, Start = 0, End = len(seq)"""
    chrom, start, end = [], [], []
    for ID, seq in zip(df.seq_id, df.fwd):
        chrom.append(ID)
        start.append(0)
        end.append(len(seq))
    df['Chromosome'] = chrom
    df['Start'] = start
    df['End'] = end
    return(df)

def convertToBed(df, isPalindrome):
    """
    Takes a dataframe with the following columns:
    Chromosome, Start, centerPlus, scorePlus, centerMinus, scoreMinus
    Outputs a dataframe in a bed format 
    """
    chrom, start, orient,scores = [],[],[],[]
    if isPalindrome == False:
        for row in zip(df['Chromosome'],df['Start'],df['centerPlus'],df['scorePlus'],df['centerMinus'],df['scoreMinus']):
            if row[2]:
                for centerP, score in zip(row[2], row[3]): # + sites
                    chrom.append(row[0])
                    start.append(row[1] + centerP)
                    orient.append('+')
                    scores.append(score)
            if row[4]:
                for centerN, score in zip(row[4],row[5]): # - sites
                    chrom.append(row[0])
                    start.append(row[1] + centerN)
                    orient.append('-')
                    scores.append(score)
    else:
        for row in zip(df['Chromosome'],df['Start'],df['centerPlus'],df['scorePlus']):
            if row[2]:
                for centerP, score in zip(row[2], row[3]): # + sites
                    chrom.append(row[0])
                    start.append(row[1] + centerP)
                    orient.append('+')
                    scores.append(score)
    bedDF = pd.DataFrame({'chrom':chrom,'start':start, 'end':start,'orient':orient, 'score':scores})
    bedDF['end'] = bedDF['end'] + 1 # exclusive end position
    return(bedDF)
##### Read in kmer data and process ####
kmer = pd.read_csv(kmerFile, sep = '\t')
k = len(kmer['kmer'][0])
coreLen = core[1] - core[0]
# Find the kPositions required, any would be sufficient to call
if k > coreLen: 
    searchEnd = core[1]
    checkK = 0
    ReqKpos = set() #
    while checkK != core[0]:
        checkK = searchEnd - k
        if checkK <= core[0]:
            ReqKpos.add(checkK)
            searchEnd = searchEnd + 1
# Or find the group of all kPositions that are needed, all or none
else:
    searchStart = core[0]
    checkK = 0
    ReqKpos = set()
    while searchStart + k <= core[1]:
        ReqKpos.add(searchStart)
        searchStart = searchStart + 1
# Determine flanks of ReqKPos for threshold score reporting
ScoredKpos = ReqKpos.copy()
if k >= coreLen:
    ScoredKpos.add(min(ReqKpos) - 1)
    ScoredKpos.add(max(ReqKpos) + 1)        

# Generate dictionary for quick retreaval of relevant kmers
thrKmers = kmer[(kmer['Escore'] > threshold) & (kmer['kposition'].isin(ScoredKpos))]
kDict = dict(zip(thrKmers['kmer'],zip(thrKmers['kposition'],thrKmers['Escore'])))


##### Generate a dataframe from the input sequence information ####
if isFasta == True: # if given a FASTA file of sequences    
    peakDF = readFasta(fastaFile)
    peakDF = parseIDFASTA(peakDF)
    del peakDF['seq_id']
    peakDF[['seq_len', 'Start','End']] = peakDF[['seq_len', 'Start','End']].astype(int) 
else: # if given a bed and genome file
    peaks = pd.read_csv(peakFile, sep = '\t',header = None, usecols=[0,1,2])
    # find out how long the peaks need to be to use them
    if coreLen > k:
        minLen = coreLen
    else:
        minLen = k + 1
    peaks = peaks[peaks[2]-peaks[1] > (minLen)].drop_duplicates() # filter for short sequences the caller would have trouble with
    # Generate peaks using bedtools
    peakBed = BedTool.from_dataframe(peaks)
    peakSequence = peakBed.sequence(fi = genomeFile)
    peakDF = readFasta(peakSequence.seqfn)
    peakDF = parseID(peakDF)
    del peakDF['seq_id']
    peakDF[['seq_len', 'Start','End']] = peakDF[['seq_len', 'Start','End']].astype(int) 


##### KOMPAS Scoring and Calling Functions #####   
    
def kmerMatch(seq):
    """
    Returns matched positions in the sequence and their kpositions
    Input: sequence
    Output: consecutive positions, kpositions, and scores above threshold
    """
    # Get the kposition and kscore for the peak, save a numpy array
    kpos,kscore = [], []
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        if window in kDict:
            kpos.append(kDict[window][0])
            kscore.append(kDict[window][1])
        else:
            kpos.append(0)
            kscore.append(-0.5)
    kpos = np.array(kpos)
    kscore = np.array(kscore)
    # Get consecutive positions, kpositions, and score via numpy operations
    if k >= coreLen:
        position = list(filter(lambda x: len(x) != 1,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
        kpos = list(filter(lambda x: len(x) != 1,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
    elif k < coreLen:
        reqLen = len(ReqKpos)
        position = list(filter(lambda x: len(x) == reqLen,np.split(np.r_[:len(kpos)], np.where(np.diff(kpos) != 1)[0]+1)))
        kpos = list(filter(lambda x: len(x) == reqLen,np.split(kpos, np.where(np.diff(kpos) != 1)[0]+1)))
    kScore = []
    for pos in position:
        kScore.append(kscore[pos])
    return(zip(position, kpos, kScore))

def findCenter(zippedMatch, orient, seqLen):
    """
    Given a zip of match position, kposition, and kscore
    Returns the center sites and threshold kscore
    """
    centerSites = []
    siteScores = []
    for pos, kpos, kScore in zippedMatch:
        centerSite = (centerPos - kpos[0]) + pos[0]
        if orient == 'rc':
            centerSite = (seqLen - centerSite) -1
        centerSites.append(centerSite)
        if k >= coreLen:
            score = threshold
            for score1, score2 in zip(kScore, kScore[1:]):
                caniScore = sorted([score1, score2])[0]
                if caniScore > score:
                    score = caniScore
            siteScores.append(score)
        elif k < coreLen:
            siteScores.append(min(kScore))
    return(pd.Series([centerSites, siteScores]))
             
###### Run Calling Functions on Sequence Data ######
if isPalindrome == True:
    peakDF[["centerPlus","scorePlus"]] = peakDF.apply(lambda peakDF: findCenter(kmerMatch(peakDF['fwd']), 'fwd', len(peakDF['fwd'])), axis = 1) 
    if scoreTable == True:
        peakDF.to_csv(f'{wkdir}CallTable_{threshold}_KOMPAS.tsv', sep = '\t')
    peakDF = peakDF[~peakDF['centerPlus'].isna()]
else:
    peakDF[["centerPlus","scorePlus"]] = peakDF.apply(lambda peakDF: findCenter(kmerMatch(peakDF['fwd']), 'fwd',len(peakDF['fwd'])), axis = 1) 
    peakDF[["centerMinus","scoreMinus"]] = peakDF.apply(lambda peakDF: findCenter(kmerMatch(peakDF['rev_comp']), 'rc',len(peakDF['rev_comp'])), axis = 1)
    if scoreTable == True:
        peakDF.to_csv(f'{wkdir}CallTable_{threshold}_KOMPAS.tsv', sep = '\t')
    peakDF = peakDF[~peakDF['centerPlus'].isna() | ~peakDF['centerMinus'].isna()]
    
# Convert to Bed format and save
finalBed = convertToBed(peakDF, isPalindrome)
finalBed = finalBed.sort_values(by='score', ascending = False).drop_duplicates(['chrom', 'start', 'end', 'orient'],keep = 'first')
finalBed.to_csv(f'{wkdir}Centered_{threshold}_KOMPAS.bed', sep = '\t', header = None, index = False)

####### Optional Output #######

# Score Table Output
if scoreTable == True:
    def ScoreTableGenerate(seq):
        """
        Returns matched positions in the sequence and their kpositions
        Input: sequence 
        Output: kpositions, and scores above threshold
        """
        # Get the kposition and kscore for the peak, save a numpy array
        kpos,kscore = [], []
        for i in range(len(seq) - k + 1):
            window = seq[i:i+k]
            if window in kDict:
                kpos.append(kDict[window][0])
                kscore.append(kDict[window][1])
            else:
                kpos.append(0)
                kscore.append(-0.5)
        kpos = kpos
        kscore = kscore
        return(kpos, kscore)

    fwdScores = peakDF.apply(lambda peakDF: kmerMatchScoreTable(peakDF['fwd']), axis = 1) 
    rcScores = peakDF.apply(lambda peakDF: kmerMatchScoreTable(peakDF['rev_comp']), axis = 1) 
    pd.DataFrame({'fwd':fwdScores, 'rc':rcScores}).to_csv(f'{wkdir}ScoreTable_{threshold}.tsv', sep = '\t')

# Log File
if logFile == True:
    ##################
    # Log the output #
    ##################
    f = open(wkdir + "/KOMPASLog.txt", "a")
    f.write("##### Parameters ##### \n")
    if isFasta == True:
        f.write(f"Fasta file: {fastaFile}"+ "\n")
    else:
        f.write(f"Peak file: {peakFile}"+ "\n")
        f.write(f"Genome file: {genomeFile}"+ "\n")
    f.write(f"kmer file: {kmerFile}"+ "\n")
    f.write(f"Core kPositions: {core}"+ "\n")
    f.write(f"Center kPositions: {centerPos}"+ "\n")
    f.write(f"Threshold: {threshold}"+ "\n")
    f.write("#Summary# \n")
    f.write(f"Final # of called sequences: {len(finalBed)}"+ "\n")
    f.close()

