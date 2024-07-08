#!/usr/bin/env python
## KMerCountRabbit ('KCR')-- Very Simple and Quick Reference-based read and kmer counting for NGS Datasets
##
## Version ak08 07_04_2024
##
## OVERVIEW OF KCR
## - You have a bunch of datasets from next generation sequencing and a reference, and want to count matches between them
## - KCR is a tool providing tables of read match counts between sequence datan(fasta or fastq read files) and reference genome(s)
## - KCR has some shared capabilities with standard tools like BLAST/Bowtie/Kallisto/BWA/Kraken and to many home-grown tools
## - The goals with KCR are to provide a potentially simpler worflow for analysis
## -   What goes in: The names of your reference file and the location of sequencing datafiles
## -   What comes out: A table of match counts (Kmers and/or Reads) and (for preview purposes) a very quick visualization
## -   Important Features: Quick, Simple, Consistent handling of sequences shared between references 
## - KCR is K-mer based
## -   Any read that matches at least one k-mer in a reference genome is assigned to that genome
## -   If only one reference genome is matched, a k-mer or read is considered to be "uniquely" mapped
## -   Reads (or K-mers) that match multiple references are considered to be "public"
## -   A third (optional) category called "local repeat" describes kmers present multiple times, but all in one reference sequence
## -   With the default K-mer length of 31, the specificity of matches is such that false positive should be vanishingly rare
## - KCR is primarily a tool for counting matches (not for many other tasks):
## -   KCR is not designed to provided strandedness data, detect distant homology, or assemble sequencing data into contigs.
## -   KCR is not designed for (or helpful in) detecting mutations, rearrangements, or termini
## -   KCR provides a provisional (and customizable) graphic output, but this graphic is intended for
## -     quick experimental examination, while publication may require additional tools.
## - Specific features of KCR overlap with local tools REVA,HQAlign,Polybench, Jazz18Heap, and many public tools
## - If everything is working properly, KCR can be much faster to implement and/or run than any of these
## - KCR is very much a work in progress- not recommended for production until we've tested it more thoroughly
## - Report any crashes or anomalous results to (afire@stanford.edu) 
##
## BASIC OPERATION
## ->Syntax
##  - python KMerCountRabbit<ver>.py  Reference=<file,file..> Data=<file,file..>  
##  -   you will need the module "VSG_ModuleFO.py" in the same folder (VSG_ModuleFN.py is on GITHUB Site)
##  - For the program to run as fast as possible, the PyPy interpreter should be used (www.pypy.org)
## ->Required Parameter 1: A Reference file with one or many 'reference' sequences in fasta format
##  - Reference=<file[s]> : Datasets to gather kmers for analyzing data
##  -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
##  -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
## ->Required Parameter 2: Data Files and Sample Names
##  - Direct FileName Option: Data=<fasta/q file[s]> : Here you input one or more sequencing data files for analysis
##  -   This can be a single file or a comma-delimited list and can include wildcards (*)
##  -   "Data" can also be a directory (in this case, all fasta/fastq/fastq.gz/fastq.gz in the directory will be analyzed).
##  -   By default the sample name will be the filename and Rabbit "AutoPair" R1/R2 data file pairs as a single ('Rx') entity   
##  - SampleToData MetaFile Option: For more flexibility, provide SampleName/ReadFileNames combinations in your own file
##  -   This is invoked with SampleToData=<your SampleDescriptorFile> instead of Data= in command line
##  -   SampleDescriptor files have simple line-by-line format where each line consists of
##  -   A sample name, followed by tab, followed bya comma-delimited list of relevant the file paths with data from that sample
##  -
## ->Some additional parameters to explore as you evaluate KCR as a tool
##  -    (all of these [and other] parameters have defaults, so you can ignore all of these and just run the program) 
##  - klen : is the Kmer length (default 31, realistically can be between 16 and read-length; values of 31 will slow this down). Default=31
##  - OutFileName : name assigned by program if not specified> : Allows user to specify base filenames for data output. Default=<autoassign> 
##  - MaxReadsPerFile : Max number to check for each Data file (default is zero (all reads).  Default=0 (process all reads)
##  - CircularReference : Setting this to true will instruct KFR to treat reference sequences as circles.  Default=False (treat as linear)
##  - ReferenceFileByFile : Setting to true treats each reference file as a single entity.  (Default=Each fasta entry is an entity)
##  - Express : Setting this to True runs the program in mode which just yields first-matching-reference read counts (Default=False)
##  - FilterList : Providing a fasta file and a cutoff number (e.g. myFilterFile.fa>5) limits the searched kmers
##  -   based on kmer copy number in a separate reference file
## -> Many other settings and options are described in the detailed instructions (start of script).
##  - If you have several different input files as data and expect the job to take more than a few minutes, you can run in multithreaded mode
##  -   Setting Threaded=True instructs rabbit to decide how many threads to use).  Of you can specify a number (e.g. Threads=3 uses three threads)
## -> Some of KCR is still being debugged-- so caveat mappor
## End Help

## User settings

## Data Input - Generally a set of .fastq or .fasta files resulting from sequencing experiments.  Can be gzipped (.gz)
#### The list of datafiles can be a simple (comma-delimited) list or can involve a wildcard (e.g. '*.fastq.gz').
#### Any list or wildcard-based specification needs to be quoted (the entire list enclosed in a single quoted block)
DataFiles1 = '' ## (data=) The input data file list.  The list should be enclosed by quotes (e.g. 'Run1.fastq,Run2.fastq,Run3.fastq')
SampleToData1 = '' ## (sampletodata=) An alternative to a simple input data file or file list.  A file here specifies correspondence between samples and filenames.  Format: Each line should have SampleName <tab> comma-delimmited list of filenames
NoNs1 = False ## ( nons=)Setting this to true ignores any read with "N" bases 
AutoPair1 = True ## (autopair=) Setting this to true attempts to pair read1 and read2 files
R1Trim5 = 0 ## (r1trim5=) Length to trim from 5' end in R1
R2Trim5 = 0 ## (r2trim5=) Length to trim from 5' end in R2
R1Trim3 = 0 ## (r1trim3=) Length to trim from 3' end in R1
R2Trim3 = 0 ## (r2trim3=) Length to trim from 3' end in R2
MaxReadLen1 = 0 ## (maxreadlength=) After filtering, reads will be truncated to this length. ==0 for no maximum read length
EndLinker1 = '' ## (endlinker=) A linker at the end of the read will be chewed off using an rfind command (for short RNA reads one such is 'TGGAAT')
RequireEndLinker1 = False ## (requireendlinker=) Setting this to true removes any read without this linker
MaxReadsPerFile1 = 0 ## Set to a positive number to limit the number of reads per file

## Reference Input
#### The reference list is generally a .fasta file or list of files (can be gzipped) with reference sequence to match against.
ReferenceFiles1 = '' ## (reference=)(ref=) List of reference files
CircularReference1 = False  ## Setting to True will compare reads to a circular reference
ReferenceFileByFile1 = False ## Setting this to True will aggregate each input reference file into a single entity

## Settings for sequence comparison
klen1 = 31 ##(klen=) Kmer length for counting matched subsequences 
ReportGranularity1 = 1000000 ## How frequently to report updates
NoDup1 = False ## (nodup=) Setting NoDup to True instructs Rabbit to not to count k-mers in a paired end read that are also in R1.  This provides somewhat more nuanced K-mer counts but will slow the program down substantially; Read counts are not affected.

## Settings for has-based pre-screen (these don't affect the final output but can speed up the program in some cases)
TurboRareK1 = 'Auto' ## A KMer length to be used for initial prescreening.  This is particularly useful for datasets where most reads fail to match the reference
TurboRareLimit1 = 0 ## Setting this to a positive number n insists that there be a matching number within the first n bases of any read for it to be considered 
EstimatedAverageReadLen1 = 300 ## (expectedaveragereadlen=) This only affects the speed-optimization (by guiding the choice of turboKlen when it is set to "Auto"). Can change this if longer or shorter reads are expected.


## Settings to avoid counting of k-mers that match a preset file of sequencing adaptors
TetritisF1 = '' ## (tetritisf=) Name of tetritis file for optional on-the-fly masking of linker bits (file illuminatetritis1223wMultiN.fa, distributed with this program, can be used)
                ## Note that the distributed file has many of the publically available Illumina linkers as well as homopolymers of the four nucleotides
TetritisK1 = 16 ## (tetritisk=) Length of K-mer for tetritis filtering

## Settings for detailed handling of multimapping/repeated sequences
FirstDibs1 = False ## (firstdibs=) Setting this to true assigns each k-mer and read that is shared to only the first matching reference (disallows multimapping)
PrivateOnly1 = False ## (privateonly=) Setting this to true ignores any k-mer present in more than one reference
AggregateUniqueAndLocalRepeats1 = False ## (aggregateuniqueandlocalrepeats=) Setting this to true provides a single output that sums both unique and local repeats

## Setting to filter the index of reference K-mers used to generate counts
DustFilter1 = False ## Setting this to true implements a simple filter to avoid low complexity search kmers
MinRefCov1 = 0 ## setting this to a positive floating point ignores all reference segments where the coverage reported by MegaHit or Spades is lower than this value
MinRefLen1 = 0 ## setting this to a positive integer ignores all reference segments that are smaller than this value
FilterList1 = '' ## This can be used to filter k-mers based on their copy number in an arbitrary fasta file (or group of files). Adding FilterList='MyFile.fa>5' to the command line will only use kmers with >6 instances in MyFile.fa.
                 ## Using the syntax FilterList1=Reference>5 filters k-mers based on their copy number in the Reference files (in this case only counting reads present >5 times)
                 ## Multiple files can be separated by commas, and multiple filters seprated by semicolons
                 ## Thus as a complex example FilterList='reference>5;myFileA.fa>3;myFileB,myFileC<2' will keep and analyze only kmers that meet the criteria
                 ##    More than 5 copies of the kmer the aggregate reference files specified in ReferenceFiles=
                 ##    More than 3 copies of the together in myFileA.fa file
                 ##    Less than 2 copies (total) in the files myFileB myFileC
## KMerThinning
##   For very large reference datasets it may be impractical to store all Kmers in memory, even in condensed form
##   In such cases, a "IndexCadence" can be set, which will use K-mers periodically sampled from the reference sequence
##   As an example, setting IndexCadence=8 will use only every 8th k-mer in searching for matches.
##   This option can lose some ability to distinguish reads from highly related reference sequences, and could miss some reads
##   But for large reference sets, it should provide a consistent and reflective set of metrics with a much small memory footprint
##   Note that KMer thinning is of most benefit where many reads match the reference and the reference is a complex set of sequences (many Kmers)
##   Otherwise the standard TurboK filtering should be fine.
##   And note that Kmer thinning and TurboK filtering don't synergize (the two together would rarely be better than either alone)
ReferenceCadence1 = 0 ## Setting this to a value of n>1 will limit the k-mers index from the reference to only every x'th k-mer 
DataCadence1 = 0 ## Setting this to a value of n>1 will limit the k-mers searched in the data to only every x'th k-mer
## KMER Recusal Filtering
## KMER Recusal filtering provides an approach toward sensitive and accurate identification of sample-to-sample count differences
## kmers that have an usually high or low prevelance can greatly skew relative counts and coverage respectively,
## This can potentially be problematic in (i) yielding spurious sample-to-sample differences and (ii) masking real differences
## Such k-mers can come from matching to other sequences not in the relevant reference, or can come from artefacts of asssembly or sequencing bias
## To provide counts that avoid such confounding factors, Rabbit has a number of options to "recuse" K-mers based on frequency
## In general the recusal filters work by obtaining a counts-per kmer distribution for each reference segment in the actual data
## Once these simple K-mer counts are obtained, a filtering step recuses K-mers that are outside of a range that is set by the user
## Two options for setting thresholds for  are provided:
##   i. Rank-based recusal
##    The user can set fractional thresholds based on the actual distribution of kmer counts for each reference
##     The relevant parameters are ColdKMerF1 and HotKMerF1
##     Setting ColdKMerF1 to a non-zero value x will recuse all kmers in the bottom x fraction as sorted by sample counts
##       So an example: Setting ColdKMerF1 = 0.2 will recuse the bottom 20% of K-mers by total count from each reference segment
##     Setting HotKMerF1 to a value <1 recuses all but the indicated fraction of kmers sorted by frequency in samples
##       So setting HotKMerF1=0.8 will recuse the top 20% of kmers for each reference based on count
##   ii. Count-relative-to-median based recusal
##    The user can set a threshold that will be based on median k-mer counts for each reference
##     Setting ColdKRelToMedian1 to a value (generally <1) recuses all K-mers whose counts are below ColdKRelToMedian*MedianKCoverage for that reference
##     Setting HotKRelToMedian1 to a value (generally >1) recuses all K-mers whose counts are above HotKRelToMedian*MedianKCoverage for that reference
## Note that gathering the recusal list entails a complete analysis of the data to obtain counts and identify k-mers for recusal
##   The recusal run is also done in single-core mode, so this will certainly slow the program down overall.
##      Note that you can save the recused K-mer list with the WriteIndex=<myIndexFileName> option
##      (and use this recused index in a future run with ReadIndex=<myIndexFileName> option)
HotKMerF1 = 0    ## Setting HotKMerF1<1 instructs Rabbit to recuse any kmer above this fractional rank amongst most abundant K-mers in each reference sequence
                 ## HotKMerF1 is of particular utility in cases where a small number of k-mers in a reference are mapping frequently to some other related (or unrelated) sequence
                 ## Example: Setting HotKMerF1=0.98 will ignore the most abundant 2% of k-mers in every reference   
ColdKMerF1 = 0   ## Sets a floor on Kmer counts based on rank 
HotKRelToMedian1 = 0  ## Sets a ceiling on Kmer counts based on median count over all datasets for kmers in that reference sequece
ColdKRelToMedian1 = 0 ## Sets a floor on Kmer counts based on median count over all datasets for kmers in that reference sequece
                   
## Settings for output
OutFileName1 = 'default' ## (outfile=)(outfilebase=)
SkipBlankRefs1 = True ## (skipblankrefs=) Setting this to True will mask output for reference sequences that are not detected
SkipBlankData1 = True ## (skipblankdata=) Setting this to True will mask output for datasets where no reference match is detected
Delimiter1 = '\r' ## (delimeter=) New line delimiter for output

## Some additional user-specifiable values
Express1 = False ## (express=) Setting this to True provides an express version of the program yielding only read counts and with ambiguous reads assigned to first match
FastQDumpProgram1 = 'fasterq-dump' ## Full path of fastq-dump/fasterq-dump program to download from SRA if needed.
PrimaryIndex1 = 'default' ## (primaryindex=) (index=) (readindex=) Other settings can offer speedup or memory management advantages for larger files.  PrimaryIndex can be from "Data" or "Reference" or a file name with a previous index
WriteIndex1 = '' ## Setting this to a filename will write the index that is compiled so it can be used for a future run.  For large reference files (>100s of MB) using this file can save time in future runs
Threads1 = 0 ## (threaded=) (threading=) Number of threads for multiprocessing (zero for single core execution, "Auto" or "True" to set this automatically, or provide a specific integer of how many cores to use)
MaxMemPerThread1 = 0.5 ## Maximum fraction of total memory that will be allocated per Thread (default=0.5-- half of total)

## Some settings for graphic and .tdv output
vscale1 = 14 ## (vscale=) Vertical scale for graphic
hscale1 = 42 ## (hscale=) Horizontal scale for graphic
MinCountsForDisplay1 = 1 ## (mincountsfordisplay=) References with fewer average reads-per-sample will be skipped in display
SortReferencesByReadCount1 = False ## (sortreferencesbyreadcount=) Setting this to true will sort references from low to higher counts
ClusterReferences1 = False ## (clusterferences=) Setting clusterferences to true implements a very rudimentary references clustering for the graphic display 
MaxRefToDisplay1 = 0 ## (maxreftodisplay=) Limits the number of references displayed on the illustrative graph (zero is unlimited).  Priority for display is based on kmer counts (read counts in express mode)
ClusterSamples1 = False ## (clusterferences=) Setting clustersamples to true implements a very rudimentary samples clustering for the graphic display 
GraphTitle1 = 'default' ## (graphtitle=) Sets a graph title.  Default will be OutputFile Name
SuppressGraphDetails1 = False ## setting this to true suppresses the detailed listing of program parameters on the output graphic
## Choice of output fields
#     'KmerCount'               :  7   ## Raw Kmer Counts 
#     'KmerCount2'              :  8   ## Sum of Raw Kmer counts squared
#     'ReadCount'               :  9   ## Read Counts
#     'KmersCovered'            : 10   ## Number of kmers covered
#     'KmerDepthAve'            : 11   ## Average depth of coverage 
#     'KmerDepthStdDev'         : 12   ## Standard deviation of coverage depth
#     'KmerCoverage'            : 13   ## Fraction of Kmers covered in data
#     'KmerSampleFraction'      : 14   ## Fraction of pass-filter kmers in sample assigned to this partition
#     'ReadSampleFraction'      : 15   ## Fraction of pass-filter reads in sample assigned to this partition
#     'KPKM'                    : 16   ## Kmers mapped per thousand reference kmers per million pass-filter kmers
#     'RPKM'                    : 17   ## Reads mapped per thousand reference kmers per million pass-filter reads
#     'KmerEnrichment'          : 18   ## Enrichment of kmers mapped to this partition over observed for all datasets
#     'ReadEnrichment'          : 19   ## Enrichment of reads mapped to this partition over observed for all datasets

HeatMapColorField1 = 'default'  ## (heatmapcolorfield=) Sets the raw or processed value used to color rectangles in the summary heatmap.  See options above
TDVFields1 = 'default'      ## (tdvfield=) (tdvfields=) Specifying one or more of the "processed" measurement values above (#11-19 above) in this list expands the .tdv output to include these 
HeatMapTextFields1 = 'default' ## (heatmaptextfield=) (heatmaptextfields=) Specifying one or more values in this list will add the relevant values to the heatmap
## Note that KMer counts are not calculated in Express mode, so any kmer-reliant value will be recorded as zero in Express mode
## If you encounter this either run with Express=False or choose read-count (or RPKM or ReadSampleFraction) as the displayed value
## ->Here are a few more details
##  -
##  - Other Notes
##  - Any version of python beyond 3.7 can be used.  For speed, Python 3.11 (or better, PyPy3.9+) is highly recommended
##  - For SRR/ERR/DRR datasets not present locally, KFR will try to download the files from NCBI (linux/mac only).
##  -   You can provide the location of fastq-dump with FastQDumpProgram= (otherwise, KFR will look for it, which can be slow)
##  - Memory and multiprocess management
##  -     i. You can set a maximum number of CPU threads for the analysis with the command line option Threads=
##  -     ii. If Reference complexity is extensive (>500M different kmers), while Data is a smaller number (e.g. small set of reads), 
##  -       then setting PrimaryIndex = 'Reference' can avoid memory issues (avoids preassembling a preliminary kmer Set from Datas)
## End list of user settings

##  -
##  DESCRIPTION OF MATCHING CRITERIA USED TO GENERATE READ COUNTS FOR EACH REFERENCE
##  - The following discussion is a detailed description of the sausage making involved in assigning reads.
##  - With a long-ish k-mer  (e.g. the default, which is 31) and reference sequences that are not highly degererate, the following considerations are not substantial
##  - For pools of references that are very similar to each other, the discussion becomes more relevant
##  
##  - The simplest category of Kmers are those that are present just once and in one reference ('Unique' k-mers)
##  - Any single k-mer shared between a read and a reference nominates the read as matching that reference
##  - In cases where there is no ambiguity (ie no match to any other reference from the read), the nomination is counted as a "private match"
##  -    This means that a read with any number of k-mer matches to a reference is counted as a match to the reference
##  - Reads that match unique "private" k-mers from more than one reference are assigned arbitrarily to the last matching reference in the list
##
##  - Kmers present in more than one reference are a second category and are called 'public'
##  - Reads with 'public' kmers can still be assigned to an individual reference if there are 'private' k-mers matching a unique reference
##  - In cases where all kmers in a read match only public k-mers, the read is assigned to the public category
##  - A public-category read will still be assigned to one or more references; so each reference can have its own "public" category
##  - Assignment of a public read is based on the most restricted Kmer present in that read (kmer present in fewest references)
##  - In standard mode, public reads are assigned to every reference matching the assigning k-mer.
##  -     Turning on the option FirstDibs=True alters this so that only the first reference with the assigning k-mer is credited
##  - Turning on the option PrivateOnly=True tells KCR to ignore public kmers and reads (focusing only on k-mers belonging to a single reference)
##  -
##  - A third category of k-mers ("LocalRepeat" Kmers) are those present more than once but only in a single reference.
##  -     These are reported as a separate category in Kmer counts
##  -     Reads that only match LocalRepeat Kmers are reported also as in the LocalRepeat reads category
##  -     Reference Assignment based on a unique k-mer takes precedence over assignment based on a local repeat k-mer
##  - Setting the parameter AggregateUniqueAndLocalRepeats=True instructs KCR to group local repeat and unique k-mers together
##  -     The combined Unique+LocalRepeat category is then designated "Private"
##  -
##  - Kmer counts are provided along with read counts and can be used to evaluate coverage where there is complex redundancy between references
##  - Where references are well distinguished from each other (few or no shared k-mers), the "Express" setting can be used
##  -      Express=True tells KCR to assign each read based on the first uniquely matching k-mer encountered (no need to go through the rest of the read)
##  -      Express mode doesn't count k-mers -- just reads; for datasets where most reads match the reference this can speed the program up greatly
##  -      Cases where only a small portion of reads match the reference won't be sped up in express mode

from VSG_ModuleFO import *
vcommand()
import sys,os,gzip,pickle
from collections import Counter
from array import array
from glob import glob
from itertools import zip_longest
from math import ceil
if not('pypy' in sys.version.lower()):
    vLog('\n*** ADVISORY *** KCR running in Standard Python: will be faster if run in PyPy (see PyPy.org)\n')
def myOpen1(fn,mode):
    vLog('Opening File:',os.path.basename(fn), "  mode='"+mode+"'",'vsilent')
    if fn.lower().endswith('.gz'):
        return gzip.open(fn, mode=mode)
    else:
        return open(fn, mode=mode)
def myClose1(f):
    vLog('Closing File: ',f.name,'vsilent')
    f.close()

## C1 Values represent whether a k-mer is observed in no references, several references, mutiple times in one reference, or uniquely
## 0 not observed
## -1 to -RefCount : Observed in multiple references with lowest observation in reference RefCount [RefCount+C1[x]]
## -RefCount-1 to -2*RefCount: Observed multiple times in one reference [2*RefCount+C1[x]
## -2*RefCount-1 to -3*RefCount: Observed uniquely in one reference [3*RefCount+C1[x]]
## A read will be labeled with the smallest matched kmer

mask1 = 4**(klen1-1)-1
ksam1 = 4**(klen1-1)
Tmask1 = 4**(TetritisK1-1)-1
Tmask0 = 4**TetritisK1-1
Tsam1 = 4**(TetritisK1-1)
BaseD1 = Counter({'G':0, 'A':1, 'T':2, 'C':3, 'N':0, 'g':0, 'a':1, 't':2, 'c':3, 'n':0, 'U':2, 'u':2})
BaseL1 = [0]*256
BaseA1 = [0]*256
for b1 in BaseD1:
    BaseL1[ord(b1)] = BaseD1[b1]
    BaseA1[ord(b1)] = (3-BaseD1[b1])*ksam1
def AreWeAPair1(myf1,myf2):
    '''Returns a joint f1+f2 name if myf1 and myf2 seem to be a Read1/Read2 pair, false if they are not .  With thanks to S.Sondheim'''
    if len(myf1)!=len(myf2):
        return False
    PairName = ''
    for i,(mc1,mc2) in enumerate(zip(myf1,myf2)):
        if mc1==mc2:
            PairName += mc1
            continue
        else:
            PairName += 'x'
        if (mc1,mc2)!=('1','2'): return False
        if i<1 or i+1==len(myf1): return False
        if myf1[i-1] in '_-. ' and myf1[i+1] in '_-. ': continue
        if myf1[i-1] in 'rR' and (not(myf1[i-2].isalpha()) or myf1[i-2].isupper()!=myf1[i-1].isupper()) and not(myf1[i+1].isdigit()): continue
        return False
    return PairName            
def TruncName1(n):
    if not(type(n))==str:
        newN = ''
        for cc1,cc2 in zip_longest(n[0],n[1]):
            if cc1==cc2:
                newN += cc1
            else:
                newN += 'x'
        n = newN
    return os.path.basename(n).split('.gz')[0].split('.fastq')[0].split('.fasta')[0]
def SimpleDust1(k,SubKLen=(3,6,9),SimpleFilter=(0.18, 0.05,0.02)):
    ## k is any k-mer (e.g., 32-mer), SubKLen is a list of lengths to do dust filtering at (e.g., 3,6,9), SimpleFilter is a set of thresholds (coincidence frequencies) above which a value of True  will be returned
    ## The function calculates how many sub-k=mers of each length match other sub-k-mers (total coincidences as a function of the maximum possible number
    ## a homopolymer k-mer k will have coincidences==possible-coincidences and c/mc will be 1.0.  A sequence with no sub-k-mer appearing more than once will have a c/mc of 0.0
    ## sequences of intermediate complexity will have an intermediate c/mc
    ## Heuristic exploration leads to a provisional recommendation of subkmer lengths of SubKLen=(3,6,9) for 32-mer words, with thresholds >0.18, >0.05, >0.02 for calling a given k-mer degenerate
    for z in range(len(SubKLen)):
        dustk = SubKLen[z]
        dustmask = 4**dustk-1
        c = Counter()
        n = klen1-dustk+1 ## number of SubKmers in each Kmer
        for j in range(n):
            cd1 = (k>>(2*j))&dustmask
            c[cd1]+=1
        v = sum([x**2 for x in c.values()])-n ## coincidences
        mc = n**2-n ## possible coincidences
        if v>SimpleFilter[z]*mc:
            return True
    return False

def CountWorker1(f1=0,sn='fromfeed',fL1='fromfeed',tools='',TurboS1=''):  #f1 is the sample number [0 based], sn is the sample name, fL1 is the file list, tools are inherited values for single thread calculation
    if tools: ## single core mode
        C1,C2,ReferenceFileList1,RefCount,ReadNumC,TetTrim,TSet,TurboRareK,Threads = tools
        myTask = f1,fL1,sn
    else:  ## multi core mode
        C1,C2,ReferenceFileList1,RefCount,ReadNumC,TetTrim,TSet,TurboRareK,Threads = vGrab("BasicData")
        myTask = vEat()
    maskS0 = 2**(2*TurboRareK)-1
    maskS1 = 2**(2*TurboRareK-3)-1
    shiftS1 = 2*TurboRareK-3
    if TurboRareK and Threads:
        TurboS1 = array('B',[0]*(2**(2*TurboRareK-3)))
        for F1 in ReferenceFileList1:
            for n11,s11 in vFastAToDict(F1,upper=True).items():
                if MinRefCov1:
                    if 'multi=' in n11:
                        myMulti1 = float(n11.split('multi=')[1].split()[0])
                        if myMulti1<MinRefCov1: continue
                    if '_cov_' in n11:
                        myLen1 = float(n11.split('_cov_')[1].split('_')[0])
                        if myLen1<MinRefCov1: continue
                if MinRefLen1:
                    if 'len=' in n11:
                        myLen1 = float(n11.split('len=')[1].split()[0])
                        if myLen1<MinRefLen1: continue
                    if '_length_' in n11:
                        myLen1 = float(n11.split('_length_')[1].split('_')[0])
                        if myLen1<MinRefLen1: continue
                    if len(s11)<MinRefLen1: continue
                if CircularReference1: s11 += s11[:klen1-1]
                a11 = vantisense(s11)
                for p1 in range(len(s11)-TurboRareK+1):
                    myHash = hash(s11[p1:p1+TurboRareK])&maskS0
                    myIndex = myHash & maskS1
                    myBit = myHash >> shiftS1
                    TurboS1[myIndex] |= 1<<myBit
                    myHash = hash(a11[p1:p1+TurboRareK])&maskS0
                    myIndex = myHash & maskS1
                    myBit = myHash >> shiftS1
                    TurboS1[myIndex] |= 1<<myBit
    while myTask:
        f1,fL1,sn = myTask
        f1 = int(f1)
        if type(fL1)==str:
            fL1 = [fL1,]
        myReadL = [0 for y in range(RefCount*3+1)]
        myKmerL = [0 for y in range(RefCount*3+1)]
        myKmerC = Counter()
        myReadC = Counter()
        TotalReads = 0 ## Includes reads that were not filtered out
        TotalKmers = 0 ## Includes reads that were not filtered out
        PassFilterReads = 0
        PassFilterKmers = 0 ## Includes kmers from reads that were not filtered out 
        for F1  in fL1:
            if type(F1)==tuple:
                Paired = True
                F2 = F1[1]
                F1 = F1[0]
            else:
                F2 = None
                Paired = False
            matchCount = 0
            if len(fL1)>1:
                myfn1 = 'File:'+TruncName1(F1)+' '  ## used to narrate filename if there are multiple files per sample
            else:
                myfn1 = ''
            DataCycle1 = 2 ## 2 for fasta, 4 for fastq
            if F1.lower().endswith('.fastq') or F1.lower().endswith('.fastq.gz'): DataCycle1 = 4
            LastLine1 = DataCycle1*(MaxReadsPerFile1-1)
            if  not(Paired) and ReadNumC[F1]==2:
                Trim5 = R2Trim5
                Trim3 = R2Trim3
            else:
                Trim5 = R1Trim5
                Trim3 = R1Trim3            
            v0 = 0
            a0 = 0
            ReportLineGranularity1 = DataCycle1*ReportGranularity1
            try:
                if Paired:
                    myIter = zip_longest(myOpen1(F1, mode='rt'),myOpen1(F2, mode='rt'))
                else:
                    myIter = myOpen1(F1, mode='rt')
                i1 = 0
                for i1,L1 in enumerate(myIter):
                    if i1%DataCycle1==1:
                        if MaxReadsPerFile1 and i1>LastLine1+1: break
                        if i1>1 and i1%ReportLineGranularity1==1:
                            vLog('Counting Progress ' +myfn1+'Sample:'+sn+' Reads:'+str(i1//DataCycle1)+' Matched:'+str(matchCount))
                        TotalReads += 1
                        if Paired:
                            L2 = L1[1].strip()
                            L1 = L1[0].strip()
                            TotalKmers += len(L1)+len(L2)-2*klen1+2
                            if NoNs1 and (('N' in L2) or ('N' in L1)): continue
                            if R1Trim5: L1 = L1[R1Trim5:] 
                            if R1Trim3: L1 = L1[:-R1Trim3] 
                            if R2Trim5: L2 = L2[R2Trim5:] 
                            if R2Trim3: L1 = L2[:-R2Trim3]
                            len1 = len(L1)
                            len2 = len(L2)
                            if len1<klen1 and len2<klen1: continue
                            if MaxReadLen1:
                                L1 = L1[:MaxReadLen1]
                                L2 = L2[:MaxReadLen1]
                        else:
                            L1 = L1.strip()
                            TotalKmers += len(L1)-klen1+1
                            if NoNs1 and 'N' in L1: continue
                            if Trim5: L1 = L1[Trim5:] 
                            if Trim3: L1 = L1[:-Trim3]
                            len1 = len(L1)
                            if len1<klen1: continue
                            if MaxReadLen1:
                                L1 = L1[:MaxReadLen1]
                        if EndLinker1:
                            LinkPos1 = L1.rfind(EndLinker1)
                            if LinkPos1>=0:
                                L1 = L1[:LinkPos1]
                            else:
                                if RequireEndLinker1:
                                    continue
                            if MaxReadLen1: L1 = L1[:MaxReadLen1]
                            len1 = len(L1)
                            if len1<klen1: continue
                        PassFilterKmers += len1-klen1+1
                        PassFilterReads += 1
                        if TurboRareK:
                            KeepMe1 = False
                            p1 = klen1-TurboRareK
                            while not(KeepMe1) and p1+TurboRareK<=len1:
                                myHash = hash(L1[p1:p1+TurboRareK])&maskS0
                                myIndex = myHash & maskS1
                                myBit = myHash >> shiftS1
                                if TurboS1[myIndex]&(1<<myBit):
                                    KeepMe1 = True
                                else:
                                    p1 += klen1-TurboRareK+1
                                    if TurboRareLimit1 and p1+TurboRareK>TurboRareLimit1: break
                            p1 = klen1-TurboRareK
                            if Paired:
                                while not(KeepMe1) and p1+TurboRareK<=len2:
                                    myHash = hash(L2[p1:p1+TurboRareK])&maskS0
                                    myIndex = myHash & maskS1
                                    myBit = myHash >> shiftS1
                                    if TurboS1[myIndex]&(1<<myBit):
                                        KeepMe1 = True
                                    else:
                                        p1 += klen1-TurboRareK+1
                                        if TurboRareLimit1 and p1+TurboRareK>TurboRareLimit1: break
                            if not(KeepMe1):
                                continue
                        if NoDup1 and Paired: L1Set = set()
                        myRef = 0
                        BestRefCountForRead = 0
                        StartMatching1 = klen1-1
                        for j1,c1 in enumerate(L1):
                            if c1 == 'N':
                                StartMatching1 = j1+klen1
                            v0 = ((v0&mask1)<<2)+BaseL1[ord(c1)]
                            a0 = (a0>>2)+BaseA1[ord(c1)]
                            if j1>=StartMatching1:
                                if DataCadence1 and j1%DataCadence1: continue
                                vamin0 = min(a0,v0)
                                u = C1[vamin0]
                                if u==0: continue
                                if NoDup1 and Paired: L1Set.add(vamin0)
                                if not(Express1):
                                    myKmerC[vamin0] += 1
                                    myKmerL[-u] += 1
                                if u<-RefCount:
                                    if u<myRef:
                                        myRef =  u
                                        myk = vamin0
                                        if u<-2*RefCount and Express1:
                                            break
                                elif u>=-RefCount and myRef>=-RefCount and not(PrivateOnly1):
                                    myRefCount = C2[u]
                                    if not(BestRefCountForRead) or myRefCount<BestRefCountForRead:
                                        BestRefCountForRead = myRefCount
                                        myRef = u
                                        myk = vamin0
                        if Paired:
                            StartMatching1 = klen1-1
                            for j1,c1 in enumerate(L2):
                                if c1 == 'N':
                                    StartMatching1 = j1+klen1
                                v0 = ((v0&mask1)<<2)+BaseL1[ord(c1)]
                                a0 = (a0>>2)+BaseA1[ord(c1)]
                                if j1>=StartMatching1:
                                    if DataCadence1 and j1%DataCadence1: continue
                                    vamin0 = min(a0,v0)
                                    if NoDup1 and (vamin0 in L1Set): continue
                                    PassFilterKmers += 1                                        
                                    u = C1[vamin0]
                                    if u==0: continue
                                    if not(Express1):
                                        myKmerC[vamin0] += 1
                                        myKmerL[-u] += 1
                                    if u<-RefCount:
                                        if u<myRef:
                                            myRef =  u
                                            myk = vamin0
                                            if u<-2*RefCount and Express1:
                                                break
                                    elif u>=-RefCount and myRef>=-RefCount:
                                        myRefCount = C2[u]
                                        if not(BestRefCountForRead) or myRefCount<BestRefCountForRead:
                                            BestRefCountForRead = myRefCount
                                            myRef = u
                                            myk = vamin0
                        if myRef:
                            if 0>myRef>=-RefCount:
                                if Express1:
                                    myReadL[-myRef] += 1
                                else:
                                    myReadC[myk] += 1
                            else:
                                myReadL[-myRef] += 1
                            matchCount += 1
                vLog('Counting Completed '+myfn1+'Sample:'+sn+' Reads:'+str(i1//DataCycle1)+' Matched:'+str(matchCount))            
            except:
                vLog('Failure of kMerCounting for Data File:'+os.path.basename(F1)+' Sample:'+sn+' Reads:'+str(i1//DataCycle1)+' Matched :'+str(matchCount))
        if tools:
            return myKmerL,myReadL,TotalReads,TotalKmers,myKmerC,myReadC,PassFilterReads,PassFilterKmers
        else:
            vDeposit([myKmerL,myReadL,TotalReads,TotalKmers,myKmerC,myReadC,PassFilterReads,PassFilterKmers])
            myTask = vEat()
if vSpawned and vSpawnFunction=='CountWorker1':
    CountWorker1()
    quit()

print()
print('******************')
print('******************')
print('Very Preliminary Version of KMerCountRabbit-- this version for testing only (no guarantees of anything working or being correct/useful')
print('******************')
print('******************')
print()
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])
def dis1(x,d=3,m=0,c=9):
    '''diaplay a float x to allow d significant digits for distinctions from zero and 1, and mimimum m significant digits; max c characters'''
    if type(x)==str:
        return x
    if type(x)==int:
        return str(x)
    if type(x)==dict or type(x)==Counter:
        s = ''
        for k in x:
            s+='GATCN'[k]+':'+dis1(x[k])+', '
        return s[:-2]
    a = abs(x)
    if a>=10**(c-5):
        return f"{x:.3g}"
    if a==0 or a>=10**(d-m):
        sd = m
    elif a < 0.1:
        l = abs(int(-log(a,10)))
        sd = max(m,l+d)
    elif 0.9<=a<1.0:
        l = abs(int(-log(1-a,10)))
        sd = max(m,l+d)
    elif 0.1<=a<0.9:
        sd = max(m,d)
    elif 1.0<=a<10**(d-m):
        l = abs(int(log(a,10)))+1
        sd = max(m,d-l)
    else:
        sd = m
    if x<0:
        return '-'+(("%."+str(sd)+"f") % x)
    else:
        return ("%."+str(sd)+"f") % x
def ListToFileSpec1(myFiles):
    if type(myFiles)==str:
        myFiles = [myFiles,]
    mySpec = ''
    for m in myFiles:
        newSpec = ''
        for c in m:
            if c in '._': break
            if c=='*':
                newSpec+='x'
            if c.isalnum():
                newSpec+=c
        newSpec+='_'
        if not newSpec in mySpec:
            mySpec += newSpec        
    return mySpec[:100]
def multiglob1(x):
    if not(type(x))==str:
        x = ','.join(x)
    l = []
    for n in x.split(','):
        n = n.strip("'").strip('"')
        l.extend(sorted(list(glob(n))))
    return l
def findfastqdump(candidate):
    ver = 0
    try:
        ver = subprocess.check_output([candidate,'-V'])
        return candidate
    except:
        pass
    if os.path.isfile('~/FastQDumpLocation.txt'):
        candidate = open('~/FastQDumpLocation.txt', mode='rt').read()
        try:
            ver = subprocess.check_output([candidate,'-V'])
            return candidate
        except:
            pass
    vLog('Looking for fasterq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    vLog('or reinstall and allow execution (chmod +X <path>)')
    vLog('Note finding fast(er)q-dump can take some real time, get a cup of coffee or tea')
    NewCands = find1('fasterq-dump')
    NewItems = []
    for candidate in NewCands:
        try:
            ver = subprocess.check_output([candidate,'-V'])
            NewItems.append([versioner1(ver),candidate])
        except:
            pass
    if NewItems:
        candidate = sorted(NewItems, reverse=True)[0][1]
        try:
            open('~/FastQDumpLocation.txt', mode='w').write(candidate)
        except:
            pass
        return candidate
    vLog('Unable to find fast-q dump.  Recommend that you reinstall this from ncbi or download fastq/fasta files directly')
    return '' 
def FileListMnemonic1(FL):
    mn  = ''
    FLx = [os.path.basename(x) for x in FL]
    for i,c in enumerate(FLx[0]):
        for n in FLx:
            if len(n)<=i or n[i]!=c or n[i]=='.':
                if len(FL)>1:
                    return mn+'x'
                else:
                    return mn
        mn+=c
    return mn
def GuessReadNum1(f,l):
    '''guess read1 or read2 based on a filename and list of files; note that this is a "best guess"-- for definitive id, the command line should have
       a defined list of R1;R2 pairs separated by commas'''
    if type(l)==str:
        l = [l,]
    f = os.path.basename(f)
    l = set([os.path.basename(l0) for l0 in l])
    DefinitiveList = ('_R1.','_R1_','_R1-','-R1.','-R1_','-R1-','.R1.','.R1_','.R1-','-R1','.R1','_R1','R1.','R1_','R1-')
    PossibleList1 = ('R1',)
    for x1 in DefinitiveList+PossibleList1:
        x2 = x1.replace('1','2')
        if x1 in f and f.replace(x1,x2) in l:
            return 1,f.replace(x1,x2)
        if x2 in f and f.replace(x2,x1) in l:
            return 2,f.replace(x2,x1)
    for x1 in DefinitiveList+PossibleList1:
        x1 = x1.replace('R','r')
        x2 = x1.replace('1','2')
        if x1 in f and f.replace(x1,x2) in l:
            return 1,f.replace(x1,x2)
        if x2 in f and f.replace(x2,x1) in l:
            return 2,f.replace(x2,x1)
    for x1 in DefinitiveList:
        x1 = x1.replace('R','')
        x2 = x1.replace('1','2')
        if x1 in f and f.replace(x1,x2) in l:
            return 1,f.replace(x1,x2)
        if x2 in f and f.replace(x2,x1) in l:
            return 2,f.replace(x2,x1)
    if '_R1_' in f and '_1.' in f and f.replace('_R1_','_R2_').replace('_1.','_2.') in l:
        return 1,f.replace('_R1_','_R2_').replace('_1.','_2.')
    if '_R2_' in f and '_2.' in f and f.replace('_R2_','_R1_').replace('_2.','_1.') in l:
        return 2,f.replace('_R2_','_R1_').replace('_2.','_1.')
    return 0,''

## OutMenu1 can include 'Positions,FractionCovered,Counts,Depth,Incidence,KPKM'
ReferenceFileList1 = []
for s11 in ReferenceFiles1.split(','):
    if os.path.isfile(s11):
        ReferenceFileList1.append(s11)
    elif list(glob(s11)):
        ReferenceFileList1.extend(multiglob1(s11))
    elif ';' in s11:
        (sf1,sf2) = s11.split(';',1)
        if os.path.isfile(sf1):
            ReferenceFileList1.extend([sf1,sf2])
        else:
            ReferenceFileList1.extend(multiglob1(sf1)+multiglob1(sf2))
    if not(ReferenceFileList1):
        ReferenceFileList1.append(s11)
if not(ReferenceFileList1):
    vLog('No FastA dataset files found as reference [you must specify Reference=<a set of fasta or fastq files>]')
    quit()
ReferenceFileD1 = dict.fromkeys(ReferenceFileList1)
for Fn1 in list(ReferenceFileD1.keys()):
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)) and not('.' in Fn1) and Fn1[:3].isalpha() and Fn1[3:].isdigit():
        if not(os.path.isfile(Fn1+'_1.fastq')):
            vLog(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,Fn1])
            vLog("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        del(ReferenceFileD1[Fn1])
        if os.path.isfile(Fn1+'_1.fastq'):
            ReferenceFileD1[Fn1+'_1.fastq'] =  None
        else:
            vLog('Failed to download',Fn1,'- Will continue for now though')
        if os.path.isfile(Fn1+'_2.fastq'):
            ReferenceFileD1[Fn1+'_2.fastq'] = None
ReferenceFileList1 = list(ReferenceFileD1.keys())
vLog('Running ',os.path.basename(sys.argv[0]),'with parameters:','\r  '.join(sys.argv),'\r  #Python Version',sys.version)
C1 = Counter()             ## A large dictionary where keys are kmer numerical representation; values are owning reference sequence index, count (or zero if not encountered)
C2 = Counter()             ## Public repeat K mers, keys are kmer numerical representation; values are the numbers of hits in reference datasets 

SampleD1 = {}  ## Keys are sample names values are ordinal value
ReadNumC1 = Counter() ## Keys are file names, values are 1 for R1, 2 for R2, 0 for not known
DataFileList1 = []
if not(DataFiles1) and not(SampleToData1):
    vLog('No data files specified, will assume all .fasta/.fastq files (+/-.gz) in current directory:'+os.getcwd()+' ; use DATA=xxx in command line to specify files')
    DataFiles1  = './'
for s11 in multiglob1(DataFiles1):
    if not(s11): continue
    s11 = s11.strip().strip('"').strip("'")
    if os.path.isdir(s11):
        DataFileList1.extend(multiglob1(os.path.join(s11,'*.fasta'))+
                             multiglob1(os.path.join(s11,'*.fastq'))+
                             multiglob1(os.path.join(s11,'*.fasta.gz'))+
                             multiglob1(os.path.join(s11,'*.fastq.gz')))
    elif os.path.isfile(s11):
        DataFileList1.append(s11)
    elif list(glob(s11)):
        DataFileList1.extend(sorted(list(glob(s11))))
    elif ';' in s11:
        (sf1,sf2) = s11.split(';',1)
        sf1 = sf1.strip("' "+'"')
        sf2 = sf2.strip("' "+'"')
        if os.path.isfile(sf1) and os.path.isfile(sf2):
            DataFileList1.extend([sf1,sf2])
            ReadNumC1[sf1] = 1
            ReadNumC1[sf2] = 2
        else:
            DataFileList1.extend(sorted(list(glob(sf1)+list(glob(sf2)))))
            for sf11 in glob(sf1):
                ReadNumC1[sf11] = 1
            for sf22 in glob(sf2):
                ReadNumC1[sf22] = 2
    if not(DataFileList1):
        vLog('No data files were specified oroperly in the command line; check your files and indicate datafiles to process with Data=<filename(s)> in the command line')
        exit()
## A rather convoluted multi-pass set of routines to obtain an orderly list of input files
DataFileD1 = {}
for df1 in DataFileList1:
    DataFileD1[df1] = None
for Fn1 in list(DataFileD1.keys()):
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)) and not('.' in Fn1) and Fn1[:3].isalpha() and Fn1[3:].isdigit():
        if not(os.path.isfile(Fn1+'_1.fastq')):
            vLog(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,Fn1])
            vLog("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        del(DataFileD1[Fn1])
        if os.path.isfile(Fn1+'_1.fastq'):
            DataFileD1[Fn1+'_1.fastq'] = None
            ReadNumC1[Fn1+'_1.fastq'] = 1
        else:
            vLog('Failed to download',Fn1,'- Will continue for now though')
        if os.path.isfile(Fn1+'_2.fastq'):
            DataFileD1[Fn1+'_2.fastq'] = None
            ReadNumC1[Fn1+'_2.fastq'] = 2
SampleLabelColor1 = Counter()
SampleCategory1 = Counter()

if SampleToData1:
    s2d1 = myOpen1(SampleToData1,mode='rt').readlines()
    s2d11 = s2d1[0].replace(',',' ').replace(';','').replace('\t',' ').split()
    NoFiles1 = all([not(os.path.isfile(x) for x in s2d11)])
    if NoFiles1:
        HeaderExtension1 = s2d11.strip()
    for L1 in s2d1:
        mySample1 = L1.split('\t')[0].strip()
        for d1 in DataFileD1:
            if d1 in L1:
                DataFileD1[d1] = mySample1            
        for f1 in L1.replace(',',' ').replace('; ',';').replace('\t',' ').split():
            if ';' in f1:
                f1,f2 = f1.split(';',1)
                f1 = f1.strip(' "\t,'+"'")
                f2 = f2.strip(' "\t,'+"'")
                if (f1.lower().endswith('fasta.gz') or f1.lower().endswith('fa.gz') or f1.lower().endswith('fastq.gz') or \
                    f1.lower().endswith('fasta') or f1.lower().endswith('fa') or f1.lower().endswith('fastq')) and os.path.isfile(f1) and \
                    (f2.lower().endswith('fasta.gz') or f2.lower().endswith('fa.gz') or f2.lower().endswith('fastq.gz') or \
                    f2.lower().endswith('fasta') or f2.lower().endswith('fa') or f2.lower().endswith('fastq')) and os.path.isfile(f2):
                    DataFileD1[(f1,f2)] = mySample1
                    ReadNumC1[f1] = 1
                    ReadNumC1[f2] = 2
            else:
                f1 = f1.strip(' "\t,'+"'")
                if (f1.lower().endswith('fasta.gz') or f1.lower().endswith('fa.gz') or f1.lower().endswith('fastq.gz') or \
                    f1.lower().endswith('fasta') or f1.lower().endswith('fa') or f1.lower().endswith('fastq')) and os.path.isfile(f1):
                    DataFileD1[f1] = mySample1
        for f1 in L1.split('\t'):
            f1 = f1.strip()
            if f1.lower().startswith('rgb(') and f1.endswith(')'):
                SampleLabelColor1[mySample1] = f1.lower()
            elif f1.lower().startswith('color(') and f1.endswith(')'):
                SampleLabelColor1[mySample1] = f1.lower().split('(')[-1].split(')')[0]
            elif f1.lower().startswith('category(') and f1.endswith(')'):
                SampleCategory1[mySample1] = f1.lower().split('(')[-1].split(')')[0]
                    
DataFileList1 = list(DataFileD1.keys())
if AutoPair1:
    DataFileD2 = {}
    i1 = 0
    while i1<len(DataFileList1):
        if i1==len(DataFileList1)-1:
            if not(DataFileList1[i1]):
                DataFileD2[DataFileList1[i1]] = TruncName1(DataFileList1[i1])
            else:
                DataFileD2[DataFileList1[i1]] = DataFileD1[DataFileList1[i1]]
            break
        myPair1 = (DataFileList1[i1],DataFileList1[i1+1])
        myPairName1 = AreWeAPair1(*myPair1)
        if i1<len(DataFileList1)-1 and DataFileD1[DataFileList1[i1]]==DataFileD1[DataFileList1[i1+1]] and myPairName1:
            if DataFileD1[DataFileList1[i1]]:
                DataFileD2[myPair1] = DataFileD1[DataFileList1[i1]]
            else:
                DataFileD2[myPair1] = TruncName1(myPairName1)
            ReadNumC1[DataFileList1[i1]] = 1
            ReadNumC1[DataFileList1[i1+1]] = 1
            i1 += 2
        else:
            if DataFileD1[DataFileList1[i1]]:
                DataFileD2[DataFileList1[i1]] = DataFileD1[DataFileList1[i1]]
            else:
                DataFileD2[DataFileList1[i1]] = TruncName1(DataFileList1[i1])
            i1 += 1
    DataFileD1 = DataFileD2
DataFileList1 = list(DataFileD1.keys())  ##   A list where paired files are 2-ples
DataFileList2 = []                     ##   A flat list of each file
for d1 in DataFileList1:
    if type(d1)==str:
        DataFileList2.append(d1)
    else:
        DataFileList2.extend(d1)
    if not(DataFileD1[d1]):
        DataFileD1[d1] = TruncName1(d1)
    if type(d1)==str and not (d1 in ReadNumC1):
        ReadNumC1[d1] = GuessReadNum1(d1,DataFileList2)[0]
## Now generate every possible conversion between sample numbers, sample names, and sets of files.
## MetaSampleList1 is the key result with a list of files for each sample
Samples1 = sorted(list(set(DataFileD1.values())))
SampleNum1 = len(Samples1)
SampleToIndexC1 = {}
DataFileToIndexC1 = {}
for i1,Sam1 in enumerate(Samples1):
    SampleToIndexC1[Sam1] = i1
for df1 in DataFileD1:
    DataFileToIndexC1[df1] = SampleToIndexC1[DataFileD1[df1]]
MetaSampleList1 = [[] for i in range(len(Samples1))]  ## Keeps paired files as individual entries
MetaSampleList2 = [[] for i in range(len(Samples1))]  ## Each file is one entry
for df1 in DataFileD1:
    MetaSampleList1[DataFileToIndexC1[df1]].append(df1)
    if type[df1]==str:
        MetaSampleList2[DataFileToIndexC1[df1]].append(df1)
    else:
        MetaSampleList2[DataFileToIndexC1[df1]].extend(df1)
def myGetSize(Fn):
    if os.path.isfile(Fn):
        return os.path.getsize(Fn)
    else:
        return len(Fn)
TotalReferenceSize1 = sum(map(myGetSize,ReferenceFileList1))
vLog('Total Reference File Size =',TotalReferenceSize1)
TotalDataSize1 = sum(map(os.path.getsize,DataFileList2))
vLog('Total Data File Size =',TotalDataSize1)
if PrimaryIndex1=='default':
    if WriteIndex1:
        PrimaryIndex1 = 'Reference'
    elif TotalReferenceSize1<TotalDataSize1:
        PrimaryIndex1 = 'Reference'
        vLog("Primary Index mode set to 'Reference' based on file sizes. Rerunning with PrimaryIndex='Data' might help in some (rare) cases if you get a memory error")
    else:
        PrimaryIndex1 = 'Data'
        vLog("Primary Index mode set to 'Data' based on file sizes. Rerunning with PrimaryIndex='Reference' might help in some (rare) cases if you get a memory error")
if OutFileName1 == 'default':
    Chaser1 = '.tdv'
    pm1 = '-MatchCounts-'
    OutFileBase1 = FileListMnemonic1(DataFileList2)
    ReferenceMnemonic1 = FileListMnemonic1(ReferenceFileList1)
    OutFileName1 = FileListMnemonic1(DataFileList2)+pm1+FileListMnemonic1(ReferenceFileList1)+'_'+vnow+Chaser1            

FieldNameD0={
    'KmerCount'               :  7,   ## Raw Kmer Counts 
    'KmerCount2'              :  8,   ## Sum of Raw Kmer counts squared
    'ReadCount'               :  9,   ## Read Counts
    'KmersCovered'            : 10,   ## Number of kmers covered
    'KmerDepthAve'            : 11,   ## Average depth of coverage 
    'KmerDepthStdDev'         : 12,   ## Standard deviation of coverage depth
    'KmerCoverage'            : 13,   ## Fraction of Kmers covered in data
    'KmerSampleFraction'      : 14,   ## Fraction of pass-filter kmers in sample assigned to this partition
    'ReadSampleFraction'      : 15,   ## Fraction of pass-filter reads in sample assigned to this partition
    'KPKM'                    : 16,   ## Kmers mapped per thousand reference kmers per million pass-filter kmers
    'RPKM'                    : 17,   ## Reads mapped per thousand reference kmers per million pass-filter reads
    'KmerEnrichment'          : 18,   ## Enrichment of kmers mapped to this partition over observed for all datasets
    'ReadEnrichment'          : 19 } ## Enrichment of reads mapped to this partition over observed for all datasets
FieldNameD1={y:x for x,y in FieldNameD0.items()} 
for n1 in list(FieldNameD1.values()):
    FieldNameD1[n1] = n1
    FieldNameD1[str(n1)] = n1
for n0 in list(FieldNameD0.values()):
    FieldNameD0[n0] = n0
    FieldNameD0[str(n0)] = n0
for n1 in list(FieldNameD1.keys()):
    FieldNameD1[str(n1)] = FieldNameD1[n1]
for n1 in list(FieldNameD1.keys()):
    FieldNameD0[str(n1)] = FieldNameD0[n1]

if HeatMapTextFields1=='default':
    if Express1:
        HeatMapTextFields1 = 'ReadCount','RPKM'
    else:
        HeatMapTextFields1 = 'ReadCount','KmerDepthAve','KmerCoverage'
if TDVFields1=='default':
    if Express1:
        TDVFields1 = 'ReadCount'
    else:
        TDVFields1 = 'ReadCount','KmerCount','KmersCovered'
if not(type(TDVFields1) in (list,tuple)):
    TDVFields1 = [TDVFields1,]
if not(type(HeatMapTextFields1) in (list,tuple)):
    HeatMapTextFields1 = [HeatMapTextFields1,]
if HeatMapColorField1 == 'default':
    HeatMapColorField1 = 'RPKM'
HeatMapColorField0 = FieldNameD0[HeatMapColorField1] ## will be a field number
HeatMapColorField1 = FieldNameD1[HeatMapColorField1] ## will be a field name
TDVFields0 = [FieldNameD0[x] for x in TDVFields1]
TDVFields1 = [FieldNameD1[x] for x in TDVFields1]
HeatMapTextFields0 = [FieldNameD0[x] for x in HeatMapTextFields1]
HeatMapTextFields1 = [FieldNameD1[x] for x in HeatMapTextFields1]
AllFields1 = [HeatMapColorField0,]+TDVFields0+HeatMapTextFields0
KeepCounts2 = (8 in AllFields1) or (12 in AllFields1)  ## keeps track of whether we'll be calculating sum-of-squares for kmer counts

TSet1 = set()  ## For trimming reads on the fly
if TetritisF1:
    TetTrim1 = True
else:
    TetTrim1 = False
for FnT1 in multiglob1(TetritisF1):
    FD1 = vFastAToDict(FnT1,upper=True)
    v0 = 0
    a0 = 0
    for L1 in FD1.values():
        for j1,c1 in enumerate(L1.strip()):
            cd1 = BaseL1[ord(c1)]
            v0 = ((v0&Tmask1)<<2)+cd1
            a0 = (a0>>2)+(3-cd1)*Tsam1
            if j1>=klen1-1:
                TSet1.add(v0)
                TSet1.add(a0)
if (PrimaryIndex1=='Data'):
    for F1 in DataFileList2:
        DataCycle1 = 2 ## 2 for fasta, 4 for fastq
        if '.fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
        LastLine1 = DataCycle1*(MaxReadsPerFile1-1)
        if ReadNumC1[F1]==2:
            Trim5 = R2Trim5
            Trim3 = R2Trim3
        else:
            Trim5 = R1Trim5
            Trim3 = R1Trim3            
        v0 = 0
        a0 = 0
        ReportLineGranularity1 = DataCycle1*ReportGranularity1
        F1o = myOpen1(F1, mode='rt')
        for i1,L1 in enumerate(F1o):
            if i1%DataCycle1==1:
                if EndLinker1:
                    LinkPos1 = L1.rfind(EndLinker1)
                    if LinkPos1>=0:
                        L1 = L1[:LinkPos1]
                    else:
                        if RequireEndLinker1:
                            continue
                if MaxReadsPerFile1 and i1>LastLine1: break
                L1 = L1.strip()
                if Trim5:
                    L1 = L1[Trim5:] 
                if Trim3:
                    L1 = L1[:-Trim3]
                if MaxReadLen1: L1 = L1[:MaxReadLen1]
                if NoNs1 and 'N' in L1: continue
                if i1%ReportLineGranularity1==1:
                    vLog('Data File PreScreen:',os.path.basename(F1),' Read:',1+i1//DataCycle1, ' Kmers In Primary:', len(C1))
                StartMatching1 = klen1-1
                for j1,c1 in enumerate(L1):
                    if c1 == 'N':
                        StartMatching1 = j1+klen1
                    cd1 = BaseL1[ord(c1)]
                    v0 = ((v0&mask1)<<2)+cd1
                    a0 = (a0>>2)+(3-cd1)*ksam1
                    vamin0 = min(v0,a0)
                    if TetTrim1 and (Tmask0&v0 in TSet1):
                        StartMatching1 = j1+klen1
                    if j1>=StartMatching1:
                        C1[vamin0] = 0
        myClose1(F1o)
        vLog('Finishing Data File PreScreen:',os.path.basename(F1),' Reads:',i1//DataCycle1, ' Kmers In Primary:', len(C1))
SkipDataPreIndex1 = False
if PrimaryIndex1.lower() != 'data':
    SkipDataPreIndex1 = True

def RefLenCovFilter1(n,s):
    '''input is name, sequence from a fastA assembly from spades or megahit; output is True (Filter this reference out) or False (use it for analysis)''' 
    if MinRefCov1:
        if 'multi=' in n:
            myMulti1 = float(n.split('multi=')[1].split()[0])
            if myMulti1<MinRefCov1: return True
        if '_cov_' in n:
            myLen1 = float(n.split('_cov_')[1].split('_')[0])
            if myLen1<MinRefCov1: return True
    if MinRefLen1:
        if 'len=' in n:
            myLen1 = float(n.split('len=')[1].split()[0])
            if myLen1<MinRefLen1: return True
        if '_length_' in n:
            myLen1 = float(n.split('_length_')[1].split('_')[0])
            if myLen1<MinRefLen1: return True
        if len(s)<MinRefLen1: return True

## Generate index from the input reference files
if PrimaryIndex1.lower() != 'reference' and PrimaryIndex1.lower() != 'data':
    pf = open(PrimaryIndex1, mode='rb')
    pp = pickle.Unpickler(pf)
    C1,C2,D11,RefName1,RefFile1,RefLen1,RefCount1 = pp.load()
    pf.close()
else:
    D11 = {}
    RefLen1 = [0]
    RefName1 = ['NoMatch']
    RefFile1 = ['None']
    Started1 = False
    ## We'll need a count of reference sequences
    if ReferenceFileByFile1:
        RefCount1 = len(ReferenceFileList1)
    else:
        RefCount1 = 0
        for F1 in ReferenceFileList1:
            for L1 in chain(vOpen(F1,mode='rt'),['>',]):
                if L1.startswith('>'):
                    if Started1:
                        if MinRefCov1:
                            if 'multi=' in L1:
                                myMulti1 = float(L1.split('multi=')[1].split()[0])
                                if myMulti1<MinRefCov1: continue
                            if '_cov_' in L1:
                                myLen1 = float(L1.split('_cov_')[1].split('_')[0])
                                if myLen1<MinRefCov1: continue
                        if MinRefLen1:
                            if 'len=' in L1:
                                myLen1 = float(L1.split('len=')[1].split()[0])
                                if myLen1<MinRefLen1: continue
                            if '_length_' in L1:
                                myLen1 = float(L1.split('_length_')[1].split('_')[0])
                                if myLen1<MinRefLen1: continue
                            if ls1<MinRefLen1: continue
                        RefCount1 += 1
                    ls1 = 0
                    Started1 = True
                else:
                    ls1 += len(L1.strip())
    for refrun1 in (0,1):
        if refrun1 and not(ReferenceCadence1): continue
        refnum1 = 0
        for F1 in ReferenceFileList1:
            vLog('Starting Reference File PreScreen: (run '+str(refrun1)+')',os.path.basename(F1))
            if ReferenceFileByFile1:
                D1 = vFastAToDict(F1,OneReference=os.path.basename(F1),upper=True)
            else:
                D1 = vFastAToDict(F1,upper=True)
            for n11,s11 in D1.items():
                if MinRefCov1:
                    if 'multi=' in n11:
                        myMulti1 = float(n11.split('multi=')[1].split()[0])
                        if myMulti1<MinRefCov1: continue
                    if '_cov_' in n11:
                        myLen1 = float(n11.split('_cov_')[1].split('_')[0])
                        if myLen1<MinRefCov1: continue
                if MinRefLen1:
                    if 'len=' in n11:
                        myLen1 = float(n11.split('len=')[1].split()[0])
                        if myLen1<MinRefLen1: continue
                    if '_length_' in n11:
                        myLen1 = float(n11.split('_length_')[1].split('_')[0])
                        if myLen1<MinRefLen1: continue
                    if len(s11)<MinRefLen1: continue
                if CircularReference1: s11 += s11[:klen1-1]
                if refrun1 == 0: 
                    RefLen1.append(len(s11))
                    RefName1.append(n11.replace('\t',' '))
                    RefFile1.append(os.path.basename(F1))
                    D11[(n11,F1)] = s11
                refnum1 += 1        
                v1 = 0; a1 = 0
                PublicSet1 = set()
                StartMatching1 = klen1-1
                for i1,c1 in enumerate(s11):
                    if c1 == 'N':
                        StartMatching1 = i1+klen1
                    cd1 = BaseL1[ord(c1)]
                    v1 = ((v1&mask1)<<2)+cd1
                    a1 = (a1>>2)+(3-cd1)*ksam1
                    if TetTrim1 and (Tmask0&v1 in TSet1):
                        StartMatching1 = i1+klen1
                    if i1>=StartMatching1:
                        if ReferenceCadence1 and (not(refrun1)!=(not(i1%ReferenceCadence1))): continue
                        vamin1 = min(a1,v1)
                        if SkipDataPreIndex1 or (vamin1 in C1):
                            u1 = C1[vamin1]
                            if not(u1) and not(refrun1):
                                C1[vamin1] = -refnum1-RefCount1*2
                            elif u1<-2*RefCount1:
                                if u1==-refnum1-2*RefCount1:
                                    C1[vamin1] = u1 + RefCount1
                                else:
                                    C1[vamin1] = max(-refnum1,u1 + 2*RefCount1)
                            elif u1<-RefCount1 and u1!=(-refnum1-RefCount1):
                                C1[vamin1] = max(-refnum1,u1 + RefCount1)
                            elif not(PrivateOnly1) and (not(refrun1) or not(ReferenceCadence1)):
                                if not(vamin1 in PublicSet1):
                                    C2[vamin1] +=1
                                    PublicSet1.add(vamin1)
        vLog('Finishing Reference File PreScreen:',os.path.basename(F1))
    if PrivateOnly1:
        C1temp1 = Counter()
        for x in C1:
            if C1[x]<-RefCount1:
                C1temp1[x] = C1[x]
        C1 = C1temp1    
    if FilterList1:
        C1temp1 = Counter()
        C1temp2 = Counter()
        FilterList1 = [x.replace('"','').replace(' ','').replace("'",'') for x in FilterList1.split(';')]    
        for g1 in FilterList1:
            if '>' in g1:
                if g1.split('>',1)[0].isdigit:
                    gn1 = int(g1.split('>',1)[0])
                    gf1 = g1.split('>',1)[1]
                    gc1 = '>'
                else:
                    g1 = g1[::-1]
                    gn1 = int(g1.split('>',1)[0][::-1])
                    gf1 = g1.split('>',1)[1][::-1]
                    gc1 = '<'                
            elif '<' in g1:
                if g1.split('<',1)[0].isdigit:
                    gn1 = int(g1.split('<',1)[0])
                    gf1 = g1.split('<',1)[1]
                    gc1 = '<'
                else:
                    g1 = g1[::-1]
                    gn1 = int(g1.split('<',1)[0][::-1])
                    gf1 = g1.split('<',1)[1][::-1]
                    gc1 = '>'                
            if gf1.lower()=='reference':
                myseqD1 = [D11,]
            else:
                myseqD1 = map(vFastAToDict1,multiglob1(gf1))
            for Dx1 in myseqD1:
                for s11 in Dx1.values():
                    StartMatching1 = klen1-1
                    for i1,c1 in enumerate(s11):
                        if c1 == 'N':
                            StartMatching1 = j1+klen1
                        cd1 = BaseL1[ord(c1)]
                        v1 = ((v1&mask1)<<2)+cd1
                        a1 = (a1>>2)+(3-cd1)*ksam1
                        if i1>=StartMatching1:
                            vamin1 = min(a1,v1)
                            if vamin1 in C1:
                                C1temp1[vamin1] += 1
            if gc1=='<':
                for x1 in C1:
                    if gn1<C1temp1[x1]:
                        C1temp2[x1] = C1[x1]
            elif gc1=='>':
                for x1 in C1:
                    if gn1>C1temp1[x1]:
                        C1temp2[x1] = C1[x1]
            C1 = C1temp2
    if DustFilter1:
        vLog('Starting Dust Filtering')
        C1 = Counter({k:C1[k] for k in C1 if not(SimpleDust1(k))})
        vLog('Finished Dust Filtering')       
    FrequencyFilter1 = ColdKMerF1 or HotKMerF1 or ColdKRelToMedian1 or HotKRelToMedian1 
    if FrequencyFilter1:
        C0 = Counter()
        for F1 in DataFileList2:
            DataCycle1 = 2 ## 2 for fasta, 4 for fastq
            if '.fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
            LastLine1 = DataCycle1*(MaxReadsPerFile1-1)
            if ReadNumC1[F1]==2:
                Trim5 = R2Trim5
                Trim3 = R2Trim3
            else:
                Trim5 = R1Trim5
                Trim3 = R1Trim3            
            v0 = 0
            a0 = 0
            ReportLineGranularity1 = DataCycle1*ReportGranularity1
            F1o = myOpen1(F1, mode='rt')
            for i1,L1 in enumerate(F1o):
                if i1%DataCycle1==1:
                    if EndLinker1:
                        LinkPos1 = L1.rfind(EndLinker1)
                        if LinkPos1>=0:
                            L1 = L1[:LinkPos1]
                        else:
                            if RequireEndLinker1:
                                continue
                    if MaxReadsPerFile1 and i1>LastLine1: break
                    L1 = L1.strip()
                    if Trim5:
                        L1 = L1[Trim5:] 
                    if Trim3:
                        L1 = L1[:-Trim3]
                    if MaxReadLen1: L1 = L1[:MaxReadLen1]
                    if NoNs1 and 'N' in L1: continue
                    if i1%ReportLineGranularity1==1:
                        vLog('KMer Rarefication Run:',os.path.basename(F1),' Read:',1+i1//DataCycle1, ' Kmers In Primary:', len(C1))
                    StartMatching1 = klen1-1
                    for j1,c1 in enumerate(L1):
                        if c1 == 'N':
                            StartMatching1 = j1+klen1
                        cd1 = BaseL1[ord(c1)]
                        v0 = ((v0&mask1)<<2)+cd1
                        a0 = (a0>>2)+(3-cd1)*ksam1
                        vamin0 = min(v0,a0)
                        if TetTrim1 and (Tmask0&v0 in TSet1):
                            StartMatching1 = j1+klen1
                        if j1>=StartMatching1:
                            if vamin0 in C1:
                                C0[vamin0] += 1
            myClose1(F1o)
        xD0 = {x1:[] for x1 in C1.values()} ## Key is the C1[k] value, value is a list of counts
        xD2 = Counter() ## MaxCount Values allowed
        xD1 = Counter() ## MinCount Values allowed
        for k1 in C0:
            xD0[C1[k1]].append(C0[k1])
        for x1 in xD0:
            xD0[x1].sort()
            lx1 = len(xD0[x1])
            if lx1%2==0:
                median1 = (xD0[x1][lx1//2-1]+xD0[x1][lx1//2])/2
            else:
                median1 = xD0[x1][lx1//2]
            xD1[x1] = 0
            xI1 = int(len(xD0[x1])* ColdKMerF1)
            if xI1>0: xD1[x1] = max(xD1[x1],xD0[x1][xI1-1])
            if ColdKRelToMedian1: xD1[x1] = max(xD1[x1],ColdKRelToMedian1*median1) ## If both median and rank-based filtering are chosen, keep only Kmers that satisfy both filters
            xD2[x1] = 2**63-1
            xI2 = int(len(xD0[x1])* HotKMerF1)-1
            if xI2>0:xD2[x1] = xD0[x1][xI2]
            if HotKRelToMedian1: xD2[x1] = min(xD2[x1],HotKRelToMedian1*median1) ## If both median and rank-based filtering are chosen, keep only Kmers that satisfy both filters
        C1 = Counter({k1:C1[k1] for k1 in C1 if xD1[C1[k1]]<=C0[k1]<=xD2[C1[k1]]})
        del(xD0); del(C0)
    
if WriteIndex1:
    if WriteIndex1 == True or WriteIndex1.lower()=='reference':
        WriteIndex1 = 'KCRIndex_'+vnow+'.pck'
    pf = open(WriteIndex1, mode='wb')
    pp = pickle.Pickler(pf)
    pp.dump([C1,C2,D11,RefName1,RefFile1,RefLen1,RefCount1])
    pf.close()                
            
if TurboRareK1 == 'Auto':
    if DataCadence1 or ReferenceCadence1:
        TurboRareK1 = 0
        ConservativeInstanceLimit1 = max(1,os.cpu_count()-1)
    else:
        SampleReadLen1 = EstimatedAverageReadLen1
        if TurboRareLimit1:
            SampleReadLen1 = min(TurboRareLimit1,EstimatedAverageReadLen1)
        TurboRareK1 = ceil(log(len(C1)*20*SampleReadLen1/(klen1-15),2)/2) ## ~~ 95% of non-matching reads should be triaged
        if TurboRareK1 < 14:
            TurboRareK1 = 14
        try:
            MyMem1 = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') 
        except:
            MyMem1 = 2**34 ## Getting total memory for a windows system here can't work with the same code as used for linux-mac-- thus just assuming 16G for now..  Could provide a better routine for this
        MemPerInstance1 = (4**TurboRareK1)//8+60*len(C1)+10**7 ## Extremely rough estimate of memory per instance
        if MemPerInstance1 > MaxMemPerThread1*MyMem1:
            TurboRareK1 -= 1
            MemPerInstance1 = (4**TurboRareK1)//8+60*len(C1)+10**7 
        if MemPerInstance1 > MaxMemPerThread1*MyMem1:
            TurboRareK1 -= 1
            MemPerInstance1 = (4**TurboRareK1)//8+60*len(C1)+10**7 
        if MemPerInstance1 > MaxMemPerThread1*MyMem1:
            TurboRareK1 -= 0
            MemPerInstance1 = 60*len(C1)+10**7
        if MemPerInstance1>MyMem1:
            vLog('************** Projecting Insufficient Memory; Will try with single thread, but stop program if system starts to hang')
            Threads1 = 0
            TurboRareK1 -= 0
        vLog('TurboRareK setting: '+str(TurboRareK1))
        ConservativeInstanceLimit1 = max(1,(MyMem1-2**30)//MemPerInstance1)
else:
    ConservativeInstanceLimit1 = max(1,os.cpu_count()-1)
if Threads1 == 'Auto' or Threads1==True:
    Threads1 = min(ConservativeInstanceLimit1,os.cpu_count(),SampleNum1)
    if Threads1 == 1:
        Threads1 = 0
    vLog('Thread Number Setting: '+str(Threads1))
elif not(Threads1):
    Threads1 = 0
if TurboRareK1 and Threads1==0:
    TurboS1 = array('B',[0]*(2**(2*TurboRareK1-3)))
    maskS0 = 2**(2*TurboRareK1)-1
    maskS1 = 2**(2*TurboRareK1-3)-1
    shiftS1 = 2*TurboRareK1-3
    for s11 in D11.values():
        a11 = vantisense(s11)
        for p1 in range(len(s11)-TurboRareK1+1):
            myHash = hash(s11[p1:p1+TurboRareK1])&maskS0
            myIndex = myHash & maskS1
            myBit = myHash >> shiftS1
            TurboS1[myIndex] |= 1<<myBit
            myHash = hash(a11[p1:p1+TurboRareK1])&maskS0
            myIndex = myHash & maskS1
            myBit = myHash >> shiftS1
            TurboS1[myIndex] |= 1<<myBit
else:
    TurboS1 = Counter()
## Do global calculations of different classes of Kmers (unique, local repeat, global) in the reference sequences
if FirstDibs1 or PrivateOnly1:
    D11 = Counter({x:'' for x in D11})
RefUniqueKmerSpecies1 = Counter()
RefUniqueKmerPositions1 = Counter()
RefPublicKmerPositions1 =  Counter()
RefPublicKmerSpecies1 = Counter()
if AggregateUniqueAndLocalRepeats1:
    RefLocalRepeatKmerSpecies1 = RefUniqueKmerSpecies1
    RefLocalRepeatKmerPositions1 =  RefUniqueKmerPositions1
else:
    RefLocalRepeatKmerSpecies1 = Counter()
    RefLocalRepeatKmerPositions1 =  Counter()

    
refnum1 = 0
C1Counts1 = Counter(C1.values())

for (n11,f11),s11 in D11.items():
    PublicList1 = []
    refnum1 += 1
    RefUniqueKmerSpecies1[refnum1] += C1Counts1[-refnum1-2*RefCount1]
    RefLocalRepeatKmerSpecies1[refnum1] += C1Counts1[-refnum1-RefCount1]
    if not(AggregateUniqueAndLocalRepeats1) or not(PrivateOnly1):
        v1 = 0; a1 = 0
        StartMatching1 = klen1-1
        for i1,c1 in enumerate(s11):
            if c1 == 'N':
                StartMatching1 = i1+klen1
            v1 = ((v1&mask1)<<2)+BaseL1[ord(c1)]
            a1 = (a1>>2)+BaseA1[ord(c1)]
            if (i1>=StartMatching1):
                vamin1 = min(a1,v1)
                u1 = C1[vamin1]
                if 0>u1>=-RefCount1:
                    PublicList1.append(vamin1)
        RefPublicKmerSpecies1[refnum1] = len(set(PublicList1))
        RefPublicKmerPositions1[refnum1] = len(PublicList1)
        RefLocalRepeatKmerPositions1[refnum1] = max(0,len(s11)-klen1+1-RefUniqueKmerSpecies1[refnum1]- RefPublicKmerPositions1[refnum1])

def SamplesSum1(JCounter):
    result = []
    for i in range(SampleNum1):
        result.append(sum([JCounter[(r,i)] for r in range(1,RefCount1+1)]))
    return result

JUniqueReads1 = Counter()
JUniqueKmerCounts1 = Counter()
JUniqueKmerCounts2 = Counter()
JUniqueKmerCovered1 = Counter()
JPublicReads1 = Counter()
JPublicKmerCounts1 = Counter()
JPublicKmerCounts2 = Counter()
JPublicKmerCovered1 = Counter()
if AggregateUniqueAndLocalRepeats1:
    JLocalRepeatReads1 = JUniqueReads1
    JLocalRepeatKmerCounts1 = JUniqueKmerCounts1
    JLocalRepeatKmerCounts2 = JUniqueKmerCounts2
    JLocalRepeatKmerCovered1 = JUniqueKmerCovered1
else:
    JLocalRepeatReads1 = Counter()
    JLocalRepeatKmerCounts1 = Counter()
    JLocalRepeatKmerCounts2 = Counter()
    JLocalRepeatKmerCovered1 = Counter()

RefMatches1 = [0] * (RefCount1+1) ## Zero if no matches, positive value if there are matches
SamplePassFilterKmers1 =  [0] * SampleNum1
SamplePassFilterReads1 =  [0] * SampleNum1
SampleTotalReads1 =  [0] * SampleNum1
SampleTotalKmers1 =  [0] * SampleNum1

if Threads1 == 0:
    for f1,fL1 in enumerate(MetaSampleList1):
        myKmerL1,myReadL1,TotalReads1,TotalKmers1,myKmerC1,myReadC1,PassFilterReads1,PassFilterKmers1 = CountWorker1(f1,fL1=fL1,sn=Samples1[f1],tools=[C1,C2,ReferenceFileList1,RefCount1,ReadNumC1,TetTrim1,TSet1,TurboRareK1,Threads1],TurboS1=TurboS1)
## Start of "duplicated block"
        ## Here is the indent for the block
        SamplePassFilterKmers1[f1] = PassFilterKmers1
        SamplePassFilterReads1[f1] = PassFilterReads1
        SampleTotalReads1[f1] = TotalReads1
        SampleTotalKmers1[f1] = TotalKmers1
        if not(myKmerC1) and not(myReadC1) and sum(myKmerL1)==0 and sum(myReadL1)==0: continue
        if not(Express1):
            CoverageByRef1 = Counter([-C1[x] for x in myKmerC1])
            FoundPublicK1 = sum([CoverageByRef1[i] for i in range(1,RefCount1+1)])
        refnum1 = 0
        if KeepCounts2:
            myKmerCounts2 = Counter()
            for k1,n1 in myKmerC1.items():
                myKmerCounts2[C1[k1]] += n1**2
        for (n11,f11) in D11:
            refnum1 += 1
            JUniqueReads1[(refnum1,f1)] += myReadL1[refnum1+2*RefCount1]
            JLocalRepeatReads1[(refnum1,f1)] += myReadL1[refnum1+RefCount1]
            JPublicReads1[(refnum1,f1)] += myReadL1[refnum1]
            if not(Express1):
                JUniqueKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1+2*RefCount1]
                JLocalRepeatKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1+RefCount1]
                JUniqueKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1+2*RefCount1]
                JLocalRepeatKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1+RefCount1]
                if KeepCounts2:
                    JUniqueKmerCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1-2*RefCount1]
                    JLocalRepeatKmerCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1-RefCount1]
                if FoundPublicK1 and not(PrivateOnly1): 
                    if FirstDibs1:
                        JPublicKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1]
                        JPublicKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1]
                        if KeepCounts2:
                            JPublicCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1]
            if not(FirstDibs1) and FoundPublicK1 and RefPublicKmerPositions1[refnum1] and not(PrivateOnly1):
                s11 = D11[(n11,f11)]
                v1 = 0; a1 = 0
                StartMatching1 = klen1-1
                PublicSet1 = set()
                for i1,c1 in enumerate(s11):
                    if c1 == 'N':
                        StartMatching1 = i1+klen1
                    v1 = ((v1&mask1)<<2)+BaseL1[ord(c1)]
                    a1 = (a1>>2)+BaseA1[ord(c1)]
                    if (i1>=StartMatching1):
                        vamin1 = min(a1,v1)
                        u1 = C1[vamin1]
                        if 0>u1>=-RefCount1:
                            JPublicReads1[(refnum1,f1)] += myReadC1[vamin1]
                            if not(Express1):
                                myc1 = myKmerC1[vamin1]
                                if not(myc1): continue
                                if not(vamin1 in PublicSet1):
                                    PublicSet1.add(vamin1)
                                    JPublicKmerCounts1[(refnum1,f1)] +=  myc1
                                    JPublicKmerCovered1[(refnum1,f1)] += 1
                                    if KeepCounts2:
                                        JPublicKmerCounts2[(refnum1,f1)] +=  myc1**2
            if Express1:
                RefMatches1[refnum1] += JUniqueReads1[(refnum1,f1)]
                if not(AggregateUniqueAndLocalRepeats1):
                    RefMatches1[refnum1] += JLocalRepeatReads1[(refnum1,f1)]
                if not(PrivateOnly1):
                    RefMatches1[refnum1] += JPublicReads1[(refnum1,f1)]
            else:
                RefMatches1[refnum1] += JUniqueReads1[(refnum1,f1)]
                if not(AggregateUniqueAndLocalRepeats1):
                    RefMatches1[refnum1] += JLocalRepeatReads1[(refnum1,f1)]
                if not(PrivateOnly1):
                    RefMatches1[refnum1] += JPublicReads1[(refnum1,f1)]

 ## End of "duplicated block"
else:
    TaskList1 = []
    vShare([C1,C2,ReferenceFileList1,RefCount1,ReadNumC1,TetTrim1,TSet1,TurboRareK1,Threads1],"BasicData")
    for pt1,ptL1 in enumerate(MetaSampleList1):
        TaskList1.append([pt1,ptL1,Samples1[pt1]])
    CoreNum1 = min(Threads1,len(TaskList1),max(1,os.cpu_count()-1))
    ProcessList1 = []
    for i1 in range(CoreNum1):
        ProcessList1.append(vProcess('CountWorker1',[],i1))
        ProcessList1[-1].spawn(CoreNum1)
    ProcessStatus1 = [0]*CoreNum1  ## 0 for not run, 1 for running, 2 for complete
    ProcessTask1 = [-1]*CoreNum1
    nextTask1 = 0
    while min(ProcessStatus1)<2:
        sleep(0.02)
        for i11 in range(CoreNum1):
            if ProcessStatus1[i11]==0:
                ProcessList1[i11].feed(TaskList1[nextTask1])
                ProcessTask1[i11] = nextTask1
                nextTask1 += 1
                ProcessStatus1[i11] = 1
            elif ProcessStatus1[i11]==1 and ProcessList1[i11].test():
                myTask1 = ProcessTask1[i11]
                myKmerL1,myReadL1,TotalReads1,TotalKmers1,myKmerC1,myReadC1,PassFilterReads1,PassFilterKmers1 = ProcessList1[i11].harvest()
                if nextTask1<len(TaskList1):
                    ProcessList1[i11].feed(TaskList1[nextTask1])
                    ProcessTask1[i11] = nextTask1
                    nextTask1 += 1
                else:
                    ProcessList1[i11].cleanup()
                    ProcessList1[i11].feed('')
                    ProcessStatus1[i11] = 2
                f1,fL1,sn1 = TaskList1[myTask1]
## Start of "duplicated block"
                ## Here is the indent for the block
                SamplePassFilterKmers1[f1] = PassFilterKmers1
                SamplePassFilterKmers1[f1] = PassFilterKmers1
                SamplePassFilterReads1[f1] = PassFilterReads1
                SampleTotalReads1[f1] = TotalReads1
                SampleTotalKmers1[f1] = TotalKmers1
                if not(myKmerC1) and not(myReadC1) and sum(myKmerL1)==0 and sum(myReadL1)==0: continue
                if not(Express1):
                    CoverageByRef1 = Counter([-C1[x] for x in myKmerC1])
                    FoundPublicK1 = sum([CoverageByRef1[i] for i in range(1,RefCount1+1)])
                refnum1 = 0
                if KeepCounts2:
                    myKmerCounts2 = Counter()
                    for k1,n1 in myKmerC1.items():
                        myKmerCounts2[C1[k1]] += n1**2
                for (n11,f11) in D11:
                    refnum1 += 1
                    JUniqueReads1[(refnum1,f1)] += myReadL1[refnum1+2*RefCount1]
                    JLocalRepeatReads1[(refnum1,f1)] += myReadL1[refnum1+RefCount1]
                    JPublicReads1[(refnum1,f1)] += myReadL1[refnum1]
                    if not(Express1):
                        JUniqueKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1+2*RefCount1]
                        JLocalRepeatKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1+RefCount1]
                        JUniqueKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1+2*RefCount1]
                        JLocalRepeatKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1+RefCount1]
                        if KeepCounts2:
                            JUniqueKmerCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1-2*RefCount1]
                            JLocalRepeatKmerCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1-RefCount1]
                        if FoundPublicK1 and not(PrivateOnly1): 
                            if FirstDibs1:
                                JPublicKmerCounts1[(refnum1,f1)] += myKmerL1[refnum1]
                                JPublicKmerCovered1[(refnum1,f1)] += CoverageByRef1[refnum1]
                                if KeepCounts2:
                                    JPublicCounts2[(refnum1,f1)] += myKmerCounts2[-refnum1]
                    if not(FirstDibs1) and FoundPublicK1 and RefPublicKmerPositions1[refnum1] and not(PrivateOnly1):
                        s11 = D11[(n11,f11)]
                        v1 = 0; a1 = 0
                        StartMatching1 = klen1-1
                        PublicSet1 = set()
                        for i1,c1 in enumerate(s11):
                            if c1 == 'N':
                                StartMatching1 = j1+klen1
                            v1 = ((v1&mask1)<<2)+BaseL1[ord(c1)]
                            a1 = (a1>>2)+BaseA1[ord(c1)]
                            if (i1>=StartMatching1):
                                vamin1 = min(a1,v1)
                                u1 = C1[vamin1]
                                if 0>u1>=-RefCount1:
                                    JPublicReads1[(refnum1,f1)] += myReadC1[vamin1]
                                    if not(Express1):
                                        myc1 = myKmerC1[vamin1]
                                        if not(myc1): continue
                                        if not(vamin1 in PublicSet1):
                                            PublicSet1.add(vamin1)
                                            JPublicKmerCounts1[(refnum1,f1)] +=  myc1
                                            JPublicKmerCovered1[(refnum1,f1)] += 1
                                            if KeepCounts2:
                                                JPublicKmerCounts2[(refnum1,f1)] +=  myc1**2
                    if Express1:
                        RefMatches1[refnum1] += JUniqueReads1[(refnum1,f1)]
                        if not(AggregateUniqueAndLocalRepeats1):
                            RefMatches1[refnum1] += JLocalRepeatReads1[(refnum1,f1)]
                        if not(PrivateOnly1):
                            RefMatches1[refnum1] += JPublicReads1[(refnum1,f1)]
                    else:
                        RefMatches1[refnum1] += JUniqueReads1[(refnum1,f1)]
                        if not(AggregateUniqueAndLocalRepeats1):
                            RefMatches1[refnum1] += JLocalRepeatReads1[(refnum1,f1)]
                        if not(PrivateOnly1):
                            RefMatches1[refnum1] += JPublicReads1[(refnum1,f1)]


    vCleanup("BasicData")

## Calculate global match numbers for each reference and Sample
RefPublicReadsFound1 = [sum([JPublicReads1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
RefLocalRepeatReadsFound1 = [sum([JLocalRepeatReads1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
RefUniqueReadsFound1 = [sum([JUniqueReads1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
RefPublicKmersFound1 = [sum([JPublicKmerCounts1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
RefLocalRepeatKmersFound1 = [sum([JLocalRepeatKmerCounts1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
RefUniqueKmersFound1 = [sum([JUniqueKmerCounts1[(r1,f1)] for f1 in range(SampleNum1)]) for r1 in range(RefCount1+1)]
AllRefUniqueKmer1 = sum(RefUniqueKmerSpecies1.values())-RefUniqueKmerSpecies1[0]
AllRefLocalRepeatKmer1 = sum(RefLocalRepeatKmerPositions1.values())-RefLocalRepeatKmerPositions1[0]
AllRefPublicKmer1 = sum(RefPublicKmerPositions1.values())-RefPublicKmerPositions1[0]
AllRefLen1 = sum(RefLen1)-RefLen1[0]
AllRefKmer1 = sum([max(0,x-klen1+1) for x in RefLen1])

SampleUniqueReadMatches1 = SamplesSum1(JUniqueReads1)
SampleLocalRepeatReadMatches1 = SamplesSum1(JLocalRepeatReads1)
SamplePublicReadMatches1 = SamplesSum1(JPublicReads1)
SampleUniqueKmerMatches1 = SamplesSum1(JUniqueKmerCounts1)
SampleLocalRepeatKmerMatches1 = SamplesSum1(JLocalRepeatKmerCounts1)
SamplePublicKmerMatches1 = SamplesSum1(JPublicKmerCounts1)
SampleUniqueKmerCovered1 = SamplesSum1(JUniqueKmerCovered1)
SampleLocalRepeatKmerCovered1 = SamplesSum1(JLocalRepeatKmerCovered1)
SamplePublicKmerCovered1 = SamplesSum1(JPublicKmerCovered1)
SampleMatches1 = []
for f1 in range(SampleNum1):
    if Express1:
        if SampleUniqueReadMatches1[f1] or  SampleLocalRepeatReadMatches1[f1] or  SamplePublicReadMatches1[f1]:
            SampleMatches1.append(True)
        else:
            SampleMatches1.append(False)
    else:        
        if SampleUniqueKmerMatches1[f1] or  SampleLocalRepeatKmerMatches1[f1] or  SamplePublicKmerMatches1[f1]:
            SampleMatches1.append(True)
        else:
            SampleMatches1.append(False)
def zad1(a,b):
    '''zero-aware-division'''
    if b==0:
        return 0.0
    return a/b
E9 = 10**9

OutFile1 = myOpen1(OutFileName1, mode='wt')
TaskHeader1 =  '<!--Rabbit_Task_Header: '+OutFileName1+'-->'+Delimiter1
TaskHeader1 += '<!--Command_Line: '+' '.join(sys.argv)+'-->'+Delimiter1
TaskHeader1 += '<!--PythonVersion: '+','.join(sys.version.splitlines())+'-->'+Delimiter1
TaskHeader1 += '<!--Rabbit_Version: '+FileInfo1(sys.argv[0])+'-->'+Delimiter1
TaskHeader1 += '<!--RunTime: '+vnow+'-->'+Delimiter1
TaskHeader1 += '<!--RunDirectory: '+os.getcwd()+'-->'+Delimiter1
TaskHeader1+=Delimiter1
AbbrevHeader1 = ''.join(TaskHeader1.splitlines()[:-1])+'<!--RabbitTableHeader-->'  ##ending with ':RabbitTableHeader' identifies a line as a row of table headers
def HeaderTranspose1(hT2):
    hT0 = '<!--\tOutput_Key\t\t-->'+Delimiter1
    hT0 += '<!--\tColumnNumber\tColumnHeader\t-->'+Delimiter1
    for iT2,nT2 in enumerate(hT2):
        nT2 = nT2.strip()
        if not(nT2.startswith('<!')):
            hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'\t-->'+Delimiter1
    return hT0+'<!--\t\t\t-->'+Delimiter1
Columns1 = ['RefName',
                 'RefPortion',  ## unique, LocalRepeat_Repeat, Public_Multimap
                 'RefNum',
                 'RefFile',
                 'RefLen',
                 'RefKCount']

sL1 = ['Sample_Info\tFiles\t\t\t\t\t',
                    'Sample_Info\tInternalIDNum\t\t\t\t\t',
                    'Sample_Info\tTotal\t0\t\t'+str(AllRefLen1)+'\t'+str(AllRefKmer1)+'\t',
                    'Sample_Info\tPassFilter\t0\t\t'+str(AllRefLen1)+'\t'+str(AllRefKmer1)+'\t',
                    'All_Reference\tUnique\tAll\t\t\t'+str(AllRefUniqueKmer1)+'\t']
if not(AggregateUniqueAndLocalRepeats1):
    sL1.append('All_Reference\tLocalRepeat\tAll\t\t\t'+str(AllRefLocalRepeatKmer1)+'\t')
if not(PrivateOnly1):
    sL1.append('All_Reference\tPublic\tAll\t\t\t'+str(AllRefPublicKmer1)+'\t')
dD1 = [[RefUniqueKmerSpecies1,
                RefUniqueKmerSpecies1,
                RefUniqueKmersFound1,
                RefUniqueReadsFound1,
                SampleUniqueKmerCovered1,
                SampleUniqueKmerMatches1,
                SampleUniqueReadMatches1,
                JUniqueKmerCounts1,
                JUniqueKmerCounts2,
                JUniqueReads1,
                JUniqueKmerCovered1],
             [RefLocalRepeatKmerPositions1,
                RefLocalRepeatKmerSpecies1,
                RefLocalRepeatKmersFound1,
                RefLocalRepeatReadsFound1,
                SampleLocalRepeatKmerCovered1,
                SampleLocalRepeatKmerMatches1,
                SampleLocalRepeatReadMatches1,
                JLocalRepeatKmerCounts1,
                JLocalRepeatKmerCounts2,
                JLocalRepeatReads1,
                JLocalRepeatKmerCovered1],
             [RefPublicKmerPositions1,
                RefPublicKmerSpecies1,
                RefPublicKmersFound1,
                RefPublicReadsFound1,
                SamplePublicKmerCovered1,
                SamplePublicKmerMatches1,
                SamplePublicReadMatches1,
                JPublicKmerCounts1,
                JPublicKmerCounts2,
                JPublicReads1,
                JPublicKmerCovered1]]    
E9 = 10**9
AllKmerAllSamples1 = sum(SamplePassFilterKmers1)
AllReadsAllSamples1 = sum(SamplePassFilterReads1)


for i1,D1 in enumerate(dD1):
    D1.append(Counter({ x : zad1(D1[7][x],D1[0][x[0]]) for x in D1[7]}))
    D1.append(Counter({ x : zad1((D1[0][x[0]]*D1[8][x]-D1[7][x]**2)**0.5,D1[0][x[0]]) for x in D1[8]}))   
    D1.append(Counter({ x : zad1(D1[10][x],D1[0][x[0]]) for x in D1[10]}))
    D1.append(Counter({ x : zad1(D1[7][x],D1[0][x[0]]) for x in D1[7]}))
    D1.append(Counter({ x : zad1(D1[9][x],SamplePassFilterReads1[x[1]]) for x in D1[9]}))
    D1.append(Counter({ x : zad1(E9*D1[7][x],D1[0][x[0]]*SampleTotalKmers1[x[1]]) for x in D1[7]}))
    D1.append(Counter({ x : zad1(E9*D1[9][x],D1[0][x[0]]*SampleTotalReads1[x[1]]) for x in D1[9]}))
    D1.append(Counter({ x : zad1(D1[7][x]*AllKmerAllSamples1,D1[2][x[0]]*SamplePassFilterKmers1[x[1]]) for x in D1[7]}))
    D1.append(Counter({ x : zad1(D1[9][x]*AllReadsAllSamples1,D1[3][x[0]]*SamplePassFilterReads1[x[1]]) for x in D1[9]}))
        
for f1,fL1 in enumerate(MetaSampleList2):
    sn1 = Samples1[f1]
    if not(SkipBlankData1) or SampleMatches1[f1]:
        Columns1.extend([x1+'__'+sn1 for x1 in TDVFields1])
        for x1 in TDVFields1:
            sL1[0] +=  ','.join([os.path.basename(x) for x in fL1])+'\t'
            sL1[1] += str(f1+1)+'\t'
            if x1 == 'ReadCount':
                sL1[2] += str(SampleTotalReads1[f1])+'\t'
                sL1[3] += str(SamplePassFilterReads1[f1])+'\t'
                sL1[4] += str(SampleUniqueReadMatches1[f1])+'\t'
                if not(AggregateUniqueAndLocalRepeats1):
                    sL1[5] += str(SampleLocalRepeatReadMatches1[f1])+'\t'
                if not(PrivateOnly1):
                    sL1[6] += str(SamplePublicReadMatches1[f1])+'\t'
            elif x1 == 'KmerCount':
                sL1[2] += str(SampleTotalKmers1[f1])+'\t'
                sL1[3] += str(SamplePassFilterKmers1[f1])+'\t'
                sL1[4] += str(SampleUniqueKmerMatches1[f1])+'\t'
                if not(AggregateUniqueAndLocalRepeats1):
                    sL1[5] += str(SampleLocalRepeatKmerMatches1[f1])+'\t'
                if not(PrivateOnly1):
                    sL1[6] += str(SamplePublicKmerMatches1[f1])+'\t'
            elif x1 == 'KmersCovered':
                TotalCovered1 = SampleUniqueKmerCovered1[f1]+SampleLocalRepeatKmerCovered1[f1]+SamplePublicKmerCovered1[f1]
                sL1[2] += str(TotalCovered1)+'\t'
                sL1[3] += str(TotalCovered1)+'\t'
                sL1[4] += str(SampleUniqueKmerCovered1[f1])+'\t'
                if not(AggregateUniqueAndLocalRepeats1):
                    sL1[5] += str(SampleLocalRepeatKmerCovered1[f1])+'\t'
                if not(PrivateOnly1):
                    sL1[6] += str(SamplePublicKmerCovered1[f1])+'\t'
            else:
                sL1[2] += '\t'
                sL1[3] += '\t'
                sL1[4] += '\t'
                if not(AggregateUniqueAndLocalRepeats1):
                    sL1[5] += '\t'
                if not(PrivateOnly1):
                    sL1[6] += '\t'
            
Header1= '\t'.join(Columns1)+'\t'+AbbrevHeader1+'\t '+Delimiter1
OutFile1.write(TaskHeader1)
OutFile1.write(HeaderTranspose1(Columns1))
OutFile1.write(Header1)
for sL0 in sL1:
    OutFile1.write(sL0+Delimiter1)
MyTDVFields0 = TDVFields0

if AggregateUniqueAndLocalRepeats1:
    PrivateUnique1 = 'Private'
else:
    PrivateUnique1 = 'Unique'
for r1 in range(1,RefCount1+1):
    if not(SkipBlankRefs1) or RefMatches1[r1]:
        if RefUniqueKmerSpecies1[r1]:
            OutItems1 = [RefName1[r1], PrivateUnique1,r1,RefFile1[r1],RefLen1[r1],RefUniqueKmerSpecies1[r1]]
            for f1 in range(SampleNum1):
                if not(SkipBlankData1) or SampleMatches1[f1]:
                    for fi1 in MyTDVFields0:
                        OutItems1.append(dD1[0][fi1][(r1,f1)])
            OutFile1.write('\t'.join(map(str,OutItems1))+'\t'+Delimiter1)
        if not(AggregateUniqueAndLocalRepeats1) and RefLocalRepeatKmerSpecies1[r1]:
            OutItems1 = [RefName1[r1], 'LocalRepeat',r1,RefFile1[r1],RefLen1[r1],RefLocalRepeatKmerPositions1[r1]]
            for f1 in  range(SampleNum1):
                if not(SkipBlankData1) or SampleMatches1[f1]:
                    for fi1 in MyTDVFields0:
                        OutItems1.append(dD1[1][fi1][(r1,f1)])
            OutFile1.write('\t'.join(map(str,OutItems1))+'\t'+Delimiter1)
        if not(PrivateOnly1) and RefPublicKmerSpecies1[r1]:
            OutItems1 = [RefName1[r1], 'Public',r1,RefFile1[r1],RefLen1[r1],RefPublicKmerPositions1[r1]]
            for f1 in  range(SampleNum1):
                if not(SkipBlankData1) or SampleMatches1[f1]:
                    for fi1 in MyTDVFields0:
                        OutItems1.append(dD1[2][fi1][(r1,f1)])
            OutFile1.write('\t'.join(map(str,OutItems1))+'\t'+Delimiter1)
myClose1(OutFile1)

    
TempCountL1 = list(dD1[0][HeatMapColorField0].values())+list(dD1[1][HeatMapColorField0].values())+list(dD1[2][HeatMapColorField0].values())
if TempCountL1:
    MaxCount1 = max(TempCountL1)
else:
    MaxCount1 = 0
tempList1 = [x for x in TempCountL1 if x>0] 
if tempList1:
    MinCount1 = min(tempList1)
else:
    MinCount1 = 0
def mydigit1(v):
    if -1<v<1:
        myAdd = 3
    else:
        myAdd = 0
    if v==0:
        return 1
    elif v>0:
        return int(abs(log(v,10)))+1+myAdd
    else:
        return int(abs(log(v,10)))+2+myAdd
def myDisplayColor1(myC,ViewDichrome=True):
    '''A color mapping of range to colors that will be visibly distinct to trichromats as well as dichromats.  ViewDichrome gives a simulated view of the output as viewed with dichromic display'''
    if myC==0 or MinCount1==0:
        return (53,53,53)
    else:
        if MaxCount1==MinCount1:
            myV = 0.5
        else:
            myV = int(log(myC/MinCount1)*1279/log(MaxCount1/MinCount1))
        if myV<256:
            myR = 0
            myG = 0
            myB = myV
        elif myV<512:
            myR = 0
            myG = myV-256
            myB = 255
        elif myV<768:
            myR = 0
            myG = 255
            myB = 768-myV
        elif myV<1024:
            myR = myV-768
            myG = 255
            myB = 0
        else:
            myR = 255
            myG = 255
            myB = myV-1024
        if not(ViewDichrome):
            return(myR,myG,myB)
        else:
            myRG = (myR+myG)//2
            return(myRG,myRG,myB)

## Three possible distance metrics for clustering
def euclid1(v1,v2):
    s1 = sum([u1**2 for u1 in v1])**0.5; s2 = sum([u2**2 for u2 in v2])**0.5
    if not(s1) or not(s2): return 2.0
    return sum([((s11/s1)-(s22/s2))**2 for s11,s22 in zip(v1,v2)])
def kljointdistance1(v1,v2,J=1):
    '''v1 and v2 are integer n-dimensional vectors, J is a regularization parameter'''
    at = sum(v1); bt = sum(v2)
    P = (at+bt)*log(at+bt+J)-at*log(at+J)-bt*log(bt+J)
    for ai,bi in zip(v1,v2):
        P += ai*log(ai+J)+bi*log(bi+J)-(ai+bi)*log(ai+bi+J)
    return P
def worstbestcaseD1(v1,v2,m=5):
    ''' v1 and v2 are vectors of integers with the same number of dimensions, m is a regularization factor'''
    b1 = sum(v1)
    b2 = sum(v2)
    bestb1 = max(0.01,b1 - m*(b1+1)**0.5)
    worstb1 = b1 + m*(b1+1)**0.5
    bestb2 = b2 + m*(b2+1)**0.5
    worstb2 = max(0.01,b2 - m*(b2+1)**0.5)
    d = 0.0
    for a1,a2 in zip(v1,v2):
        besta1 = a1 + m*(a1+1)**0.5
        worsta1 = max(0.01,a1 - m*(a1+1)**0.5)
        besta2 = max(0.01,a2 - m*(a2+1)**0.5)
        worsta2 = a2 + m*(a2+1)**0.5
        bestR = (besta1*bestb2)/(besta2*bestb1)
        worstR = (worsta1*worstb2)/(worsta2*worstb1)
        if worstR>1:
            d += log(worstR,2.0)
        elif bestR<1:
            d += -log(bestR,2.0)
    return d
RefCommonL1 =[]
for r1 in range(1,RefCount1+1):
    if RefUniqueReadsFound1[r1]+RefLocalRepeatReadsFound1[r1]*(not(AggregateUniqueAndLocalRepeats1))+RefPublicReadsFound1[r1]*(not(PrivateOnly1))>MinCountsForDisplay1*SampleNum1:
        RefCommonL1.append(r1)
if SortReferencesByReadCount1:
    RefCommonL1 = sorted(RefCommonL1, key=lambda x:RefUniqueReadsFound1[r1]+RefLocalRepeatReadsFound1[r1]*(not(AggregateUniqueAndLocalRepeats1))+RefPublicReadsFound1[r1]*(not(PrivateOnly1)))
nRefCommon1 = len(RefCommonL1)
JAllReads1 = Counter()
JAllReads1 += JUniqueReads1
if not(AggregateUniqueAndLocalRepeats1):
    JAllReads1 += JLocalRepeatReads1
if not(PrivateOnly1):
    JAllReads1 += JPublicReads1

if ClusterSamples1:
    distC1 = Counter()
    BestD1 = 99999999999999999
    for f1 in range(SampleNum1):
        myvec1 = [JAllReads1[(r1,f1)] for r1 in RefCommonL1]
        for f2 in range(SampleNum1):
            if f1==f2: continue
            myvec2 = [JAllReads1[(r1,f2)] for r1 in RefCommonL1]
            distC1[(f1,f2)] = worstbestcaseD1(myvec1,myvec2) #euclid1(myvec1,myvec2)  #kljointdistance1(myvec1,myvec2)
            if distC1[(f1,f2)]<BestD1:
                Bestf1 = f1
                BestD1 = distC1[(f1,f2)]
    treeFList1 = [Bestf1,]

    while len(treeFList1)<SampleNum1:
        bestDist1 = 99999999999999999
        for f1 in range(SampleNum1):
            if not(f1 in treeFList1) and distC1[(f1,treeFList1[-1])]<bestDist1:
                Bestf1 = f1
                bestDist1 = distC1[(f1,treeFList1[-1])]
        treeFList1.append(Bestf1)
else:
    treeFList1 = list(range(SampleNum1))
if ClusterReferences1:
    distC2 = Counter()
    BestD1 = 99999999999999999
    for r1 in RefCommonL1:
        myvec1 = [JAllReads1[(r1,f1)] for f1 in range(SampleNum1)]
        for r2 in RefCommonL1:
            if r1==r2: continue
            myvec2 = [JAllReads1[(r2,f1)] for f1 in range(SampleNum1)]
            distC2[(r1,r2)] = worstbestcaseD1(myvec1,myvec2) #euclid1(myvec1,myvec2)  #kljointdistance1(myvec1,myvec2)
            if distC2[(r1,r2)]<BestD1:
                Bestr1 = r1
                BestD1 = distC2[(r1,r2)]
    treeRList1 = [Bestr1,]

    while len(treeRList1)<len(RefCommonL1):
        bestDist1 = 99999999999999999
        for r1 in RefCommonL1:
            if not(r1 in treeRList1) and distC2[(r1,treeRList1[-1])]<bestDist1:
                Bestr1 = r1
                bestDist1 = distC2[(r1,treeRList1[-1])]
        treeRList1.append(Bestr1)
else:
    treeRList1 = RefCommonL1
if MaxRefToDisplay1:
    RefToDisplayS1 = set(sorted(treeRList1, key=lambda x: RefMatches1[x], reverse=True)[:MaxRefToDisplay1])
    treeRList1 = [x for x in treeRList1 if x in RefToDisplayS1]
x11 = 0
fontsize1 = str(vscale1)
vscale2 = vscale1
increment1 = vscale1
for i1 in range(1,len(HeatMapTextFields0)):
    increment1 = int(round(increment1/1.33))
    vscale2 += increment1
fontsize2 = str(int(vscale2*0.7))

if AggregateUniqueAndLocalRepeats1:
    RefLocalRepeatReadsFound1 = Counter()
if PrivateOnly1:
    RefPublicReadsFound1 = Counter()
    
for r1 in treeRList1[::-1]:
    for i1,myRef1 in [2,RefPublicReadsFound1[r1]],[1,RefLocalRepeatReadsFound1[r1]],[0,RefUniqueReadsFound1[r1]]:
        if not(myRef1): continue
        DigitNeed1 = mydigit1(myRef1)
        xdel1 = max(hscale1,DigitNeed1*(0.67+vscale1/2)+1.5)
        x11 -= xdel1
        y11 = 0    
        for f1 in treeFList1:
            v1 = dD1[i1][HeatMapColorField0][(r1,f1)]
            myColor1 = myDisplayColor1(v1)
            myTextColor1 = (0,0,0)
            if myColor1[1]<128: myTextColor1 = (255,255,255)
            myLabelColor1 = ['black','gray60','gray80'][i1]
            vrect(x1=x11,y1=y11,x2=x11+xdel1,y2=y11+vscale2,fill=myColor1,stroke='black',strokewidth=0.5)
            y22 = y11+vscale2+vscale1//4
            vclearance1 = int(round(vscale1))
            for p1,j1 in enumerate(HeatMapTextFields0):
                myValT1 = dis1(dD1[i1][j1][(r1,f1)])
                fontsize3 = min(vclearance1,int(1.5*(xdel1-1.5)/(len(myValT1)+0.5)))
                vtext(text=myValT1,xc=x11+xdel1/2,y1=y22-vclearance1,color=myTextColor1,font='dejavusans'+str(fontsize3))
                y22 -= vclearance1
                vclearance1 = int(round(vclearance1/1.33))
            y11 -= vscale2
        y11 -= vscale2
        vrect(x1=x11,x2=x11+xdel1,y1=y11,y2=y11+vscale2,fill='white',stroke='black',strokewidth=0.5)
        vtext(text=str(myRef1),xc=x11+xdel1/2,yc=y11+vscale2/2,color=black,font='dejavusans '+fontsize1)
        myRefName1 = RefName1[r1]
        if i1==2: myRefName1+='_Public'
        if i1==1: myRefName1+='_LocalRepeat'
        vtext(text=myRefName1,x1=x11+xdel1/2-0.22*vscale2,y1=y11-4,color=myLabelColor1,font='dejavusansbold'+fontsize2,rotate=90)
y11 = -vscale2*(len(treeFList1)+1)
y22 = y11+0.5*vscale2
x11 = 0
vclearance1 = int(round(vscale1))
myV1 = y11-2.5*sum([vclearance1*1.2**(-x) for x in range(len(HeatMapTextFields0))])
vrect(x1=x11+hscale1,x2=x11+7*hscale1,y1=y11+0.5*vscale2,y2 = myV1,fill='lightyellow',stroke='black',strokewidth=2.0)
vrect(x1=x11+hscale1,x2=x11+7*hscale1,y1=myV1-1.5*vscale2,y2=myV1-3.5*vscale2,fill='white',stroke='black',strokewidth=2.0)
vtext(text='Sample/Reference Read Counts',xc=x11+4*hscale1,yc=myV1-2.5*vscale2,color=black,font='dejavusans '+str(int(1.2*vclearance1)))
for p1,j1 in enumerate(HeatMapTextFields0):    
    fontsize3 = str(2*vclearance1)
    vtext(text=HeatMapTextFields1[p1],xc=x11+4*hscale1,y1=y22-2*vclearance1,color=black,font='dejavusansbold '+fontsize3)
    y22 -= 2.5*vclearance1
    vclearance1 = int(round(vclearance1/1.2))
y11 = 0
x11 = vscale1//2
vm1 = max([SamplePassFilterReads1[x] for x in range(SampleNum1)])
DigitNeed1 = mydigit1(vm1)
xdel1 = max(hscale1,DigitNeed1*(0.67+vscale1/2)+1.5)
if SampleCategory1:
    maxCategoryLen1 = max([len(x) for x in SampleCategory1.values()])
    CategoryBoxWidth1 = fontsize1*0.7*maxCategoryLen1
for f1 in treeFList1: 
    v1 = SamplePassFilterReads1[f1]
    vrect(x1=x11,x2=x11+xdel1,y1=y11,y2=y11+vscale2,fill='white',stroke='black',strokewidth=0.5)
    vtext(text=str(v1),xc=x11+xdel1/2,yc=y11+vscale2/2,color=black,font='dejavusans '+fontsize1)
    StartLabel1 = x11+xdel1+3
    if SampleCategory1:
        vrect(x1=x11+xdel1+3,x2=x11+xdel1+CategoryBoxWidth1+3,y1=y11,y2=y11+vscale2,fill='white',stroke='black',strokewidth=0.5)
        vtext(text=str(v1),xc=x11+xdel1+3+CategoryBoxWidth1/2,yc=y11+vscale2/2,color=SampleLabelColor1[Samples1[f1]],font='dejavusansbold '+fontsize1)
        StartLabel1 += CategoryBoxWidth1+3
    vtext(text=Samples1[f1],StartLabel1=x11+xdel1+3,yc=y11+vscale2/2,color=SampleLabelColor1[Samples1[f1]],font='dejavusansbold '+fontsize2)
    y11 -= vscale2    
vlegend(text=' ')
if MinCount1:
    myMaxKey1 = 1.0001*10**ceil(log(MaxCount1,10))
    myMinKey1 = 0.9999*10**floor(log(MinCount1,10))
    LogDivisions1 = 1+ceil(log(myMaxKey1/myMinKey1,10))
    KeyWidth1 = (2/7)*(VSG.xmax-VSG.xmin)
    AvailWidth1 = int(min(8*vclearance1, 0.3*KeyWidth1/LogDivisions1))
    fontsizeCK1 = str(max(int(fontsize1),AvailWidth1))
    vcolorkey(logmode=True, SpectrumTitle=HeatMapColorField1, mincolorindex=myMinKey1, maxcolorindex=myMaxKey1, colorvalue=myDisplayColor1, ckfont='dejavusansbold '+fontsizeCK1)
if GraphTitle1 == 'default':
    GraphTitle1 = 'SimpleSummary_'+OutFileBase1+'_'+vnow
vtitle(text=' ', font="dejasansbold "+str(2*vscale2))
vtitle(text=GraphTitle1, font="dejasansbold "+str(2*vscale2))
if not(SuppressGraphDetails1):
    vlegend(text=' ')
    vlegend(text=vSysLogInfo1, font='dejavusans 4')
vdisplay('SimpleSummary_'+OutFileBase1+'_'+vnow+'.svg')
## Copyright 2024 Andrew Fire, Stanford University, All Rights Reserved
## With thanks for ideas/prompts/encouragement/healthy-speticism (in no particular order) to K.Artiles, C.Benko, S.U.Enam, D. Galls, E.Greenwald, O.Ilbay, D.E. Jeong, J.Kim, D.Lipman, M.McCoy, M.Shoura, I.Zheludev
## Version Adjustments starting 03_24_24
## 03_24_24 aa_00 first version
