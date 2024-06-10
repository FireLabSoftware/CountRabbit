KMerCountRabbit ('KCR')-- Very Simple and Quick Reference-based read and kmer counting for NGS Datasets

Version aj00 06_09_2024

OVERVIEW OF KCR
- You have a bunch of datasets from next generation sequencing and a reference, and want to count matches between them
- KCR is a tool providing tables of read match counts between sequence datan(fasta or fastq read files) and reference genome(s)
- KCR has some shared capabilities with standard tools like BLAST/Bowtie/Kallisto/BWA/Kraken and to many home-grown tools
- The goals with KCR are to provide a potentially simpler worflow for analysis
-   What goes in: The names of your reference file and the location of sequencing datafiles
-   What comes out: A table of match counts (Kmers and/or Reads) and (for preview purposes) a very quick visualization
-   Important Features: Quick, Simple, Consistent handling of sequences shared between references 
- KCR is K-mer based
-   Any read that matches at least one k-mer in a reference genome is assigned to that genome
-   If only one reference genome is matched, a k-mer or read is considered to be "uniquely" mapped
-   Reads (or K-mers) that match multiple references are considered to be "public"
-   A third (optional) category called "local repeat" describes kmers present multiple times, but all in one reference sequence
-   With the default K-mer length of 31, the specificity of matches is such that false positive should be vanishingly rare
- KCR is primarily a tool for counting matches (not for many other tasks):
-   KCR is not designed to provided strandedness data, detect distant homology, or assemble sequencing data into contigs.
-   KCR is not designed for (or helpful in) detecting mutations, rearrangements, or termini
-   KCR provides a provisional (and customizable) graphic output, but this graphic is intended for
-     quick experimental examination, while publication may require additional tools.
- Specific features of KCR overlap with local tools REVA,HQAlign,Polybench, Jazz18Heap, and many public tools
- If everything is working properly, KCR can be much faster to implement and/or run than any of these
- KCR is very much a work in progress- not recommended for production until we've tested it more thoroughly
- Report any crashes or anomalous results to (afire@stanford.edu) 

BASIC OPERATION

->Syntax
 - python KMerCountRabbit<ver>.py  Reference=<file,file..> Data=<file,file..>  
 -   you will need the module "VSG_ModuleFM.py" in the same folder (VSG_ModuleFM.py is on GITHUB Site)
 - For the program to run as fast as possible, the PyPy interpreter should be used (www.pypy.org)

->Required Parameter 1: A Reference file with one or many 'reference' sequences in fasta format
 - Reference=<file[s]> : Datasets to gather kmers for analyzing data
 -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
 -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.

->Required Parameter 2: Data Files and Sample Names
 - Direct FileName Option: Data=<fasta/q file[s]> : Here you input one or more sequencing data files for analysis
 -   This can be a single file or a comma-delimited list and can include wildcards (*)
 -   "Data" can also be a directory (in this case, all fasta/fastq/fastq.gz/fastq.gz in the directory will be analyzed).
 -   By default the sample name will be the filename and Rabbit "AutoPair" R1/R2 data file pairs as a single ('Rx') entity   
 - SampleToData MetaFile Option: For more flexibility, provide SampleName/ReadFileNames combinations in your own file
 -   This is invoked with SampleToData=<your SampleDescriptorFile> instead of Data= in command line
 -   SampleDescriptor files have simple line-by-line format where each line consists of
 -   A sample name, followed by tab, followed bya comma-delimited list of relevant the file paths with data from that sample

->Some additional parameters to explore as you evaluate KCR as a tool
 -    (all of these [and other] parameters have defaults, so you can ignore all of these and just run the program) 
 - klen : is the Kmer length (default 31, realistically can be between 16 and read-length; values of 31 will slow this down). Default=31
 - OutFileName : name assigned by program if not specified> : Allows user to specify base filenames for data output. Default=<autoassign> 
 - MaxReadsPerFile : Max number to check for each Data file (default is zero (all reads).  Default=0 (process all reads)
 - CircularReference : Setting this to true will instruct KFR to treat reference sequences as circles.  Default=False (treat as linear)
 - ReferenceFileByFile : Setting to true treats each reference file as a single entity.  (Default=Each fasta entry is an entity)
 - Express : Setting this to True runs the program in mode which just yields first-matching-reference read counts (Default=False)
 - FilterList : Providing a fasta file and a cutoff number (e.g. myFilterFile.fa>5) limits the searched kmers
 -   based on kmer copy number in a separate reference file
-> Many other settings and options are described in the detailed instructions (start of script).
 - Some still being debugged-- so caveat mappor
