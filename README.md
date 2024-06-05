# KMerCountRabbit
Fast Read Positioning and KMer counting program
########################
## KMerCountRabbit ('KCR')-- Very Simple and Quick Reference-based read and kmer counting for NGS Datasets
##
## OVERVIEW OF KCR
## - You have a bunch of datasets from next generation sequencing and a reference, and want to count matches between them
## - KCR is a tool providing tables of read match counts between sequence datan(fasta or fastq read files) and reference genome(s)
## - KCR has some shared capabilities with standard tools like BLAST/Bowtie/Star/BWA and to many home-grown tools
## - The goals with KCR are to provide a potentially simpler worflow for analysis
## -   What goes in: The names of your reference file and the location of sequencing datafiles
## -   What comes out: A table of match counts (Kmers and Reads)
## -   Important Features: Quick, Simple, Consistent handling of sequences shared between references 
## - KCR is K-mer based
## -   Any read that matches at least one k-mer in a reference genome is assigned to that genome
## -   If no other reference genome is matched, a k-mer or read is considered to be "uniquely" mapped
## -   Reads (or K-mers) that match multiple references are considered to be "public" and can match several references
## -   With the default K-mer length of 31, the specificity of matches is such that false positive should be vanishingly rare
## - KCR is primarily a tool for counting matches (not for many other tasks):
## -   KCR is not designed to provided strandedness data, detect distant homology, or assemble sequencing data into contigs.
## -   KCR is not designed for (or helpful in) detecting mutations, rearrangements, or termini
## -   KCR provides a provisional (and customizable) graphic output, but this graphic is intended for quick experimental examination, while publication may require additional tools.
## - Specific features of KCR overlap with local tools REVA, HQAlign, Polybench, and Jazz18Heap, and many public tools
## - If everything is working properly, KCR can be much faster to implement and/or run than any of these
## - KCR is very much a work in progress- not recommended for production until we've tested it more thoroughly
## - Report any crashes or anomalous results to (afire@stanford.edu) 
##
## BASIC OPERATION
## Version ah00 06_04_2024
## ->Fuction: KCR Takes a reference (fasta) file and a set of data (fasta or fastq) files with multiple contigs or segments and returns a table of k-mer counts
##
## ->Syntax
##  - python KMerCountRabbit<ver>.py  Data=<file,file..>  Reference=<file,file..>
##  -   you will need the module "VSG_ModuleFL.py" in the same folder (also on GITHUB Site)
##  - For the program to run as fast as possible, the PyPy interpreter should be used (www.pypy.org)
##
## ->Required Parameters
##  - Data=<file[s]> : Starting dataset(s) (reads from these files are analyzed to generate your output).  This can also include a wildcard ('*') or be a directory (in which case all fasta/fastq files, gzipped or not will be analyzed).
##  -   Rabbit will try to pair R1 and R2 files using the names, reporting aggregate results for each pair.  To turn off pairing, indicate AutoPair=False
##  -   As an alternative to a set of individual files, the user can specify a file that lists sample names and associated files in a line-by-line format
##  -   The syntax for this is "SampleToData=<myfile>" where <myfile> is a simple text file with each line starting with a sample name, followed by a tab and a comma-delimited list of relevant file paths for sequence files with data for that sample
##  - Reference=<file[s]> : Datasets to gather kmers for analyzing data
##  -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
##  -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
##
## ->Some of the many parameters that may be worth exploring sooner rather than later (all have defaults, so you can ignore all of these and just run the program) 
##  - klen : is the Kmer length (default 31, realistically can be between 16 and read-length; values of 31 will slow this down). Default=31
##  - OutFileName : name assigned by program if not specified> : Allows user to specify base filenames for data output. Default=<autoassigned name> 
##  - MaxReadsPerFile : Max number to check for each Data file (default is zero (all reads).  Default=0 (process all reads)
##  - CircularReference : Setting this to true will instruct KFR to treat reference sequences as circles.  Default=False (treat as linear)
##  - ReferenceFileByFile : Each reference file is treated as a single entity (even if there are many contigs in the file-- e.g. for a partial assembly).  Default=False (each fasta entry treated as a separate entity)
##  - Express : Setting this to True runs the program in an "Express" mode which quickly provides basic first-matching-reference read counts much more quickly (Default=False)
##  - FilterList : Providing a fasta file and a cutoff number (e.g. myFilterFile.fa>5) limits the searched kmers to those whose copy number in the filter file matches the inequality (in this case >5)
## -> Many other settings and options are described in the detailed instructions (start of script). Some still being debugged-- so caveat mappor
#################
