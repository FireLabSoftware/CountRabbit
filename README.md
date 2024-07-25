II.  CountRabbit: 30,000 foot view
- What goes in: 
	- The names of your reference file 
	- The location of sequencing datafiles
	- <Optional> choices of parameters for 
- What comes out: 
	-A table of match counts (K-mers and/or Reads) 
	-A very quick visualization
- Important Features: Quick, Simple, Consistent handling of sequences shared between references 

-CountRabbit is K-mer based
	A-ny read that matches at least one k-mer in a reference genome is assigned to that genome
	-If only one reference genome is matched, a k-mer or read is considered to be "uniquely" mapped
	-Reads (or K-mers) that match multiple references are considered to be "public"
	-A third (optional) category "local repeat" describes K-mers 
		=present multiple times
		=but all in one reference sequence
	-With the default K-mer length of 31, specificity of matches is such that false positive will be very rare

-CountRabbit is primarily a tool for counting matches (not for many other tasks):
	-CountRabbit is not designed to 
		-provide strandedness data
		-detect distant homology, or 
		-assemble sequencing data into contigs.
	-CountRabbit is not designed for (or helpful in) detecting mutations, rearrangements, or termini
	-CountRabbit provides a provisional (and customizable) graphic output
		-But this graphic is intended for quick experimental examination
				-publication may require additional tools.
	-Features of CountRabbit overlap with
		-Many excellent and free public tools
		-Several local tools: REVA,HQAlign,Polybench, Jazz18Heap
	-If everything is working properly, CountRabbit can be much faster to implement and/or run
	-CountRabbit is a work in progress
		-For now, any results (or lack thereof) should be checked/verified with other tools
	-Please report any bugs/crashes to afire@stanford.edu
 -For more detail, see the CountRabbitDocumentationFile (CountRabbit_Docs_ba03_072424.docx) in this repository

III.  Count Rabbit: Basic Operation
Syntax
	python CountRabbit<ver>.py  Reference=<file,file..> Data=<file,file..>  
	you will need the module "VSG_ModuleFP.py" in the same folder (can obtain from Github site)
	For the program to run as fast as possible, the PyPy interpreter should be used (www.pypy.org)

Required Parameter 1: A Reference file with one or many 'reference' sequences in fasta format
	Reference=<file[s]> : Datasets to gather K-mers for analyzing data
	File lists can be comma separated (no spaces) and/or specified with a wildcard (*) 
		Wildcards require quotes around name
		Input files can be gzipped or not (or mixture). Please no spaces,commas,semicolons in filenames.
Required Parameter 2: Data Files and Sample Names
	Direct FileName Option: Data=<fasta/q file[s]>
		Here you input one or more sequencing data files for analysis
		This can be a single file or a comma-delimited list and can include wildcards (*)
	
"Data" can also be a directory
		In this case, all fasta/fastq/fastq.gz/fastq.gz in the directory will be analyzed).
	By default the sample name will be the filename 
	By default, R1/R2 data file pairs are identified by CountRabbit and treated as a single ('Rx') entity   
	SampleToData MetaFile Option
		For more flexibility, provide SampleName/ReadFileNames combinations in your own file
		This is invoked with SampleToData=<your SampleFile> instead of Data= in the command line
		SampleDescriptor files have simple line-by-line format where each line consists of
			A sample name, followed by tab, followed bya comma-delimited list of relevant the file paths 

"Most-likely-to-be-useful" parameters to explore as you evaluate CountRabbit as a tool
(all of these [and other] parameters have defaults, so you can ignore these and just run the program) 
klen: K-mer length (default=31) can be anything >~16; values >31 will slow CountRabbit down
OutFileName: Default=<assigned by program > : Allows user to specify base filenames for output. 
MaxReadsPerFile: Max reads to check for each Data file (default is none: all reads processed)
CircularReference: (Default=False) Setting to true will treat reference sequences as circles.  
ReferenceFileByFile: (Default=False) Setting to true treats each reference file as a single entity.  
Express  (Default=False)=True runs in mode that identifies just the first-K-mer match for each read 
FilterList: Providing a fasta file and a cutoff number (e.g. myFilterFile.fa>5) limits the searched K-mers based on K-mer copy number in a separate reference file
Threaded: (Default=SingleThread)  Threaded=True instructs rabbit to use a multi-processor mode and decide how many threads to use.  Of you can specify a number (e.g. Threads=3 uses three threads).  Multithreading is only useful if you have multiple input samples.
SingleCount: (Default=False)  This will run CountRabbit in a very limited but quick mode where only a single read count is reported for each reference and dataset.

IV. Packing List:
A full package for this project consists of a zipped folder (CountRabbitFolder) with the following contents:
1.  The CountRabbit python script (CountRabbit_ver_date.py)
2.  The needed VSG_ModuleFP.py script
3.  This documentation (CountRabbit_Docs_ba03_072324.docx)
4.  A useful junk-linker file (illuminadetritis1223wMultiN.fa)
5.  A subdirectory containing a set of test data and references sequences (CountRabbitTestData: should contain 8 .fastq.gz files plus a cerep.fa reference file)
The following script, if run in the CountRabbitFolder directory should produce usable test output and a nice graphic within a few seconds.  If not, give us a holler.

pypy3 CountRabbit_ba03_072324.py Reference='./CountRabbitTestData/celrep.fa' Data='./CountRabbitTestData/*.fastq.gz' DetritisF='./illuminadetritis1223wMultiN.fa'
