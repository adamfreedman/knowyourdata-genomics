Getting to know your data
===================

Learning Objectives:
-------------------
#### What's the goal for this lesson?
* You should be able to identify common features of "data" and data formats and what those features imply
* You should be able to look at a file and be able to identify the following types of data: fasta, fastq, fastg, sra, sff, vcf, sam, bam, bed, etc...
* You should be able to identify where your data came from (provenance)

#### At the end of this lesson you should be able to:
* Recognize file formats
* Whether your data files are zipped or unzipped. You should have an understanding of file compression.
* How to determine the file size
* How many units are in the file (i.e. nucleotides, lines of data, sequence reads, etc.)
* How to look at data structure using the shell -- does it agree with the file extension?
* Identify a binary, flat text file? 
 


**The following table shows some of the common types of data files and includes some information about them:**

| File Extension |	Type of Data |	Format |	Example(s) | 
| :------------- | :------------- | :---------------- | :----------------| :----------------| :---------------|
| txt | multi-format | Text | study metadata, tab-delimited data | <a href="https://en.wikipedia.org/wiki/Text_file" txt</a> | |
| fastq	| nucleotide  | Text |	sequencing reads |<a href="https://en.wikipedia.org/wiki/FASTQ_format" fastq </a> |  |
| fasta	| nucleotide, protein | Text | the human genome | <a href="https://en.wikipedia.org/wiki/FASTA" fasta</a>| |
| sff	| nucleotide	| Binary |	Roche/454 sequencing data |	<a href="http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=format#sff" sff</a> |	|
| vcf | multi-format | Text	 |	variation/SNP calls |	|  |
| sam | alignment | Text  |	reads aligned to a reference  | <a href="https://samtools.github.io/hts-specs/SAMv1.pdf" sam </a> |	 |
| bam | alignment	| Binary  |	reads aligned to a reference | <a href="https://samtools.github.io/hts-specs/SAMv1.pdf" bam </a> |	 |
| bed | metadata / feature definitions  | Binary  | genome coverage | <a href="http://www.ensembl.org/info/website/upload/bed.html" bed </a> |  |
| h5 | binary hierarchical | Binary | PacBio sequencing data | <a href="https://en.wikipedia.org/wiki/Hierarchical_Data_Format" h5 </a>| |
| pileup | alignment | Text | mpileup, SNP and indel calling| <a href="https://en.wikipedia.org/wiki/Pileup_format" pileup</a>| |

##Looking at the data for this workshop

In this workshop, there are a few bioinformatics-related data types we will focus on (beyond simple text files - although in principle many of the files are text). First let's consider the definition/documentation for these file types:

**Unmapped read data (fastq)**

NGS reads from a sequencing run are stored in fastq (fasta with qualities) files. For each read, there is on separate lines:
* a read identifier containing information about the instrument,sequencing run, flow cell coordinates,lane, etc.
* the nucleotide sequence
* ascii-encoded phred-scaled quality scores for each called base (phred quality = -10*log<sub>10</sub> probability of an error)
* multiple (typically 4) lines per read

To look at the structure of a fastq file, go to /n/regal/datac/fastq (NEED TO MAKE THIS MORE SPECIFIC), and look at the contents of the first fastq file with head
```
head -8  $filename.fq
```
Using this command, you will be viewing the information for the first two reads in the file.

*Important things to note:*
1. There are different character encodings of qualities. If fastq files come from different instruments/timepoints/SRA accessions,they may not have the same quality scheme, and file conversions will be necessary.
2. Quality scores vary in a semi-systematic way that depends upon things such as instrument, position along read, nucleotide context. Tomorrow, we will explore quality score distributions in the E. coli data and trim low quality bases.

More details on the fastq format are covered in the [lesson on QC of Sequence Read Data](https://github.com/devbioinfoguy/wrangling-genomics-HPC/blob/gh-pages/lessons/00-readQC.md#details-on-the-fastq-format)<br>

For the full documentation on fastq, go to [fastq wiki page](https://en.wikipedia.org/wiki/FASTQ_format)

**Aligned reads (sam)**

Alignment of fastq reads to a reference genome can be conducted with a dizzying array of alignment tools. However, these tools all use the sam format standard. For each read, a sam file (or a binary, compressed [bam](https://www.broadinstitute.org/igv/BAM) ) will contain many details, including:
* the chromosome/scaffold to which the read aligned
* the start position of that alignment
* phred-scaled quality of the alignment
* a string summarizing matches and mismatches between the read and genome to which it was mapped (CIGAR string)

|Query Name|Flag|Reference|Position|MAPQ|CIGAR|NA|NA|NA|
| ----------|-----|-------------|--------|----|-----|---|---|---|
| SRR097977.141 | 4 | * | 0 | 0 | * | * | 0 | 0 |
| SRR097977.1 | 0 | NC_012967.1 | 1366270 | 37 | 36M | * | 0 | 0 |

<b>Continued.....</b>

|Sequence|Base Qualities| 
|----------------------------------|-------------------------------------| 
|TGCCTGACCTTTCTTATGGATTTTCATTTTTTCT|C:CCCCCCCCCCCCCCC87>CC&CC9CC,??0-?|
|TATTCTGCCATAATGAAATTCGCCACTTGTTAGTGT|CCCCCCCCCCCCCCC>CCCCC7CCCCCCACA?5A5<| 

These files typically have headers that contain important information such as the sequencing strategy, sample ID, and reference genome. We can use samtools, a valuable tool for querying and viewing the contents of a sam file.  For example, to view a header (and not the reads themselves), one can cd into /n/regal/datac/precomputed/lite/variant_calling/ , and do
```
samtools view -SH /n/regal/datac/precomputed/lite/variant_calling/sam_files/SRR097977_alignment.sam
```
where, 'S' indicates the infile is sam format, the 'H' is for show header only.

To view the actual reads, one simply removes the 'H'. However....using view all by itself will read the entire (very large) sam file to standard out. So, best to pipe to head and look at the first reads.
```
samtools view -S  /n/regal/datac/precomputed/lite/variant_calling/sam_files/SRR097977_alignment.sam | head -4
```

One can also use unix command line tools to parse what comes out of samtools view. One such tool is cut. 

```
samtools view -S SRR097977_alignment.sam | cut  -f 6 |head -20
```

will output the CIGAR strings for the first 20 reads.


Awk is antother such tool. It can give you quick control over which columns you want to access. For example, to print out only the sequences, one could do:
```
samtools view -S SRR097977_alignment.sam | awk -F"\t" '{print $10}'
```

where -F specifies what the field separator is, and $10 indicates the 10th column. To print out the entire line, your print statement would be '{print $0}' .

One can also print out more than one column, for example:

* {print $1$2$3} prints the first 3 columns separated by (default) spaces
* {print $1","$4} prints the first and 4th columns separated by a comma


If you wanted to count alignments with a mapping quality>10, one could do something like this:
```
samtools view -S SRR097977_alignment.sam | awk -F"\t" '$5>10{print $0}' | wc
```

In this case, the first element of the wc command would tell you the number of reads, which should be 3866316. 

There are also methods for using bitwise flags encoded in the FLAG field to filter on such features as to whether a read is mapped...to be explored on your own at a future date! In addition, as you might suspect, you can use "|", ">", grep, and other tools you've learned today in conjunction with samtools to perform other operations.

Samtools has a lot of bells and whistles for querying and manipulating sam (and bam) files, but piping to standard out allows access to other tools for manipulating/analyzing the content of your alignment files. 

For more info on sam files, go to [sam format documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)<br>
For more on samtools command line arguments, go to [samtools manual](http://www.htslib.org/doc/samtools.html)<br>
A handy resource for awk "one-liners" can be found at [awk one liners](http://www.pement.org/awk/awk1line.txt)<br>

**Called genotypes (vcf)**
* vcf file - [definition](https://samtools.github.io/hts-specs/VCFv4.1.pdf)


**Compressed/binary**


##Exercises 

Discuss with your neighbors in conducting the exercises below.


### A. Calculating the number of reads

* How would you go about calculating the number of reads in a fastq file vs. a sam file? 

### B. How many reads are mapped? How many are unmapped?

* Pick a few sam files and calculate these quantities using one of the tools for selecting columns.
* Is there much variation among the samples?

### C. Classifying variant types

* How many variants represent mutations from A to C?
* How would you confirm that a vcf file only contains bi-allelic sites?


