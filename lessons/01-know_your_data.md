#Getting to know your data
===================

##Learning Objectives:
-------------------

* You should be able to identify common genomics data formats
* You should be able to identify where your data came from (provenance)
* To understand in more detail fastq, sam, and vcf formats


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

##Unmapped read data (FASTQ)

NGS reads from a sequencing run are stored in fastq (fasta with qualities). Although it looks complicated  (and maybe it is), its easy to understand the [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format with a little decoding. Some rules about the format include...

|Line|Description|
|----|-----------|
|1|Always begins with '@' and then information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+' and sometimes the same info in line 1|
|4|Has a string of characters which represent the quality scores; must have same number of characters as line 2|

so for example in our data set, one complete read is:
```
$ head -n4 ~/dc_sample_data/untrimmed_fastq/SRR098281.fastq 
@SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
CNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+SRR098281.1 HWUSI-EAS1599_1:2:1:0:318 length=35
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
This is a pretty bad read. 

Notice that line 4 is:
```
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
As mentioned above, line 4 is a encoding of the quality. In this case, the code is the [ASCII](https://en.wikipedia.org/wiki/ASCII#ASCII_printable_code_chart) character table. According to the chart a '#' has the value 35 and '!' has the value 33. If only it were that simple. There are actually several historical differences in how Illumina and other players have encoded the scores. Heres the chart from wikipedia:

```
  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
  0........................26...31.......40                                
                           -5....0........9.............................40 
                                 0........9.............................40 
                                    3.....9.............................40 
  0.2......................26...31........41                              

 S - Sanger        Phred+33,  raw reads typically (0, 40)
 X - Solexa        Solexa+64, raw reads typically (-5, 40)
 I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
 J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
     (Note: See discussion above).
 L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
 ```
 So using the Illumina 1.8 encouding, which is what you will mostly see from now on, our first c is called with a Phred score of 0 and our Ns are called with a score of 2. Read quality is assessed using the Phred Quality Score.  This score is logarithmically based and the score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------|:---------------------------------:|-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|
|50	|1 in 100,000|	99.999%|
|60	|1 in 1,000,000|	99.9999%|

* a read identifier containing information about the instrument,sequencing run, flow cell coordinates,lane, etc.
* the nucleotide sequence
* ascii-encoded phred-scaled quality scores for each called base (phred quality = -10*log<sub>10</sub> probability of an error)
* multiple (typically 4) lines per read
* There are different character encodings of qualities. If fastq files come from different instruments/timepoints/SRA accessions, may be an issue
* Quality scores vary in a semi-systematic way that depends upon things such as instrument, position along read, nucleotide context.


##Aligned reads (sam)

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

These files typically have headers that contain important information such as the sequencing strategy, sample ID, and reference genome. We can use samtools, a valuable tool for querying and viewing
 the contents of a sam file.<br>

For example, to view a header (and not the reads themselves), cd into `$PRECOMPUTED/lite/variant_calling/samfiles/`, and try the following:
 
```bash
samtools view -SH SRR097977_alignment.sam
```
where, 'S' indicates the infile is sam format, the 'H' is for show header only.

To view the actual reads, one simply removes the 'H'. However....using view all by itself will read the entire (very large) sam file to standard out. So, best to pipe to head and look at the first reads.
```bash
samtools view -S SRR097977_alignment.sam | head -4
```

One can also use unix command line tools to parse what comes out of samtools view. One such tool is cut. 

```bash
samtools view -S SRR097977_alignment.sam | cut  -f 6 |head -20
```

will output the CIGAR strings for the first 20 reads.


Awk is antother such tool. It can give you quick control over which columns you want to access. For example, to print out only the sequences, one could do:
```bash
samtools view -S SRR097977_alignment.sam | awk -F"\t" '{print $10}'
```

where -F specifies what the field separator is, and $10 indicates the 10th column. To print out the entire line, your print statement would be '{print $0}' .

One can also print out more than one column, for example:

* {print $1$2$3} prints the first 3 columns separated by (default) spaces
* {print $1","$4} prints the first and 4th columns separated by a comma


If you wanted to count alignments with a mapping quality greater or equal to 10, one could do something like this:
```bash
samtools view -S SRR097977_alignment.sam | awk -F"\t" '$5>=10{print $0}' | wc
```

In this case, the first element of the wc command would tell you the number of reads, which should be 3866316. 

We could also determine how many sites have a mapping quality equal to 37:
```bash
samtools view -S SRR097977_alignment.sam | awk -F"\t" '$5==37{print $0}' | wc
```
Note the "==" notation. This is the logical statement version of "are equal." 

There are also methods for using bitwise flags encoded in the FLAG field to filter on such features as 
to whether a read is mapped...to be explored on your own at a future date! In addition, as you might suspect,
 you can use "|", ">", grep, and other tools you've learned today in conjunction with samtools to perform other operations.

Samtools has a lot of bells and whistles for querying and manipulating sam (and bam) files, but piping to standard out allows
 access to other tools for manipulating/analyzing the content of your alignment files. 

For more info on sam files, go to [sam format documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)<br>
For more on samtools command line arguments, go to [samtools manual](http://www.htslib.org/doc/samtools.html)<br>
A handy resource for awk "one-liners" can be found at [awk one liners](http://www.pement.org/awk/awk1line.txt)<br>

##Called genotypes (vcf)

Variant call format (vcf) files provide detailed information regarding the genomic position of called variants,
the reference nucleotide at the position of interest, the mutation, measures of confidence in the variant call,
and a diversity of other features including custom annotations. Genotyping tools can set to only output variant site
information, or to also include invariant sites (useful for conducting sliding window analyses).<br>

An important part of vcf files is the header section, that describes all of the fields. 

Go ahead and cd into `variant_calling/vcfs` and look at the header section (lines commented out with ##), and scroll down until you get to the
column labels, and the first few genotypes. Using less, you can scroll up or down with the up and down arrow keys.
```bash
less $PRECOMPUTED/lite/variant_calling/vcfs/SRR097977_alignment_varfltrd.vcf
```

* for more information on the vcf format - [definition](https://samtools.github.io/hts-specs/VCFv4.1.pdf)


##Exercises 

Discuss with your neighbors in conducting the exercises below.


### A. Calculating the number of reads

* How would you go about calculating the number of reads in a fastq file vs. a sam file? 

### B. How many reads are mapped? How many are unmapped?

* Pick a few sam files and calculate these quantities using one of the tools for selecting columns.
* Is there much variation among the samples?

### C. Classifying variant types

* How many variants represent mutations from A to C? Hint: you can link a sequence of queries on particular columns with '|' .



