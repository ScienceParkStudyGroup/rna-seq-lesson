---
title: "From fastq files to read counts"
teaching: 45
exercises: 0 
questions:
- "How do I perform a quality check of my RNA-seq fastq files with `FastQC`?"
- "How can I remove RNA-seq reads of low quality? using `trimmomatic`?"
- "How do I align my reads to a reference genome using `STAR`?"
- "What is the SAM/BAM format?"
- "How do I turn RNA-seq read genome alignments into a count table?"
objectives:
- "Be able to remove RNA-seq reads with adapters and low quality bases."
keypoints:
- "Performing read trimming ensures that no sequencing adapter is left in your final reads."
---

# Table of contents
<!-- MarkdownTOC autolink="True" levels="1,2" -->

- [1 Quality control and trimming](#1-quality-control-and-trimming)
  - [1.1 Running FastQC](#11-running-fastqc)
  - [1.2 Viewing the FastQC results](#12-viewing-the-fastqc-results)
  - [1.3 Decoding the other FastQC outputs](#13-decoding-the-other-fastqc-outputs)
  - [1.4 Working with the FastQC text output](#14-working-with-the-fastqc-text-output)
  - [1.5 Documenting our work](#15-documenting-our-work)
- [2. Trimming and filtering](#2-trimming-and-filtering)
- [3. Alignment to a reference genome](#3-alignment-to-a-reference-genome)
  - [3.1 Setting up](#31-setting-up)
  - [3.2 Index the reference genome](#32-index-the-reference-genome)
  - [3.3 Align reads to reference genome](#33-align-reads-to-reference-genome)
  - [3.4 The SAM/BAM format](#34-the-sambam-format)
- [4. Creating the counts file](#4-creating-the-counts-file)

<!-- /MarkdownTOC -->
  

<img src="../img/04-workflow-overview.png" width="500px" alt="workflow overview">

# 1 Quality control and trimming


## 1.1 Running FastQC  

We will now assess the quality of the reads that we downloaded. First, we need to make an output directory for the fastqc results to be stored. This we want to do in the 'home' directory that contains all the needed files.

~~~
$ docker run -it --name bioinfo scienceparkstudygroup/master-gls:fastq-latest

$ conda activate fastq

$ cd home

$ mkdir fastqc
~~~
{: .bash}


Next we need to get to the directory thay actually contains the the fastq files.

~~~
$ ls fq-files
~~~
{: .bash}


Running fastqc uses the following command

~~~
fastqc -o fastqc Arabidopsis_sample1.fq.gz
~~~
{: .bash}

Of course we don't want to do this for all the samples seperately so we can loop through the list of samples and run them all sequentially

With the use of echo you can start off with a "dry run"

~~~
$ for filename in  *.fq.gz
  do
    echo fastqc -o fastqc $filename
  done
~~~
{: .bash}

The echo command only prints the commands to the screen, and doesn't really run it.

~~~
fastqc -o fastqc Arabidopsis_sample1.fq.gz
fastqc -o fastqc Arabidopsis_sample2.fq.gz
fastqc -o fastqc Arabidopsis_sample3.fq.gz
fastqc -o fastqc Arabidopsis_sample4.fq.gz
~~~
{: .output}

If it looks good remove the echo and go for it.

~~~
$ for filename in *.fq.gz
  do
    fastqc -o fastqc $filename
  done
~~~
{: .bash}

You will see an automatically updating output message telling you the
progress of the analysis. It shoul look something like this:

~~~
Started analysis of Arabidopsis_sample1.fq.gz
Approx 5% complete for Arabidopsis_sample1.fq.gz
Approx 10% complete for Arabidopsis_sample1.fq.gz
Approx 15% complete for Arabidopsis_sample1.fq.gz
Approx 20% complete for Arabidopsis_sample1.fq.gz
…
Approx 85% complete for Arabidopsis_sample4.fq.gz
Approx 90% complete for Arabidopsis_sample4.fq.gz
Approx 95% complete for Arabidopsis_sample4.fq.gz
Approx 100% complete for Arabidopsis_sample4.fq.gz
Analysis complete for Arabidopsis_sample4.fq.gz
~~~
{: .output}

In total, it should take about five minutes for FastQC to run on all
four of our zipped FASTQ files.

If the command doesn't run or you want more information on fastqc, run the following to get the help page.

~~~
$ fastqc -h
~~~
{: .bash}

But if all went right, the FastQC program will have created several new files within our
`~/home/fastqc` directory.

~~~
$  ls fastqc/
~~~
{: .bash}

~~~
Arabidopsis_sample1_fastqc.html  Arabidopsis_sample2_fastqc.zip   Arabidopsis_sample4_fastqc.html
Arabidopsis_sample1_fastqc.zip	 Arabidopsis_sample3_fastqc.html  Arabidopsis_sample4_fastqc.zip
Arabidopsis_sample2_fastqc.html  Arabidopsis_sample3_fastqc.zip
~~~
{: .output}


## 1.2 Viewing the FastQC results

For each of the samples there are two files. a .html and a .zip 

If we were working on our local computer, outside of the container, we'd be able to display each of these
HTML files as a webpage:

If we were working on our local computers, we'd be able to display each of these
HTML files as a webpage:

~~~
$ cd fastqc/
$ open Arabidopsis_sample1_fastqc.html
~~~
{: .bash}


However, if you try this in the docker container we're working in, you'll get an error:

~~~
bash: open: command not found
~~~
{: .output}

This is because the container were working in doesn't have any web
browsers installed on it, so the remote computer doesn't know how to
open the file. We want to look at the webpage summary reports, so
let's transfer them to our local computers (i.e. your laptop).

If you're also working on a remote computer you will first have to copy 
the files outside of the container using `docker cp` and next from the 
remote computer to your local computer with the help of `scp`.

First we need to exit the container and next we can transfer our HTML 
files to our local computer using `docker cp`.

~~~
$ mkdir -p ~/Desktop/fastqc_html
~~~
{: .bash}


~~~
$ docker cp bioinfo:/home/fastqc/ ~/Desktop/fastqc_html
~~~
{: .bash}

This will transfer all files in the folder home/fastqc/ to your Desktop.

bioinfo refers to the name of the container.
The second part starts with a `:` and then gives the absolute path of the files you want to transfer.
The third part of the command gives the absolute path of the location
you want to put the files. This is on your local computer and is the
directory we just created `~/Desktop/fastqc_html`.


When working on a remote computer make use of the following command

~~~
$ scp -r tbliek@genseq-cn02.science.uva.nl:~/rna_seq_lesson/fastqc/ ~/Desktop/fastqc_html
~~~
{: .bash}


As a reminder, the first part
of the command `tbliek@genseq-cn02.science.uva.nl` is
the address for your remote computer. Make sure you replace everything
after `dcuser@` with your instance number (the one you used to log in).

The second part starts with a `:` and then gives the absolute path
of the files you want to transfer from your remote computer. Don't
forget the `:`. We used a wildcard (`*.html`) to indicate that we want all of
the HTML files.

The third part of the command gives the absolute path of the location
you want to put the files. This is on your local computer and is the
directory we just created `~/Desktop/fastqc_html`.

The -r option is used to tell scp to recursively copy the source directory and its contents.

You should see a status output like this:

~~~
Arabidopsis_sample1_fastqc.html                      100%  249KB 152.3KB/s   00:01    
Arabidopsis_sample1_fastqc.html                       100%  254KB 219.8KB/s   00:01      
~~~
{: .output}

Now we can go to our new directory and open the HTML files.

~~~
$ cd ~/Desktop/fastqc_html/
$ open *.html
~~~
{: .bash}


Your computer will open each of the HTML files in your default web browser. Depending on your settings, this might be as six separate
tabs in a single window or six separate browser windows.

 
## 1.3 Decoding the other FastQC outputs
We've now looked at quite a few "Per base sequence quality" FastQC graphs, but there are nine other graphs that we haven't talked about! Below we have provided a brief overview of interpretations for each of these plots. It's important to keep in mind

+ **Per tile sequence quality**: the machines that perform sequencing are divided into tiles. This plot displays patterns in base quality along these tiles. Consistently low scores are often found around the edges, but hot spots can also occur in the middle if an air bubble was introduced at some point during the run.
+ **Per sequence quality scores**: a density plot of quality for all reads at all positions. This plot shows what quality scores are most common.
+ **Per base sequence content**: plots the proportion of each base position over all of the reads. Typically, we expect to see each base roughly 25% of the time at each position, but this often fails at the beginning or end of the read due to quality or adapter content.
+ **Per sequence GC content**: a density plot of average GC content in each of the reads.  
+ **Per base N content**: the percent of times that 'N' occurs at a position in all reads. If there is an increase at a particular position, this might indicate that something went wrong during sequencing.  
+ **Sequence Length Distribution**: the distribution of sequence lengths of all reads in the file. If the data is raw, there is often on sharp peak, however if the reads have been trimmed, there may be a distribution of shorter lengths.
+ **Sequence Duplication Levels**: A distribution of duplicated sequences. In sequencing, we expect most reads to only occur once. If some sequences are occurring more than once, it might indicate enrichment bias (e.g. from PCR). If the samples are high coverage (or RNA-seq or amplicon), this might not be true.  
+ **Overrepresented sequences**: A list of sequences that occur more frequently than would be expected by chance.
+ **Adapter Content**: a graph indicating where adapater sequences occur in the reads.

## 1.4 Working with the FastQC text output

Now that we've looked at our HTML reports to get a feel for the data,
let's look more closely at the other output files. Go back to the tab
in your terminal program that is connected to your AWS instance
(the tab lab will start with `dcuser@ip`) and make sure you're in
our results subdirectory.   

~~~
$ cd home/fastqc/
$ ls
~~~
{: .bash}

~~~
Arabidopsis_sample1_fastqc.html  Arabidopsis_sample2_fastqc.zip   Arabidopsis_sample4_fastqc.html
Arabidopsis_sample1_fastqc.zip	 Arabidopsis_sample3_fastqc.html  Arabidopsis_sample4_fastqc.zip
Arabidopsis_sample2_fastqc.html  Arabidopsis_sample3_fastqc.zip
~~~
{: .output}

Our `.zip` files are compressed files. They each contain multiple
different types of output files for a single input FASTQ file. To
view the contents of a `.zip` file, we can use the program `unzip`
to decompress these files. Let's try doing them all at once using a
wildcard.

~~~
$ conda install -c conda-forge unzip
$ unzip *.zip
~~~
{: .bash}

~~~
Archive:  Arabidopsis_sample1_fastqc.zip
caution: filename not matched:  Arabidopsis_sample2_fastqc.zip
caution: filename not matched:  Arabidopsis_sample3_fastqc.zip
caution: filename not matched:  Arabidopsis_sample4_fastqc.zip
~~~
{: .output}

This didn't work. We unzipped the first file and then got a warning
message for each of the other `.zip` files. This is because `unzip`
expects to get only one zip file as input. We could go through and
unzip each file one at a time, but this is very time consuming and
error-prone. Someday you may have 500 files to unzip!

A more efficient way is to use a `for` loop like we learned in the Shell Genomics lesson to iterate through all of
our `.zip` files. Let's see what that looks like and then we'll
discuss what we're doing with each line of our loop.

~~~
$ for filename in *.zip
> do
>  unzip $filename
> done
~~~
{: .bash}


In this example, the input is four filenames (one filename for each of our `.zip` files).
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `unzip` command.
The first time through the loop,
`$filename` is `Arabidopsis_sample1_fastqc.zip`.
The interpreter runs the command `unzip` on `Arabidopsis_sample1_fastqc.zip`.
For the second iteration, `$filename` becomes
`Arabidopsis_sample2_fastqc.zip`. This time, the shell runs `unzip` on `Arabidopsis_sample2_fastqc.zip`.
It then repeats this process for the four other `.zip` files in our directory.


When we run our `for` loop, you will see output that starts like this:

~~~
Archive:  fastqc/Arabidopsis_sample2_fastqc.zip
   creating: Arabidopsis_sample2_fastqc/
   creating: Arabidopsis_sample2_fastqc/Icons/
   creating: Arabidopsis_sample2_fastqc/Images/
  inflating: Arabidopsis_sample2_fastqc/Icons/fastqc_icon.png  
  ...
  inflating: Arabidopsis_sample4_fastqc/Images/duplication_levels.png  
  inflating: Arabidopsis_sample4_fastqc/Images/adapter_content.png  
  inflating: Arabidopsis_sample4_fastqc/fastqc_report.html  
  inflating: Arabidopsis_sample4_fastqc/fastqc_data.txt  
  inflating: Arabidopsis_sample4_fastqc/fastqc.fo 
~~~
{: .output}

The `unzip` program is decompressing the `.zip` files and creating
a new directory (with subdirectories) for each of our samples, to
store all of the different output that is produced by FastQC. There
are a lot of files here. The one we're going to focus on is the
`summary.txt` file.

If you list the files in our directory now you will see:

~~~
Arabidopsis_sample1_fastqc     	 Arabidopsis_sample2_fastqc.html  Arabidopsis_sample3_fastqc.zip
Arabidopsis_sample1_fastqc.html  Arabidopsis_sample2_fastqc.zip   Arabidopsis_sample4_fastqc
Arabidopsis_sample1_fastqc.zip	 Arabidopsis_sample3_fastqc	  Arabidopsis_sample4_fastqc.html
Arabidopsis_sample2_fastqc	 Arabidopsis_sample3_fastqc.html  Arabidopsis_sample4_fastqc.zip
~~~
{: .output}


The `.html` files and the uncompressed `.zip` files are still present,
but now we also have a new directory for each of our samples. We can
see for sure that it's a directory if we use the `-F` flag for `ls`.

~~~
$ ls -F
~~~
{: .bash}

~~~
Arabidopsis_sample1_fastqc/	 Arabidopsis_sample2_fastqc.html  Arabidopsis_sample3_fastqc.zip   
Arabidopsis_sample1_fastqc.html  Arabidopsis_sample2_fastqc.zip   Arabidopsis_sample4_fastqc/
Arabidopsis_sample1_fastqc.zip	 Arabidopsis_sample3_fastqc/	  Arabidopsis_sample4_fastqc.html
Arabidopsis_sample2_fastqc/	 Arabidopsis_sample3_fastqc.html  Arabidopsis_sample4_fastqc.zip
~~~
{: .output}

Let's see what files are present within one of these output directories.

~~~
$ ls -F Arabidopsis_sample1_fastqc/
~~~
{: .bash}

~~~
fastqc_data.txt  fastqc.fo  fastqc_report.html	Icons/	Images/  summary.txt
~~~
{: .output}

Use `less` to preview the `summary.txt` file for this sample.

~~~
$ less Arabidopsis_sample1_fastqc/summary.txt
~~~
{: .bash}

~~~
PASS    Basic Statistics        Arabidopsis_sample1.fq.gz
PASS    Per base sequence quality       Arabidopsis_sample1.fq.gz
PASS    Per sequence quality scores     Arabidopsis_sample1.fq.gz
FAIL    Per base sequence content       Arabidopsis_sample1.fq.gz
PASS    Per sequence GC content Arabidopsis_sample1.fq.gz
PASS    Per base N content      Arabidopsis_sample1.fq.gz
PASS    Sequence Length Distribution    Arabidopsis_sample1.fq.gz
WARN    Sequence Duplication Levels     Arabidopsis_sample1.fq.gz
WARN    Overrepresented sequences       Arabidopsis_sample1.fq.gz
PASS    Adapter Content Arabidopsis_sample1.fq.gz
~~~
{: .output}

The summary file gives us a list of tests that FastQC ran, and tells
us whether this sample passed, failed, or is borderline (`WARN`). Remember to quit from `less` you enter `q`.

## 1.5 Documenting our work

We can make a record of the results we obtained for all our samples
by concatenating all of our `summary.txt` files into a single file
using the `cat` command. We'll call this `full_report.txt` and move
it to `~/dc_workshop/docs`.

~~~
$ cat */summary.txt > fastqc_summaries.txt
~~~
{: .bash}

# 2. Trimming and filtering

Before we will do the alignment we need to remove sequences of low quality and sequences that are to short (below 25 bases).
Also in this case we will trim down long sequences to 100 bases, quality of the Ion-torrent reads drops the further it gets.
When making use of illumina reads this is not as much of a problem and 3'-trimming would then be a waste of data.

To start off make a directory trimmed for the output and then back to the rawReads directory.

~~~
$ cd /home/
$ mkdir trimmed
~~~
{: .bash}



The trimming and quality filtering will be done with **trimmomatic**.
In the programm the following arguments can be used.

| step   | meaning |
| ------- | ---------- |
| `SE` or `PE` | Reads are single end or paired end. |
| `ILLUMINACLIP` | Perform adapter removal |
| `SLIDINGWINDOW` | Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`  | Cut bases off the start of a read, if below a threshold quality.  |
|  `TRAILING` |  Cut bases off the end of a read, if below a threshold quality. |
| `CROP`  |  Cut the read to a specified length. |
|  `HEADCROP` |  Cut the specified number of bases from the start of the read. |
| `MINLEN`  |  Drop an entire read if it is below a specified length. |
|  `TOPHRED33` | Convert quality scores to Phred-33.  |
|  `TOPHRED64` |  Convert quality scores to Phred-64. |


To run this on a single sample it looks something like this
~~~
trimmomatic SE -phred33 -threads 2 Arabidopsis_sample1.fq.gz trimmed/Arabidopsis_sample1_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
~~~
{: .bash}


Of cource we don't want to do this for all the reads seperately so lets create a loop through all the fastq files.

When doing the fastqc only input files needed to be specified. In this case both the input and a matching output filenames need to be given.
this can be done with the help of 'basename'


~~~
$ for infile in *.fq.gz
do
 echo inputfile $infile
 outfile="$(basename $infile .fq.gz)"_qc.fq
 echo outputfile $outfile
 echo
done
~~~
{: .bash}


This be be producing the following list

~~~
inputfile Arabidopsis_sample1.fq.gz
outputfile Arabidopsis_sample1_qc.fq

inputfile Arabidopsis_sample2.fq.gz
outputfile Arabidopsis_sample2_qc.fq

inputfile Arabidopsis_sample3.fq.gz
outputfile Arabidopsis_sample3_qc.fq

inputfile Arabidopsis_sample4.fq.gz
outputfile Arabidopsis_sample4_qc.fq
~~~
{: .output}

Next we can start writing the trimmomatic loop.
Again starting with a dry run with echo.

~~~
$ for infile in *.fq.gz
do
  outfile="$(basename $infile .fq.gz)"_qc.fq
  echo "trimmomatic SE -phred33 -threads 2 $infile trimmed/$outfile ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25"
done
~~~
{: .bash}


should be producing something like this

~~~
trimmomatic SE -phred33 -threads 2 Arabidopsis_sample1.fq.gz trimmed/Arabidopsis_sample1.fq.gz_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
trimmomatic SE -phred33 -threads 2 Arabidopsis_sample2.fq.gz trimmed/Arabidopsis_sample2.fq.gz_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
trimmomatic SE -phred33 -threads 2 Arabidopsis_sample3.fq.gz trimmed/Arabidopsis_sample3.fq.gz_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
trimmomatic SE -phred33 -threads 2 Arabidopsis_sample4.fq.gz trimmed/Arabidopsis_sample4.fq.gz_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
~~~
{: .output}


If it all looks ok, rerun with out `echo`

~~~
$ for infile in *.fq.gz
do
    outfile="$(basename $infile .fq.gz)"_qc.fq
    trimmomatic SE -phred33 -threads 2 $infile ../trimmed/"$outfile ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
done

~~~
{: .bash}


The following should appear:

~~~
TrimmomaticSE: Started with arguments:
 -phred33 -threads 2 Arabidopsis_sample1.fq.gz trimmed/Arabidopsis_sample1.fq.gz_qc.fq ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
Using Long Clipping Sequence: 'GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG'
Using Long Clipping Sequence: 'CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT'
…
Skipping duplicate Clipping Sequence: 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
ILLUMINACLIP: Using 0 prefix pairs, 31 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Reads: 250000 Surviving: 248320 (99.33%) Dropped: 1680 (0.67%)
TrimmomaticSE: Completed successfully
~~~
{: .output}

It's possible to scroll up to check if the percentage of surviving & dropped is within the same range in all of the samples.



# 3. Alignment to a reference genome

<img src="../img/RNAseqWorkflow.png" height="400" >

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. In this tutorial we will be using [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) but also 
a tool like [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml) does the job.



### STAR Alignment Strategy

STAR is shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive. The algorithm achieves this highly efficient mapping by performing a two-step process:

Seed searching
Clustering, stitching, and scoring
Seed searching

For every read that STAR aligns, STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the Maximal Mappable Prefixes (MMPs):

<img src="../img/alignment_STAR_step1.png" height="350" >


The different parts of the read that are mapped separately are called ‘seeds’. So the first MMP that is mapped to the genome is called seed1.

STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be seed2.

<img src="../img/alignment_STAR_step2.png" height="275" >

This sequential searching of only the unmapped portions of reads underlies the efficiency of the STAR algorithm. STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes. Other slower aligners use algorithms that often search for the entire read sequence before splitting reads and performing iterative rounds of mapping.

If STAR does not find an exact matching sequence for each part of the read due to mismatches or indels, the previous MMPs will be extended.

<img src="../img/alignment_STAR_step3.png" height="650" >

**If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

<img src="../img/alignment_STAR_step4.png" height="300" >

**Clustering, stitching, and scoring**

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).

<img src="../img/alignment_STAR_step5.png" height="400" >

### The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


### 3.1 Setting up

## 3.2 Index the reference genome
Our first step is to index the reference genome for use by STAR. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment (index files are not exchangeable between tools).

Take note that depending on the genome size these index files produced by STAR can be pretty big. Make sure there's enough disk space available.

~~~
$ cd /home/

$ mkdir genomeIndex

$ gunzip AtChromosome1.fa.gz

$ STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles AtChromosome1.fa --runThreadN 2
~~~
{: .bash}

While the index is created, you will see output something like this:

~~~
Apr 29 16:55:14 ..... Started STAR run
Apr 29 16:55:14 ... Starting to generate Genome files
Apr 29 16:55:16 ... starting to sort  Suffix Array. This may take a long time...
Apr 29 16:55:16 ... sorting Suffix Array chunks and saving them to disk...
Apr 29 16:56:13 ... loading chunks from disk, packing SA...
Apr 29 16:56:26 ... writing Suffix Array to disk ...
Apr 29 16:56:27 ... Finished generating suffix array
Apr 29 16:56:27 ... starting to generate Suffix Array index...
Apr 29 16:56:48 ... writing SAindex to disk
Apr 29 16:57:00 ..... Finished successfully

~~~
{: .output}


The indexing should have produced 8 star index files. Use the following command to see if they're really there. 

~~~
$ ls -l genomeIndex/
~~~
{: .bash}

result should be:
~~~
-rw-r--r-- 1 tbliek genseq-local          9 Apr 29 16:55 chrLength.txt
-rw-r--r-- 1 tbliek genseq-local         14 Apr 29 16:55 chrNameLength.txt
-rw-r--r-- 1 tbliek genseq-local          5 Apr 29 16:55 chrName.txt
-rw-r--r-- 1 tbliek genseq-local         11 Apr 29 16:55 chrStart.txt
-rw-r--r-- 1 tbliek genseq-local   30670848 Apr 29 16:55 Genome
-rw-r--r-- 1 tbliek genseq-local        290 Apr 29 16:55 genomeParameters.txt
-rw-r--r-- 1 tbliek genseq-local  249672325 Apr 29 16:56 SA
-rw-r--r-- 1 tbliek genseq-local 1565873616 Apr 29 16:56 SAindex
~~~
{: .output}



## 3.3 Align reads to reference genome

In some tools like hisat2 creating the sequence alignment files (bam-files) is done in two steps. first the aligning it self. After that the alignment file will be filtered for instance to only contain the reads that actualy map to the genome. This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed). STAR on the other hand has a build in filter and also a sort function. So the output is ready to use for downstream tools.  



First of course we will need to create a directory to output the alignment files

~~~
$ cd /home/

$ mkdir mapped
~~~
{: .bash}

Running STAR to align ( or map ) the reads and optionaly filter and sort them.

In contrast to most tools, STAR does not have a help function.
running STAR -h or STAR --help will result in an error. For information on what arguments to use you can 
use have a look at the 
[STAR manual.](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).


Here are some examples of common used arguments.

| argument   | meaning |
| ------- | ---------- |
| `--runThreads` | number of threads |
| `--genomeDir` | /path/to/genomeDir |
| `--readFilesIn` | /path/to/read1 [/path/to/read2] |
| `--readFilesCommand zcat` | when making use of gzipped fastq files |
| `--outFileNamePrefix` | /path/to/output file name |
| `--outSAMtype` | BAM/SAM or None  [optional: SortedByCoordinate] |
| `--outReadsUnmapped` | [default: None] Fastx ; output in separate fasta/fastq file |
| `--outFilterMultimapNmax` | [default: 10] max number of alignments accepted |
| `--outFilterMismatchNmax` | [default: 10] max number of mismatches accepted |
| `--outFilterMismatchNoverLmax`  | [default: 0.3] max fraction of mismatches mapped length |
| `--outFilterMismatchNoverReadLmax` | [default: 1.0] max fraction of mismatches read length |
| `--alignEndsType` | EndToEnd force end-to-end alignment, don't soft-clip |



For now we will be using STAR with the following arguments
~~~
$  STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/Arabidopsis_sample1_qc.fq --outFileNamePrefix mapped/Arabidopsis_sample1_qc --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All
~~~
{: .bash}


Next we want to make a loop to do all the files

It's good again to first start with a 'dry' run with the use of echo

~~~
$ for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)"_
   echo "STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All"
 done
~~~
{: .bash}

If the commands look good, rerun but this time without the echo.

~~~
$for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)"_
   STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outSAMattributes All
 done
~~~
{: .bash}

When running the STAR command, you will see output something like this:

~~~
May 04 12:51:55 ..... Started STAR run
May 04 12:52:47 ..... Started mapping
May 04 12:55:06 ..... Started sorting BAM
May 04 12:55:59 ..... Finished successfully
~~~
{: .output}

The final.out file contains all the characteristics of the alignment.

~~~
$ less mapped/Arabidopsis_sample1_qc.final.out
~~~
{: .bash}

resulting in table containing all the alignment values.

~~~
                                 Started job on |       May 04 12:51:55
                             Started mapping on |       May 04 12:52:47
                                    Finished on |       May 04 12:55:59
       Mapping speed, Million of reads per hour |       156.56

                          Number of input reads |       8349970
                      Average input read length |       101
                                    UNIQUE READS:
                   Uniquely mapped reads number |       8047613
                        Uniquely mapped reads % |       96.38%
                          Average mapped length |       100.60
                       Number of splices: Total |       2243517
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       2223720
                       Number of splices: GC/AG |       16101
                       Number of splices: AT/AC |       479
               Number of splices: Non-canonical |       3217
                      Mismatch rate per base, % |       0.24%
                         Deletion rate per base |       0.00%
                        Deletion average length |       1.30
                        Insertion rate per base |       0.01%
                       Insertion average length |       1.11
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       300104
             % of reads mapped to multiple loci |       3.59%
        Number of reads mapped to too many loci |       385
             % of reads mapped to too many loci |       0.00%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       0.02%
                     % of reads unmapped: other |       0.00%
Arabidopsis_sample1_qcLog.final.out (END) 
~~~
{: .output}


## 3.4 The SAM/BAM format
The [SAM file](https://github.com/adamfreedman/knowyourdata-genomics/blob/gh-pages/lessons/01-know_your_data.md#aligned-reads-sam),
is a tab-delimited text file that contains information for each individual read and its alignment to the genome. While we do not
have time to go in detail of the features of the SAM format, the paper by
[Heng Li et al.](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full) provides a lot more detail on the specification.

**The compressed binary version of SAM is called a BAM file.** We use this version to reduce size and to allow for *indexing*, which enables efficient random access of the data contained within the file.

The file begins with a **header**, which is optional. The header is used to describe source of data, reference sequence, method of
alignment, etc., this will change depending on the aligner being used. Following the header is the **alignment section**. Each line
that follows corresponds to alignment information for a single read. Each alignment line has **11 mandatory fields** for essential
mapping information and a variable number of other fields for aligner specific information. An example entry from a SAM file is
displayed below with the different fields highlighted.

<img src="../img/sam_bam_1.png">

<img src="../img/sam_bam2.png">


# 4. Creating the counts file

For downstream application for each of the samples the number of reads that maps within a gene has to be determined.
Featurecounts from the subread package can do this.

FeatureCounts can count the number of reads that map within a feature. In case of the arabidopsis annotation there are three different features to choose from. Depending on the downstream applications the choice is gene, transcript or exon. In this study we are just looking for differientially expressed genes.

~~~
$ cd /home/

$ gunzip ath_annotation.gff3.gz

$ featureCounts -O -t gene -g ID -a ath_annotation.gff3 -o counts.txt mapped/*.bam
~~~

-a <string>         Name of an annotation file. GTF/GFF format by default.
-o <string>         Name of the output file including read counts.
-O                  Assign reads to all their overlapping meta-features.
-t <string>         Specify feature type in GTF annotation
-g <string>         Specify attribute type in GTF annotation. Determines the name of the features.

The output file produced by `featureCounts` is a tab-delimited file, can be opened in a program like excel.



# 5. Removal of Container and Image

If you have run this lesson locally and finished it all you might want to remove the container and the image (occupies about 4 gb of space).

If there are file you want to save copy them using `docker cp` as mentioned earlier in this lesson.

First you can stop the container and remove it.

~~~
$ docker stop bioinfo

$ docker rm bioinfo
~~~
{: .bash}

To check if all the containers using the image are gone run

~~~
$ docker ps -a
~~~
{: .bash}

This will show a list of the conainers, and the image they make use of.

If empty for the image in question, run the following to get the ID of the image.

~~~
$ docker images
~~~
{: .bash}

Copy the image id from the output
~~~
REPOSITORY                         TAG                 IMAGE ID            CREATED             SIZE
scienceparkstudygroup/master-gls   fastq-latest        daef7efb73ec        4 days ago          1.36GB
~~~
{: .output}

To remove the image run

~~~
$ docker rmi daef7efb73ec
~~~
{: .bash}