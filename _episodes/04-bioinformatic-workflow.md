---
title: "From fastq files to read counts"
teaching: 45
exercises: 0 
questions:
- ""
objectives:
- ""
keypoints:
- ""
---
# Quality control and trimming

<img src="https://github.com/ScienceParkStudyGroup/2019-03-07-rnaseq-workshop/blob/gh-pages/images/RNAseqWorkflow.png" height="400" >



## Running FastQC  

We will now assess the quality of the reads that we downloaded. First, we need to make an output directory for the fastqc results to be stored. This we want to do in the 'RNAseq070319' directory.

~~~
$ cd ~/RNAseq070319

$ mkdir fastqc
~~~


Next we need to get to the directory thay actually contains the the fastq files.

~~~
$ cd ~/RNAseq070319/rawReads
~~~



Running fastqc uses the following command

~~~
fastqc -o ../fastqc $filename
~~~


Of course we don't want to do y=this for all the samples seperately so we can loop through the list of samples and run them all sequentially

With the use of echo you can start off with a "dry run"

~~~
$ for filename in *.fastq
  do
    echo fastqc -o ../fastqc $filename
  done
~~~


The echo command only prints the commands to the screen, and doesn't really run it.

~~~
fastqc -o ../fastqc sub06.fastq
fastqc -o ../fastqc sub07.fastq
fastqc -o ../fastqc sub08.fastq
fastqc -o ../fastqc sub21.fastq
fastqc -o ../fastqc sub23.fastq
fastqc -o ../fastqc sub24.fastq
~~~


If it looks good remove the echo and go for it.

~~~
$ for filename in *.fastq
  do
    fastqc -o ../fastqc $filename
  done
~~~


You will see an automatically updating output message telling you the
progress of the analysis. It shoul look something like this:

~~~
Started analysis of sub06.fastq
Approx 5% complete for sub06.fastq
Approx 10% complete for sub06.fastq
Approx 15% complete for sub06.fastq
Approx 20% complete for sub06.fastq
…
Approx 85% complete for sub24.fastq
Approx 90% complete for sub24.fastq
Approx 95% complete for sub24.fastq
Approx 100% complete for sub24.fastq
Analysis complete for sub24.fastq
~~~


In total, it should take about five minutes for FastQC to run on all
six of our FASTQ files.

If the command doesn't run or you want more information on fastqc, run the following to get the help page.

~~~
$ fastqc -h
~~~

But if all went right, the FastQC program will have created several new files within our
`~/RNAseq070319/fastqc` directory.

~~~
$ cd ~/RNAseq070319/fastqc
$ ls
~~~


~~~
sub06_fastqc.html  sub07_fastqc.zip   sub21_fastqc.html  sub23_fastqc.zip
sub06_fastqc.zip   sub08_fastqc.html  sub21_fastqc.zip   sub24_fastqc.html
sub07_fastqc.html  sub08_fastqc.zip   sub23_fastqc.html  sub24_fastqc.zip
~~~



## Viewing the FastQC results

If we were working on our local computers, we'd be able to display each of these
HTML files as a webpage:

~~~
$ open sub06_fastqc.html
~~~


However, if you try this on our genseq instance, you'll get an error:

~~~
Couldn't get a file descriptor referring to the console
~~~


This is because the genseq instance we're using doesn't have any web
browsers installed on it, so the remote computer doesn't know how to
open the file. We want to look at the webpage summary reports, so
let's transfer them to our local computers (i.e. your laptop).

To transfer a file from a remote server to our own machines, we will
use `scp`.

First we
will make a new directory on our computer to store the HTML files
we're transfering. Let's put it on our desktop for now. Open a new
tab in your terminal program (you can use the pull down menu at the
top of your screen or the Cmd+t keyboard shortcut) and type:

~~~
$ mkdir -p ~/Desktop/fastqc_html
~~~


Now we can transfer our HTML files to our local computer using `scp`.

~~~
$ scp tbliek@genseq-cn02.science.uva.nl:~/RNAseq070319/fastqc/*.html ~/Desktop/fastqc_html
~~~


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

You should see a status output like this:

~~~
sub06_fastqc.html                      100%  249KB 152.3KB/s   00:01    
sub06_fastqc.html                      100%  254KB 219.8KB/s   00:01    
sub06_fastqc.html                      100%  254KB 271.8KB/s   00:00    
sub06_fastqc.html                      100%  251KB 252.8KB/s   00:00    
sub06_fastqc.html                      100%  249KB 370.1KB/s   00:00    
sub06_fastqc.html                      100%  251KB 592.2KB/s   00:00  
~~~


Now we can go to our new directory and open the HTML files.

~~~
$ cd ~/Desktop/fastqc_html/
$ open *.html
~~~


Your computer will open each of the HTML files in your default web
browser. Depending on your settings, this might be as six separate
tabs in a single window or six separate browser windows.

 
## Decoding the other FastQC outputs
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

## Working with the FastQC text output

Now that we've looked at our HTML reports to get a feel for the data,
let's look more closely at the other output files. Go back to the tab
in your terminal program that is connected to your AWS instance
(the tab lab will start with `dcuser@ip`) and make sure you're in
our results subdirectory.   

~~~
$ cd ~/RNAseq070319/fastqc/
$ ls
~~~


~~~
sub06_fastqc.html  sub07_fastqc.zip   sub21_fastqc.html  sub23_fastqc.zip
sub06_fastqc.zip   sub08_fastqc.html  sub21_fastqc.zip   sub24_fastqc.html
sub07_fastqc.html  sub08_fastqc.zip   sub23_fastqc.html  sub24_fastqc.zip
~~~


Our `.zip` files are compressed files. They each contain multiple
different types of output files for a single input FASTQ file. To
view the contents of a `.zip` file, we can use the program `unzip`
to decompress these files. Let's try doing them all at once using a
wildcard.

~~~
$ unzip *.zip
~~~


~~~
Archive:  SRR2584863_1_fastqc.zip
caution: filename not matched:  sub07_fastqc.zip
caution: filename not matched:  sub08_fastqc.zip
caution: filename not matched:  sub21_fastqc.zip
caution: filename not matched:  sub22_fastqc.zip
caution: filename not matched:  sub24_fastqc.zip
~~~


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
> unzip $filename
> done
~~~


In this example, the input is six filenames (one filename for each of our `.zip` files).
Each time the loop iterates, it will assign a file name to the variable `filename`
and run the `unzip` command.
The first time through the loop,
`$filename` is `sub06_fastqc.zip`.
The interpreter runs the command `unzip` on `sub06_fastqc.zip`.
For the second iteration, `$filename` becomes
`Ssub07_fastqc.zip`. This time, the shell runs `unzip` on `sun07_fastqc.zip`.
It then repeats this process for the four other `.zip` files in our directory.


When we run our `for` loop, you will see output that starts like this:

~~~
Archive:  sub06_fastqc.zip
   creating: sub06_fastqc/
   creating: sub06_fastqc/Icons/
   creating: sub06_fastqc/Images/
  inflating: sub06_fastqc/Icons/fastqc_icon.png  
  inflating: sub06_fastqc/Icons/warning.png  
  inflating: sub06_fastqc/Icons/error.png  
  inflating: sub06_fastqc/Icons/tick.png  
  inflating: sub06_fastqc/summary.txt  
  inflating: sub06_fastqc/Images/per_base_quality.png  
  inflating: sub06_fastqc/Images/per_tile_quality.png  
  inflating: sub06_fastqc/Images/per_sequence_quality.png  
  inflating: sub06_fastqc/Images/per_base_sequence_content.png  
  inflating: sub06_fastqc/Images/per_sequence_gc_content.png  
  inflating: sub06_fastqc/Images/per_base_n_content.png  
  inflating: sub06_fastqc/Images/sequence_length_distribution.png  
  inflating: sub06_fastqc/Images/duplication_levels.png  
  inflating: sub06_fastqc/Images/adapter_content.png  
  inflating: sub06_fastqc/fastqc_report.html  
  inflating: sub06_fastqc/fastqc_data.txt  
  inflating: sub06_fastqc/fastqc.fo  
~~~


The `unzip` program is decompressing the `.zip` files and creating
a new directory (with subdirectories) for each of our samples, to
store all of the different output that is produced by FastQC. There
are a lot of files here. The one we're going to focus on is the
`summary.txt` file.

If you list the files in our directory now you will see:

~~~
sub06_fastqc       sub08_fastqc       sub22_fastqc
sub06_fastqc.html  sub08_fastqc.html  sub22_fastqc.html
sub06_fastqc.zip   sub08_fastqc.zip   sub22_fastqc.zip
sub07_fastqc       sub21_fastqc       sub24_fastqc
sub07_fastqc.html  sub21_fastqc.html  sub24_fastqc.html
sub07_fastqc.zip   sub21_fastqc.zip   sub24_fastqc.zip
~~~


The `.html` files and the uncompressed `.zip` files are still present,
but now we also have a new directory for each of our samples. We can
see for sure that it's a directory if we use the `-F` flag for `ls`.

~~~
$ ls -F
~~~


~~~
sub06_fastqc/      sub08_fastqc/      sub22_fastqc/
sub06_fastqc.html  sub08_fastqc.html  sub22_fastqc.html
sub06_fastqc.zip   sub08_fastqc.zip   sub22_fastqc.zip
sub07_fastqc/      sub21_fastqc/      sub24_fastqc/
sub07_fastqc.html  sub21_fastqc.html  sub24_fastqc.html
sub07_fastqc.zip   sub21_fastqc.zip   sub24_fastqc.zip
~~~


Let's see what files are present within one of these output directories.

~~~
$ ls -F sub06_fastqc/
~~~


~~~
fastqc_data.txt  fastqc.fo  fastqc_report.html	Icons/	Images/  summary.txt
~~~


Use `less` to preview the `summary.txt` file for this sample.

~~~
$ less sub06_fastqc/summary.txt
~~~


~~~
PASS    Basic Statistics        sub06.fastq
PASS    Per base sequence quality       sub06.fastq
PASS    Per tile sequence quality       sub06.fastq
PASS    Per sequence quality scores     sub06.fastq
WARN    Per base sequence content       sub06.fastq
WARN    Per sequence GC content sub06.fastq
PASS    Per base N content      sub06.fastq
PASS    Sequence Length Distribution    sub06.fastq
PASS    Sequence Duplication Levels     sub06.fastq
PASS    Overrepresented sequences       sub06.fastq
WARN    Adapter Content sub06.fastq
~~~


The summary file gives us a list of tests that FastQC ran, and tells
us whether this sample passed, failed, or is borderline (`WARN`). Remember to quit from `less` you enter `q`.

## Documenting Our Work

We can make a record of the results we obtained for all our samples
by concatenating all of our `summary.txt` files into a single file
using the `cat` command. We'll call this `full_report.txt` and move
it to `~/dc_workshop/docs`.

~~~
$ cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
~~~


# Trimming and filtering

Before we will do the alignment we need to remove sequences of low quality and sequences that are to short (below 25 bases).
Also in this case we will trim down long sequences to 100 bases, quality of the Ion-torrent reads drops the further it gets.
When making use of illumina reads this is not as much of a problem and 3'-trimming would then be a waste of data.

To start off make a directory trimmed for the output and then back to the rawReads directory.

~~~
$ cd ~/RNAseq070319/
$ mkdir trimmed
$ cd ~/RNAseq070319/rawReads/
~~~




The trimming and quality filtering will be done with trimmomatic.
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
trimmomatic SE -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
~~~



Of cource we don't want to do this for all the reads seperately so lets create a loop through all the fastq files.

When doing the fastqc only input files needed to be specified. In this case both the input and a matching output filenames need to be given.
this can be done with the help of 'basename'


~~~
$ for infile in *.fastq
do
 echo inputfile $infile
 outfile="$(basename $infile .fastq)"_qc.fq
 echo outputfile $outfile
 echo
done
~~~


This be be producing the following list

~~~
inputfile sub06.fastq
outputfile sub06_qc.fq

inputfile sub07.fastq
outputfile sub07_qc.fq

inputfile sub08.fastq
outputfile sub08_qc.fq

inputfile sub21.fastq
outputfile sub21_qc.fq

inputfile sub23.fastq
outputfile sub23_qc.fq

inputfile sub24.fastq
outputfile sub24_qc.fq
~~~


Next we can start writing the trimmomatic loop.
Again starting with a dry run with echo.

~~~
$ for fastq in *.fastq
do
  outputFile="$(basename $fastq .fastq)"_qc.fq
  echo trimmomatic SE -phred33 -threads 2 $fastq ../trimmed/$outputFile ILLUMINACLIP:../adapters.fasta:2:30:10 LEADING:3 TRAILING:3     SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
done
~~~


should be producing something like this

~~~
trimmomatic SE -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub07.fastq ../trimmed/sub07_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub08.fastq ../trimmed/sub08_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub21.fastq ../trimmed/sub21_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub23.fastq ../trimmed/sub23_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
trimmomatic SE -phred33 -threads 2 sub24.fastq ../trimmed/sub24_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
~~~


If it al seems ok rerun with out 'echo'

~~~
$ for infile in *.fastq
do
    outfile="$(basename $infile .fastq)"_qc.fq
    trimmomatic SE -phred33 -threads 2 $infile ../trimmed/"$outfile ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
done
~~~


The following should appear:

~~~
TrimmomaticSE: Started with arguments:
 -phred33 -threads 2 sub06.fastq ../trimmed/sub06_qc.fq ILLUMINACLIP:../general/adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 CROP:100
…
Input Reads: 1000000 Surviving: 887553 (88.76%) Dropped: 112447 (11.24%)
TrimmomaticSE: Completed successfully
~~~


It's possible to scroll up to check if the percentage of surviving & dropped is within the same range in all of the samples.



# Alignment to a reference genome

<img src="../images/RNAseqWorkflow.png" height="400" >

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to
choose from and, while there is no gold standard, there are some tools that are better suited for particular NGS analyses. In this tutorial we will be using [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) but also 
a tool like [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml) does the job.

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome


# Setting up



### Index the reference genome
Our first step is to index the reference genome for use by STAR. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment (index files are not exchangeable between tools).

Take note that depending on the genome size these index files produced by STAR can be pretty big. Make sure there's enough disk space available.

~~~
$ cd ~/RNAseq070319/general

$ mkdir genomeIndex

$ STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles AtChromosome1.fa --runThreadN 8

~~~


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



The indexing should have produced 8 star index files. Use the following command to see if they're really there. 

~~~
$ ls -l genomeIndex/
~~~


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




### Align reads to reference genome

In some tools like hisat2 creating the sequence alignment files (bam-files) is done in two steps. first the aligning it self. After that the alignment file will be filtered for instance to only contain the reads that actualy map to the genome. This is done with [sam flags](https://broadinstitute.github.io/picard/explain-flags.html) in samtools view (with the '-F 4' all the unmapped reads will be removed). STAR on the other hand has a build in filter and also a sort function. So the output is ready to use for downstream tools.  



First of course we will need to create a directory to output the alignment files

~~~
$ cd ~/RNAseq070319/

$ mkdir mapped
~~~


Running STAR to align ( or map ) the reads and optionaly filter and sort them.

In contrast to most atools, STAR does not have a help function.
running STAR -h or STAR --help will result in an error. For information on what arguments to use you can or need to 
use have a look at the 
[STAR manual.](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).


Here are some examples of common used arguments.

| argument   | meaning |
| ------- | ---------- |
| `--runThreads` | number of threads |
| `--genomeDir` | /path/to/genomeDir |
| `--readFilesIn` | /path/to/read1 [/path/to/read2] |
| `--readFilesCommand zcat` | when make use of gzipped fastq files |
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

$  STAR --genomeDir genomeindex --runThreadN 2 --readFilesIn ERR1406259.fq.gz --readFilesCommand zcat --outFileNamePrefix ERR1406259 --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --outFilterMultimapNmax 1

~~~



Next we want to make a loop to do all the files

It's good again to first start with a 'dry' run with the use of echo

~~~

$ for infile in trimmed/*.fq
 do
   outfile="$(basename $infile .fq)”
   echo "STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --readFilesCommand zcat --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --alignEndsType EndToEnd --outFilterMultimapNmax 1"
 done
~~~


If the commands look good, rerun but this time without the echo.

~~~
$for infile in *.fq
 do
   outfile="$(basename $infile .fq)”
   STAR --genomeDir genomeIndex --runThreadN 2 --readFilesIn trimmed/$infile --readFilesCommand zcat --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outFilterMismatchNmax 3 --alignEndsType EndToEnd --outFilterMultimapNmax 1
 done
~~~


When running the STAR command, you will see output something like this:

~~~
May 04 12:51:55 ..... Started STAR run
May 04 12:52:47 ..... Started mapping
May 04 12:55:06 ..... Started sorting BAM
May 04 12:55:59 ..... Finished successfully
~~~


The final.put file contains all the characteristics of the alignment.

~~~
$ less ERR1406259Log.final.out
~~~

resulting in table containing all the alignment numbers.

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
ERR1406259Log.final.out (END) 
~~~




#### SAM/BAM format
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

![sam_bam_1](../images/sam_bam_1.png)

![sam_bam2](../images/sam_bam2.png)


### Creating the counts file

For downstream application for each of the samples the number of reads that maps within a gene has to be determent.
Featurecounts from the subread package can do this.


~~~
$ cd ~/RNAseqWorkshop/mapped

$ featureCounts -O -t mRNA -g ID -a ../general/annotation.all_transcripts.exon_features.ath.gff3 -o counts.txt *.bam
~~~



The file produced by featureCounts is a tab-delimited file.

