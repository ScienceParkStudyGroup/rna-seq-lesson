---
title: 
---


# Setup

<!-- MarkdownTOC autolink="True" levels="1,2" -->

- [1. Option 1 \(preferred\): using a Docker image](#1-option-1-preferred-using-a-docker-image)
	- [1.1 Installing Docker](#11-installing-docker)
	- [1.2 The `fastq-latest` image for bioinformatic steps \(episodes 03 and 04\)](#12-the-fastq-latest-image-for-bioinformatic-steps-episodes-03-and-04)
	- [1.3 The `rnaseq-latest` image for count data analysis \(episodes 05, 06 and 07\)](#13-the-rnaseq-latest-image-for-count-data-analysis-episodes-05-06-and-07)
- [2. Option 2: manual installation](#2-option-2-manual-installation)
	- [2.1 Softwares and packages](#21-softwares-and-packages)
	- [2.2 Data files](#22-data-files)
- [3. Original study](#3-original-study)
	- [3.1 Gene counts](#31-gene-counts)
	- [3.2 Experimental design table](#32-experimental-design-table)

<!-- /MarkdownTOC -->

# 1. Option 1 (preferred): using a Docker image

## 1.1 Installing Docker

The preferred option to install all softwares and packages is to use a tailor-made Docker image. See [this nice introduction to Docker here](https://aws.amazon.com/docker/).   

There are two Docker images necessary to complete this RNA-seq lesson:
1. The command-line Docker `fastq-latest` image necessary to perform all bioinformatic analyses on the sequencing files: trimming, alignment and count table generation.
2. The RStudio Docker `rnaseq-latest` image necessary to perform all count-related analyses: EDA, differential expression and downstream functional analyses.   


So first thing first, we need to install Docker itself. 

> ## Install Docker
> Unfortunately, in many common situations installing Docker on your laptop will not straightforward if you do not have a large amount of technical experience. We have helpers on hand that have worked their way through the install process but be prepared for some troubleshooting.
> Please try to install the appropriate software from the list below depending on the operating system that your laptop is running:
> ### Microsoft Windows
> **You must have admin rights to run docker!** Some parts of the lesson will work without running as admin but if you are unable to `Run as admin` on your machine some of this workshop might not work easily.
> 
> If you have Windows 10 Pro Edition:
>  - First try to install the [Docker Desktop (Windows)](https://hub.docker.com/editions/community/docker-ce-desktop-windows), or **failing that**;
> - Install the [Docker Toolbox (Windows)](https://docs.docker.com/toolbox/toolbox_install_windows/).
>
> If you have Windows 10 Home Edition:
> - Install the [Docker Toolbox (Windows)](https://docs.docker.com/toolbox/toolbox_install_windows/).
>
> ### Apple macOS
> Either:
> - First, try to install the [Docker Desktop (Mac)](https://hub.docker.com/editions/community/docker-ce-desktop-mac), or **failing that**:
> - Install the [Docker Toolbox (Mac)](https://docs.docker.com/toolbox/toolbox_install_mac/).
> 
> ### Linux
> There are too many varieties of Linux to give precise instructions here, but hopefully you can locate documentation for getting Docker installed on your Linux distribution. It may already be installed. Note that Docker do list a number of versions of the Docker Engine for different Linux distributions [here](https://hub.docker.com/search/?type=edition&offering=community). 
>
> ### Troubleshooting
> Sometimes with git-bash and Windows, you can get issues listed here:   
> `the input device is not a TTY.  If you are using mintty, try prefixing the command with 'winpty'`. This can be troubleshooted following [this blog post](https://pitman.io/posts/tips-for-using-docker-from-git-bash-on-windows/).
{: .prereq}



## 1.2 The `fastq-latest` image for bioinformatic steps (episodes 03 and 04)

This Docker image will allow you to complete the [episodes 03](https://scienceparkstudygroup.github.io/rna-seq-lesson/03-qc-of-sequencing-results/index.html) and [04](https://scienceparkstudygroup.github.io/rna-seq-lesson/04-bioinformatic-workflow/index.html) that work on `.fastq` sequencing files.

The Docker image is called `fastq-latest` and contains softwares and data required for the command-line part of the lesson. It can be found [found at the Science Park Study Group DockerHub](https://hub.docker.com/repository/docker/scienceparkstudygroup/master-gls) with the tag `fastq-latest`.

> ## Before you start
>
> Before the training, please make sure you have done the following: 
>
> 1. First, install [Docker desktop](https://www.docker.com/products/docker-desktop) for your operating system (Mac OS X or Windows).  
> 2. If needed, install Shell Bash: [follow these instructions](http://swcarpentry.github.io/shell-novice/setup.html).
> 3. Open a new Shell Bash window and navigate to a folder that will be your workspace. For instance, you could create a folder named `rnaseq-tutorial/` on your Desktop and move inside with the Shell using `cd ~/Desktop/rnaseq-tutorial/`. 
> 4. In a Shell Bash window, type the following command: `docker run -it --name bioinfo -v $PWD:/home/ scienceparkstudygroup/master-gls:fastq-latest`. This will download a Docker image for the bioinformatic part of the course, create and run a container where Bash will be running. You will enter the container directly where you can start working.     
> 5. To quit, type `exit` and you will exit the container and be on your machine file system again. The container will be stopped. 
> 6. To go back to the container, type `docker start bioinfo` and then `docker exec -it bioinfo bash`. You will enter inside the container again where you can find all softwares and data. 
> 7. Type `exit` to go back to your file system. 
{: .prereq}

__Docker command-line explanations:__  
- The `--it` starts an interactive session in which you directly start AND enter the container.       
- The `--name` gives a name to the container for easy retrieval.  
- The `-v $PWD:/home/` maps your working directory (e.g. `~/Desktop/rnaseq-tutorial`) to the container `/home/` folder. 


## 1.3 The `rnaseq-latest` image for count data analysis (episodes 05, 06 and 07)

This image is based on a [Bioconductor Docker imag release 3.10](https://hub.docker.com/r/bioconductor/bioconductor_docker/tags) image with additional packages such as `pheatmap` or `tidyverse`.

The latest image can be [found at the Science Park Study Group DockerHub](https://hub.docker.com/repository/docker/scienceparkstudygroup/master-gls) with the tag `rnaseq-latest`.


> ## Before you start
>
> Before the training, please make sure you have done the following: 
>
> 1. First, install [Docker desktop](https://www.docker.com/products/docker-desktop) for your operating system.  
> 2. If needed, install Shell Bash: [follow these instructions](http://swcarpentry.github.io/shell-novice/setup.html).
> 3. Open a new Shell Bash window and navigate to a folder that will be your workspace. For instance, you could create a folder named `rnaseq-tutorial/` on your Desktop and move inside with the Shell using `cd ~/Desktop/rnaseq-tutorial/`. 
> 4. In a Shell Bash window, type the following command: `docker run --rm --name rstudio_instance -v $PWD:/home/rstudio/ -e PASSWORD=mypwd -p 8787:8787 scienceparkstudygroup/master-gls:rnaseq-latest`. This will download a Docker image for the course, create and run a container where RStudio will be running.   
> 4. Navigate to [http://localhost:8787](http://localhost:8787) in your web browser. You should have an RStudio session running. Type `rstudio` as the user name and `mypwd` as your password. 
> 5. To quit, close the web browser window where RStudio is running and exit the Shell too. 
{: .prereq}



> ## Important note
>
> You can save files to your disk when working inside the Docker-powered R session. You need to save them as you would normally. The files (e.g. `my_plot.png`) will be where you were working (the directory from which you launched the Docker container). 
>
{: .callout}


__Docker command-line explanations:__  
- The `--rm` removes the container when it has been run. No need to store it into your computer after use.      
- The `--name` gives a name to the running container for easy retrieval.  
- The `-p 8787:8787` follow the format `-p host_port:container_port`. Therefore the port 8787 inside the container will be exposed to the outside port on the host machine. That way, the running instance of RStudio can be access through the <IP address>:port format.

# 2. Option 2: manual installation
This is the second way to install softwares and packages. It _should_ work but there is no guarantee that it _will_ work since R and packages versions on your machine might be different from the software and package versions used in this lesson. Thus, the preferred way is still to use the Docker image (option 1).  

## 2.1 Softwares and packages

> ## Before you start.
>
> Before the training, please make sure you have done the following: 
>
> 1. Download and install **up-to-date versions** of:
>    - R: [https://cloud.r-project.org](https://cloud.r-project.org).
>    - RStudio: [http://www.rstudio.com/download](http://www.rstudio.com/download). 
>    - The `DESeq2` package: [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
>    - The `tidyverse` package: [https://www.tidyverse.org/](https://www.tidyverse.org/).
>    - The `EnhancedVolcano` package: [http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html](http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
>    - The `pheatmap` package: [https://cran.r-project.org/web/packages/pheatmap/index.html](https://cran.r-project.org/web/packages/pheatmap/index.html).
>    - The `biomartr` package: [https://cran.r-project.org/web/packages/biomartr/index.html](https://cran.r-project.org/web/packages/biomartr/index.html).
>    - The `clusterProfiler` package: [https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).
>    - The `org.At.tair.db` package: [https://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html).
>    - The `biomaRt` package: [https://bioconductor.org/packages/release/bioc/html/biomaRt.html](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).
>    - The `pwr` package: [https://cran.r-project.org/web/packages/pwr/index.html](https://cran.r-project.org/web/packages/pwr/index.html).
> 2. Read the workshop [Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html) to make sure this workshop stays welcoming for everybody.
> 3. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also following this tutorial on your own.
> More instructions are available on the workshop website in the **Setup** section.
{: .prereq}


## 2.2 Data files 
https://zenodo.org/record/4530667

> ## What you need to download.
> [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4530667.svg)](https://doi.org/10.5281/zenodo.4530667)
> Please download the necessary data files for the lesson [from the Zenodo archive](https://doi.org/10.5281/zenodo.4530667).  
>
> - **Counts**: A `counts.txt` dataframe of the sample raw counts. It is a tab separated file therefore data are in tabulated separated columns.
> - **Experimental design**: the `experimental_design_modified.txt` dataframe indicates the correspondence between samples and experimental conditions (e.g. control, treated).  
- **Differentially expressed genes**: `differential_genes.tsv` dataframe contains the result of the DESeq2 analysis.  
> - Please read the original study description below and have a look at the file preview to understand their format.  
> - These `counts.txt` file was obtained by running the `v0.1.1` version of a [RNA-Seq bioinformatic pipeline](https://github.com/KoesGroup/Snakemake_hisat-DESeq/blob/master/README.md) on the [mRNA-Seq sequencing files from Vogel et al. (2016)](https://www.ebi.ac.uk/ena/data/view/PRJEB13938).
{: .prereq}

# 3. Original study
This RNA-seq lesson will make use of a dataset from a study on the model plant _Arabidopsis thaliana_ inoculated with commensal leaf bacteria (_Methylobacterium extorquens_ or _Sphingomonas melonis_) and infected or not with a leaf bacterial pathogen called _Pseudomonas syringae_. Leaf samples were collected from Arabidopsis plantlets from plants inoculated or not with commensal bacteria and infected or not with the leaf pathogen either after two days (2 dpi, dpi: days post-inoculation) or seven days (6 dpi). 

All details from the study are available in [Vogel et al. in 2016 and was published in New Phytologist](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14036).  


## 3.1 Gene counts
The dimension of this table are 33,769 rows x 49 columns.  
  * 33,769 rows: one for gene and sample names and the rest for gene counts.  
  * 49 columns: one for the gene id and the rest for sample accession identifiers (from the EBI European Nucleotide Archive).

| Geneid     | ERR1406259 | ERR1406260 | ERR1406261 | ERR1406262 | ERR1406263 | ERR1406264 | ERR1406265 | ERR1406266 | ERR1406268 | ERR1406269 | ERR1406270 | ERR1406271 | ERR1406272 | ERR1406273 | ERR1406274 | ERR1406275 | ERR1406276 | ERR1406277 | ERR1406278 | ERR1406279 | ERR1406280 | ERR1406281 | ERR1406282 | ERR1406284 | ERR1406285 | ERR1406286 | ERR1406287 | ERR1406288 | ERR1406289 | ERR1406290 | ERR1406291 | ERR1406292 | ERR1406293 | ERR1406294 | ERR1406296 | ERR1406297 | ERR1406298 | ERR1406299 | ERR1406300 | ERR1406301 | ERR1406302 | ERR1406303 | ERR1406304 | ERR1406305 | ERR1406306 | ERR1406307 | ERR1406308 | ERR1406309 |
|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|
| AT1G01010  | 59         | 81         | 40         | 51         | 57         | 110        | 93         | 87         | 99         | 131        | 80         | 79         | 142        | 216        | 102        | 76         | 92         | 116        | 100        | 126        | 151        | 249        | 61         | 189        | 161        | 92         | 80         | 125        | 77         | 106        | 90         | 86         | 164        | 71         | 64         | 83         | 100        | 86         | 91         | 214        | 142        | 76         | 84         | 123        | 91         | 69         | 75         | 85         |
| AT1G01020  | 365        | 466        | 440        | 424        | 393        | 567        | 397        | 468        | 465        | 365        | 382        | 365        | 595        | 509        | 323        | 422        | 325        | 358        | 415        | 403        | 498        | 501        | 441        | 498        | 409        | 396        | 472        | 566        | 422        | 462        | 504        | 434        | 717        | 534        | 408        | 346        | 757        | 456        | 443        | 976        | 517        | 467        | 533        | 648        | 457        | 393        | 538        | 579        |
| AT1G03987  | 8          | 16         | 13         | 19         | 13         | 20         | 19         | 24         | 8          | 10         | 10         | 14         | 11         | 13         | 10         | 9          | 11         | 20         | 14         | 10         | 10         | 8          | 14         | 25         | 14         | 13         | 18         | 17         | 19         | 4          | 12         | 14         | 29         | 15         | 19         | 47         | 28         | 6          | 21         | 20         | 5          | 5          | 8          | 17         |            |            |            |            |
| AT1G01030  | 111        | 200        | 189        | 164        | 141        | 389        | 200        | 175        | 127        | 186        | 140        | 189        | 147        | 193        | 102        | 101        | 103        | 128        | 136        | 120        | 162        | 229        | 124        | 177        | 125        | 136        | 169        | 197        | 141        | 217        | 214        | 180        | 253        | 161        | 98         | 152        | 371        | 219        | 170        | 566        | 441        | 99         | 207        | 220        | 169        | 117        | 123        | 183        |
| AT1G03993  | 131        | 179        | 169        | 157        | 114        | 156        | 138        | 184        | 193        | 143        | 135        | 155        | 218        | 236        | 159        | 194        | 149        | 156        | 168        | 128        | 174        | 269        | 183        | 215        | 176        | 165        | 171        | 247        | 179        | 181        | 177        | 199        | 313        | 236        | 154        | 169        | 313        | 201        | 202        | 332        | 169        | 218        | 203        | 250        | 190        | 188        | 223        | 218        |
| AT1G01040  | 1491       | 1617       | 1418       | 1543       | 1224       | 1635       | 1524       | 1665       | 1565       | 1566       | 1496       | 1499       | 2244       | 1881       | 1177       | 1751       | 1444       | 1631       | 1393       | 1407       | 1880       | 2311       | 1529       | 1919       | 1662       | 1537       | 1691       | 2142       | 1469       | 1733       | 1910       | 1873       | 3079       | 2179       | 1486       | 1471       | 2840       | 1891       | 1924       | 3136       | 1520       | 1901       | 1950       | 2596       | 1802       | 1851       | 2133       | 1984       |
| AT1G01046  | 35         | 30         | 48         | 32         | 28         | 50         | 51         | 56         | 36         | 26         | 29         | 38         | 48         | 30         | 15         | 44         | 23         | 31         | 22         | 27         | 33         | 51         | 41         | 35         | 48         | 38         | 41         | 49         | 27         | 36         | 39         | 50         | 57         | 49         | 41         | 30         | 54         | 41         | 43         | 85         | 42         | 42         | 59         | 65         | 49         | 64         | 50         | 46         |
| ath-miR838 | 12         | 11         | 22         | 18         | 15         | 21         | 22         | 24         | 16         | 12         | 10         | 15         | 17         | 16         | 7          | 20         | 11         | 14         | 6          | 11         | 16         | 17         | 17         | 15         | 26         | 12         | 17         | 13         | 15         | 12         | 18         | 25         | 26         | 25         | 15         | 15         | 22         | 20         | 14         | 37         | 20         | 20         | 22         | 27         | 17         | 21         | 23         | 23         |
| AT1G01050  | 1484       | 1483       | 1237       | 1544       | 1119       | 1453       | 1280       | 1256       | 1768       | 1869       | 1709       | 1649       | 2431       | 1858       | 1195       | 1518       | 1325       | 2013       | 1645       | 1666       | 2056       | 2258       | 1530       | 1834       | 1477       | 1532       | 1609       | 2220       | 1552       | 1976       | 1706       | 1807       | 2656       | 1873       | 1329       | 1512       | 2915       | 1646       | 1983       | 2687       | 1548       | 1740       | 1632       | 2330       | 1578       | 1521       | 1970       | 1977       |

... many more lines ...

## 3.2 Experimental design table

dpi: days post-inoculation. 

| sample     | growth                          | infected                    | dpi |
|------------|---------------------------------|-----------------------------|-----|
| ERR1406259 | MgCl2                           | mock                        | 2   |
| ERR1406271 | MgCl2                           | mock                        | 2   |
| ERR1406282 | MgCl2                           | mock                        | 2   |
| ERR1406294 | MgCl2                           | mock                        | 2   |
| ERR1406305 | MgCl2                           | mock                        | 7   |
| ERR1406306 | MgCl2                           | mock                        | 7   |
| ERR1406307 | MgCl2                           | mock                        | 7   |
| ERR1406308 | MgCl2                           | mock                        | 7   |
| ERR1406260 | MgCl2                           | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406261 | MgCl2                           | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406262 | MgCl2                           | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406309 | MgCl2                           | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406263 | MgCl2                           | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406264 | MgCl2                           | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406265 | MgCl2                           | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406266 | MgCl2                           | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406287 | Methylobacterium_extorquens_PA1 | mock                        | 2   |
| ERR1406288 | Methylobacterium_extorquens_PA1 | mock                        | 2   |
| ERR1406289 | Methylobacterium_extorquens_PA1 | mock                        | 2   |
| ERR1406290 | Methylobacterium_extorquens_PA1 | mock                        | 2   |
| ERR1406291 | Methylobacterium_extorquens_PA1 | mock                        | 7   |
| ERR1406292 | Methylobacterium_extorquens_PA1 | mock                        | 7   |
| ERR1406293 | Methylobacterium_extorquens_PA1 | mock                        | 7   |
| ERR1406296 | Methylobacterium_extorquens_PA1 | mock                        | 7   |
| ERR1406297 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406298 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406299 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406300 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406301 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406302 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406303 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406304 | Methylobacterium_extorquens_PA1 | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406268 | Sphingomonas_melonis_Fr1        | mock                        | 2   |
| ERR1406269 | Sphingomonas_melonis_Fr1        | mock                        | 2   |
| ERR1406270 | Sphingomonas_melonis_Fr1        | mock                        | 2   |
| ERR1406272 | Sphingomonas_melonis_Fr1        | mock                        | 2   |
| ERR1406273 | Sphingomonas_melonis_Fr1        | mock                        | 7   |
| ERR1406274 | Sphingomonas_melonis_Fr1        | mock                        | 7   |
| ERR1406275 | Sphingomonas_melonis_Fr1        | mock                        | 7   |
| ERR1406276 | Sphingomonas_melonis_Fr1        | mock                        | 7   |
| ERR1406277 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406278 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406279 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406280 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 2   |
| ERR1406281 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406284 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406285 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 7   |
| ERR1406286 | Sphingomonas_melonis_Fr1        | Pseudomonas_syringae_DC3000 | 7   |
