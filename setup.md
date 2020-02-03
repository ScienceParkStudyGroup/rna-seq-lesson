---
title: 
---

# Setup

## Softwares

> ## What you need to install.
>
> Before the training, please make sure you have done the following: 
>
> 1. Download and install **up-to-date versions** of:
>    - R: [https://cloud.r-project.org](https://cloud.r-project.org).
>    - RStudio: [http://www.rstudio.com/download](http://www.rstudio.com/download). 
>    - The `DESeq2` package: [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
>    - The `tidyverse` package: [https://www.tidyverse.org/](https://www.tidyverse.org/).
> 2. Read the workshop [Code of Conduct](https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html) to make sure this workshop stays welcoming for everybody.
> 3. Get comfortable: if you're not in a physical workshop, be set up with two screens if possible. You will be following along in RStudio on your own computer while also following this tutorial on your own.
> More instructions are available on the workshop website in the **Setup** section.
{: .prereq}

## Lesson files 

> ## What you need to download.
>
> Please download the necessary data files for the lesson [from the Zenodo archive](https://zenodo.org/record/3627021).  
>
> - **Counts**: A `counts.txt` dataframe of the sample raw counts. It is a tab separated file therefore data are in tabulated separated columns.
> - **Experimental design**: the `Experimental_design.txt` dataframe indicates the correspondence between samples and experimental conditions (e.g. control, treated).   
> - Please read the original study description below and have a look at the file preview to understand their format.  
> - These two files were obtained by running the `v0.1.1` version of a [RNA-Seq bioinformatic pipeline](https://github.com/KoesGroup/Snakemake_hisat-DESeq/blob/master/README.md) on the [mRNA-Seq sequencing files from Vogel et al. (2016)](https://www.ebi.ac.uk/ena/data/view/PRJEB13938).
{: .prereq}

### Original study
This RNA-seq lesson will make use of a dataset from a study on the model plant _Arabidopsis thaliana_ inoculated with commensal leaf bacteria (_Methylobacterium extorquens_ or _Sphingomonas melonis_) and infected or not with a leaf bacterial pathogen called _Pseudomonas syringae_. Leaf samples were collected from Arabidopsis plantlets from plants inoculated or not with commensal bacteria and infected or not with the leaf pathogen either after two days (2 dpi, dpi: days post-inoculation) or seven days (6 dpi). 

All details from the study are available in [Vogel et al. in 2016 and was published in New Phytologist](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.14036).  


### File previews

#### Gene counts
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

#### Experimental design table

| sample     | growth condition                             | infect                      | sampling time point (days) |
|------------|----------------------------------------------|-----------------------------|----------------------------|
| ERS1147158 | 10 mM MgCl2                                  | mock                        | 2                          |
| ERS1147169 | 10 mM MgCl2                                  | mock                        | 2                          |
| ERS1147182 | 10 mM MgCl2                                  | mock                        | 2                          |
| ERS1147194 | 10 mM MgCl2                                  | mock                        | 2                          |
| ERS1147204 | 10 mM MgCl2                                  | mock                        | 7                          |
| ERS1147206 | 10 mM MgCl2                                  | mock                        | 7                          |
| ERS1147207 | 10 mM MgCl2                                  | mock                        | 7                          |
| ERS1147208 | 10 mM MgCl2                                  | mock                        | 7                          |
| ERS1147159 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 2                          |
| ERS1147160 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 2                          |
| ERS1147161 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 2                          |
| ERS1147209 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 2                          |
| ERS1147162 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 7                          |
| ERS1147163 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 7                          |
| ERS1147164 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 7                          |
| ERS1147165 | 10 mM MgCl2                                  | Pseudomonas syringae DC3000 | 7                          |
| ERS1147186 | Methylobacterium extorquens PA1 colonization | mock                        | 2                          |
| ERS1147187 | Methylobacterium extorquens PA1 colonization | mock                        | 2                          |
| ERS1147188 | Methylobacterium extorquens PA1 colonization | mock                        | 2                          |
| ERS1147189 | Methylobacterium extorquens PA1 colonization | mock                        | 2                          |
| ERS1147190 | Methylobacterium extorquens PA1 colonization | mock                        | 7                          |
| ERS1147191 | Methylobacterium extorquens PA1 colonization | mock                        | 7                          |
| ERS1147192 | Methylobacterium extorquens PA1 colonization | mock                        | 7                          |
| ERS1147195 | Methylobacterium extorquens PA1 colonization | mock                        | 7                          |
| ERS1147196 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 2                          |
| ERS1147197 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 2                          |
| ERS1147198 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 2                          |
| ERS1147199 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 2                          |
| ERS1147200 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 7                          |
| ERS1147201 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 7                          |
| ERS1147202 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 7                          |
| ERS1147203 | Methylobacterium extorquens PA1 colonization | Pseudomonas syringae DC3000 | 7                          |
| ERS1147166 | Sphingomonas melonis Fr1 colonization        | mock                        | 2                          |
| ERS1147167 | Sphingomonas melonis Fr1 colonization        | mock                        | 2                          |
| ERS1147168 | Sphingomonas melonis Fr1 colonization        | mock                        | 2                          |
| ERS1147171 | Sphingomonas melonis Fr1 colonization        | mock                        | 2                          |
| ERS1147172 | Sphingomonas melonis Fr1 colonization        | mock                        | 7                          |
| ERS1147173 | Sphingomonas melonis Fr1 colonization        | mock                        | 7                          |
| ERS1147174 | Sphingomonas melonis Fr1 colonization        | mock                        | 7                          |
| ERS1147175 | Sphingomonas melonis Fr1 colonization        | mock                        | 7                          |
| ERS1147176 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 2                          |
| ERS1147177 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 2                          |
| ERS1147178 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 2                          |
| ERS1147179 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 2                          |
| ERS1147181 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 7                          |
| ERS1147183 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 7                          |
| ERS1147184 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 7                          |
| ERS1147185 | Sphingomonas melonis Fr1 colonization        | Pseudomonas syringae DC3000 | 7                          |


