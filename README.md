-----

# Introduction

This R script is to demonstrate the steps to download GSE data from NCBI
GEO database, and to obtain the differential expressed genes from GSE
data.

## A little background of GEO

The Gene Expression Omnibus (GEO) is a data repository hosted by the
National Center for Biotechnology Information (NCBI). NCBI contains all
publicly available nucleotide and protein sequences. Presently, all
records in GenBank NCBI are generated from direct submission to the DNA
sequence databases from the original authors, who volunteer their
records to make the data publicly available or do so as part of the
publication process. The NCBI GEO is intended to house different types
of expression data, covering all type of sequencing data in both raw and
processed formats.

## Example here

Here, GSE63477, which is a microarray data, will be analysed. It
contains an expression profile of prostate cancer cells (LNCaP) after
treatment of cabazitaxel or docetaxel for 16 hr. You may read the
details in NCBI GEO, under the “overall design” section. Or for better
understanding, it’s always good if we read the original paper…link here:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi>

-----

# Using GEOquery to obtain microarray data

First, set the working directory.

The GEOquery package allows you to access the data from GEO. Depending
on your needs, you can download only the processed data and metadata
provided by the depositor. In some cases, you may want to download the
raw data as well, if it was provided by the depositor.

The function to download a GEO dataset is ‘getGEO’ from the ‘GEOquery’
package.Check how many platforms used for the GSE data, usually there
will only be one platform. We set the first object in the list and gse
now is an expressionSet. You can see that it contains assayData,
phenoData, feature etc.

We can have a look at the sample information, gene annotation, and the
expression data. This allow us to have a rough idea of the information
storing in this expressionSet.

``` r
getwd()
```

    ## [1] "C:/Users/Lynn/Documents/R_GEOdata"

``` r
setwd("C:/Users/Lynn/Documents/R_GEOdata")

###https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html#Introduction
#----import the data------------------------------------
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
my_id <- "GSE63477"
gse <- getGEO(my_id)
```

    ## Found 1 file(s)

    ## GSE63477_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   ID_REF = col_double(),
    ##   GSM1550559 = col_double(),
    ##   GSM1550560 = col_double(),
    ##   GSM1550561 = col_double(),
    ##   GSM1550562 = col_double(),
    ##   GSM1550563 = col_double(),
    ##   GSM1550564 = col_double(),
    ##   GSM1550565 = col_double(),
    ##   GSM1550566 = col_double(),
    ##   GSM1550567 = col_double(),
    ##   GSM1550568 = col_double(),
    ##   GSM1550569 = col_double(),
    ##   GSM1550570 = col_double()
    ## )

    ## File stored at:

    ## C:\Users\Lynn\AppData\Local\Temp\RtmpcbmFsU/GPL16686.soft

    ## Warning: 190 parsing failures.
    ##   row col expected         actual         file
    ## 53792  ID a double AFFX-BioB-3_at literal data
    ## 53793  ID a double AFFX-BioB-3_st literal data
    ## 53794  ID a double AFFX-BioB-5_at literal data
    ## 53795  ID a double AFFX-BioB-5_st literal data
    ## 53796  ID a double AFFX-BioB-M_at literal data
    ## ..... ... ........ .............. ............
    ## See problems(...) for more details.

``` r
## check how many platforms used
length(gse)
```

    ## [1] 1

``` r
gse <-gse[[1]]
gse
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 44629 features, 12 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM1550559 GSM1550560 ... GSM1550570 (12 total)
    ##   varLabels: title geo_accession ... treatment:ch1 (35 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 16657436 16657440 ... 17118478 (44629 total)
    ##   fvarLabels: ID RANGE_STRAND ... RANGE_GB (8 total)
    ##   fvarMetadata: Column Description labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation: GPL16686

``` r
pData(gse)[1:2,] ## print the sample information
```

    ##                                                                        title
    ## GSM1550559 LNCaP CTRL treated in charcoal dextran treated serum, replicate 1
    ## GSM1550560 LNCaP CTRL treated in charcoal dextran treated serum, replicate 2
    ##            geo_accession                status submission_date last_update_date
    ## GSM1550559    GSM1550559 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550560    GSM1550560 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ##            type channel_count source_name_ch1 organism_ch1 characteristics_ch1
    ## GSM1550559  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550560  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ##                                                characteristics_ch1.1
    ## GSM1550559 treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550560 treatment: CTRL treated in charcoal dextran treated serum
    ##                                                      treatment_protocol_ch1
    ## GSM1550559 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550560 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ##                            growth_protocol_ch1 molecule_ch1
    ## GSM1550559 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550560 Cells treated at ca. 80% confluency    total RNA
    ##                                                                                                                                              extract_protocol_ch1
    ## GSM1550559 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550560 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ##            label_ch1
    ## GSM1550559    biotin
    ## GSM1550560    biotin
    ##                                                                                                                      label_protocol_ch1
    ## GSM1550559 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550560 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ##            taxid_ch1
    ## GSM1550559      9606
    ## GSM1550560      9606
    ##                                                                                                                                                                                                                                                hyb_protocol
    ## GSM1550559 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550560 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ##                                   scan_protocol
    ## GSM1550559 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550560 Affymetrix Gene ChIP Scanner 3000 7G
    ##                                                                                                                                    data_processing
    ## GSM1550559 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550560 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ##            data_processing.1 data_processing.2 platform_id   contact_name
    ## GSM1550559 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550560 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ##                                             contact_institute
    ## GSM1550559 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550560 Thomas Jefferson University - Kimmel Cancer Center
    ##                     contact_address contact_city contact_state
    ## GSM1550559 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550560 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ##            contact_zip/postal_code contact_country
    ## GSM1550559                   19107             USA
    ## GSM1550560                   19107             USA
    ##                                                                                               supplementary_file
    ## GSM1550559 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550559/suppl/GSM1550559_01_LN_CDT-CTRL_1.CEL.gz
    ## GSM1550560 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550560/suppl/GSM1550560_01_LN_CTS-CTRL_2.CEL.gz
    ##            data_row_count cell line:ch1
    ## GSM1550559          44629         LNCaP
    ## GSM1550560          44629         LNCaP
    ##                                             treatment:ch1
    ## GSM1550559 CTRL treated in charcoal dextran treated serum
    ## GSM1550560 CTRL treated in charcoal dextran treated serum

``` r
fData(gse)[1,] ## print the gene annotation
```

    ##                ID RANGE_STRAND RANGE_START RANGE_END total_probes    GB_ACC
    ## 16657436 16657436            +       12190     13639           25 NR_046018
    ##                   SPOT_ID     RANGE_GB
    ## 16657436 chr1:12190-13639 NC_000001.10

``` r
exprs(gse)[1,] ## print the expression data
```

    ## GSM1550559 GSM1550560 GSM1550561 GSM1550562 GSM1550563 GSM1550564 GSM1550565 
    ##   24.63215   21.96198   24.36674   24.28032   24.82574   22.72258   25.76430 
    ## GSM1550566 GSM1550567 GSM1550568 GSM1550569 GSM1550570 
    ##   23.03947   24.45452   23.74868   23.77131   21.67176

# Check the normalisation and scales used

We can use this command to check the normalization method, to see if the
data has already processed. So this expression data was RMA normalized
and filtered to remove low-expressing genes. RMA means Robust Multiarray
Average, it is the most common method to determine probeset expression
level for Affymetrix arrays.

The ‘summary’ function can then be used to print the distributions of
expression levels, if the data has been log transformed, typically in
the range of 0 to 16. Hmm…the values are quite big and go beyond 16.
It’s quite weird because RMA should already log2 transform the data at
the last step, but well, let’s do it on our own and move to the next
step. For a more careful analysis, we can try to run the raw data of
this dataset again, by applying RMA normalization on our own, to see if
there is any difference.

Anyway, here, let’s perform a log2 transformation. We may check the
summary of expression level again. And draw a boxplot. We can see that
the distributions of each sample are highly similar, which means the
data have been normalised.

``` r
pData(gse)$data_processing[1]
```

    ## [1] "Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes."

``` r
# For visualisation and statistical analysis, we will inspect the data to 
# discover what scale the data are presented in. The methods we will use assume 
# the data are on a log2 scale; typically in the range of 0 to 16.

## have a look on the expression value
summary(exprs(gse))
```

    ##    GSM1550559        GSM1550560        GSM1550561        GSM1550562     
    ##  Min.   :  17.42   Min.   :  17.49   Min.   :  17.54   Min.   :  17.44  
    ##  1st Qu.:  19.59   1st Qu.:  19.53   1st Qu.:  19.59   1st Qu.:  19.45  
    ##  Median :  22.29   Median :  22.34   Median :  22.39   Median :  22.43  
    ##  Mean   :  48.09   Mean   :  48.51   Mean   :  48.08   Mean   :  48.91  
    ##  3rd Qu.:  33.45   3rd Qu.:  33.90   3rd Qu.:  33.70   3rd Qu.:  33.92  
    ##  Max.   :5889.83   Max.   :6043.74   Max.   :5907.55   Max.   :6087.95  
    ##    GSM1550563        GSM1550564        GSM1550565        GSM1550566     
    ##  Min.   :  17.48   Min.   :  17.47   Min.   :  17.43   Min.   :  17.49  
    ##  1st Qu.:  19.48   1st Qu.:  19.51   1st Qu.:  19.59   1st Qu.:  19.54  
    ##  Median :  22.19   Median :  22.22   Median :  22.38   Median :  22.29  
    ##  Mean   :  49.16   Mean   :  48.72   Mean   :  48.29   Mean   :  48.83  
    ##  3rd Qu.:  33.70   3rd Qu.:  33.84   3rd Qu.:  33.49   3rd Qu.:  33.85  
    ##  Max.   :5991.22   Max.   :5827.09   Max.   :5704.80   Max.   :5938.88  
    ##    GSM1550567        GSM1550568        GSM1550569        GSM1550570     
    ##  Min.   :  17.49   Min.   :  17.32   Min.   :  17.38   Min.   :  17.38  
    ##  1st Qu.:  19.54   1st Qu.:  19.58   1st Qu.:  19.41   1st Qu.:  19.52  
    ##  Median :  22.36   Median :  22.35   Median :  22.31   Median :  22.27  
    ##  Mean   :  48.72   Mean   :  48.31   Mean   :  49.31   Mean   :  48.94  
    ##  3rd Qu.:  33.62   3rd Qu.:  33.58   3rd Qu.:  33.82   3rd Qu.:  33.71  
    ##  Max.   :6140.01   Max.   :6066.52   Max.   :6307.27   Max.   :5844.71

``` r
# From this output we clearly see that the values go beyond 16, 
# so we need to perform a log2 transformation.
exprs(gse) <- log2(exprs(gse))

# check again the summary
summary(exprs(gse))
```

    ##    GSM1550559       GSM1550560       GSM1550561       GSM1550562    
    ##  Min.   : 4.122   Min.   : 4.128   Min.   : 4.133   Min.   : 4.124  
    ##  1st Qu.: 4.292   1st Qu.: 4.288   1st Qu.: 4.292   1st Qu.: 4.282  
    ##  Median : 4.479   Median : 4.482   Median : 4.485   Median : 4.488  
    ##  Mean   : 4.887   Mean   : 4.892   Mean   : 4.887   Mean   : 4.890  
    ##  3rd Qu.: 5.064   3rd Qu.: 5.083   3rd Qu.: 5.075   3rd Qu.: 5.084  
    ##  Max.   :12.524   Max.   :12.561   Max.   :12.528   Max.   :12.572  
    ##    GSM1550563       GSM1550564       GSM1550565       GSM1550566    
    ##  Min.   : 4.128   Min.   : 4.127   Min.   : 4.123   Min.   : 4.128  
    ##  1st Qu.: 4.284   1st Qu.: 4.286   1st Qu.: 4.292   1st Qu.: 4.288  
    ##  Median : 4.472   Median : 4.474   Median : 4.484   Median : 4.478  
    ##  Mean   : 4.893   Mean   : 4.892   Mean   : 4.887   Mean   : 4.892  
    ##  3rd Qu.: 5.075   3rd Qu.: 5.081   3rd Qu.: 5.066   3rd Qu.: 5.081  
    ##  Max.   :12.549   Max.   :12.509   Max.   :12.478   Max.   :12.536  
    ##    GSM1550567       GSM1550568       GSM1550569       GSM1550570    
    ##  Min.   : 4.128   Min.   : 4.114   Min.   : 4.119   Min.   : 4.120  
    ##  1st Qu.: 4.288   1st Qu.: 4.291   1st Qu.: 4.279   1st Qu.: 4.287  
    ##  Median : 4.483   Median : 4.482   Median : 4.479   Median : 4.477  
    ##  Mean   : 4.891   Mean   : 4.888   Mean   : 4.891   Mean   : 4.892  
    ##  3rd Qu.: 5.071   3rd Qu.: 5.070   3rd Qu.: 5.080   3rd Qu.: 5.075  
    ##  Max.   :12.584   Max.   :12.567   Max.   :12.623   Max.   :12.513

``` r
boxplot(exprs(gse),outline=F)
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Inspect the clinical variables

Now we try to look into the pData for the elements that we need for the
analysis. We want to know the sample name, whether it is treatment or
control…in this dataset, the info is stored in the column of
‘characteristics\_ch1.1’.

We can use the select function to subset the column of interest. It will
be useful also to rename the column to something more shorter and
easier.

To make a column of simplified group names for each sample, ‘Stringr’ is
helpful. Two new columns are created, named “group” and “serum”. The
function ‘str\_detect’ is to detect the presence of the words, and then
fill the row accordingly. It totally depends on your dataset to make
these necessary categories in the new columns, just modify these
commands for your dataset of interest.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
sampleInfo <- pData(gse)
head(sampleInfo)
```

    ##                                                                                         title
    ## GSM1550559                  LNCaP CTRL treated in charcoal dextran treated serum, replicate 1
    ## GSM1550560                  LNCaP CTRL treated in charcoal dextran treated serum, replicate 2
    ## GSM1550561 LNCaP 16h cabazitaxel (1nM) treated in charcoal dextran treated serum, replicate 1
    ## GSM1550562 LNCaP 16h cabazitaxel (1nM) treated in charcoal dextran treated serum, replicate 2
    ## GSM1550563   LNCaP 16h docetaxel (1nM) treated in charcoal dextran treated serum, replicate 1
    ## GSM1550564   LNCaP 16h docetaxel (1nM) treated in charcoal dextran treated serum, replicate 2
    ##            geo_accession                status submission_date last_update_date
    ## GSM1550559    GSM1550559 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550560    GSM1550560 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550561    GSM1550561 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550562    GSM1550562 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550563    GSM1550563 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ## GSM1550564    GSM1550564 Public on Jan 01 2015     Nov 19 2014      Jan 01 2015
    ##            type channel_count source_name_ch1 organism_ch1 characteristics_ch1
    ## GSM1550559  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550560  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550561  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550562  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550563  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ## GSM1550564  RNA             1           LNCaP Homo sapiens    cell line: LNCaP
    ##                                                                 characteristics_ch1.1
    ## GSM1550559                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550560                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550561 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550562 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550563   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550564   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum
    ##                                                      treatment_protocol_ch1
    ## GSM1550559 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550560 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550561 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550562 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550563 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ## GSM1550564 16h treatment with 1nM cabazitaxel, docetaxel, or control (EtOH)
    ##                            growth_protocol_ch1 molecule_ch1
    ## GSM1550559 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550560 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550561 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550562 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550563 Cells treated at ca. 80% confluency    total RNA
    ## GSM1550564 Cells treated at ca. 80% confluency    total RNA
    ##                                                                                                                                              extract_protocol_ch1
    ## GSM1550559 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550560 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550561 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550562 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550563 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ## GSM1550564 Standard Trizol RNA extraction, followed by cDNA amplification using the Ovation Pico WTA-system V2 RNA amplification system (NuGen Technologies, Inc)
    ##            label_ch1
    ## GSM1550559    biotin
    ## GSM1550560    biotin
    ## GSM1550561    biotin
    ## GSM1550562    biotin
    ## GSM1550563    biotin
    ## GSM1550564    biotin
    ##                                                                                                                      label_protocol_ch1
    ## GSM1550559 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550560 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550561 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550562 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550563 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ## GSM1550564 5ug cDNA was fragmanted and chemically labeled with biotin using the FL-Ovation cDNA biotin module (NuGen Technologies, Inc)
    ##            taxid_ch1
    ## GSM1550559      9606
    ## GSM1550560      9606
    ## GSM1550561      9606
    ## GSM1550562      9606
    ## GSM1550563      9606
    ## GSM1550564      9606
    ##                                                                                                                                                                                                                                                hyb_protocol
    ## GSM1550559 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550560 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550561 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550562 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550563 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ## GSM1550564 5ug cDNA in 220ul hybridization cocktail was hybridized on the Affymetrix Human Gene 2.0 ST Array in a GeneChip Hybridization Oven 645. Target denaturation was performed at 99C for 2 minutes, 45C for 5min, followed by hybridization for 18h.
    ##                                   scan_protocol
    ## GSM1550559 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550560 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550561 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550562 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550563 Affymetrix Gene ChIP Scanner 3000 7G
    ## GSM1550564 Affymetrix Gene ChIP Scanner 3000 7G
    ##                                                                                                                                    data_processing
    ## GSM1550559 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550560 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550561 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550562 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550563 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ## GSM1550564 Data were processed with GenSpring 11.5 software. The expression data were RMA normalized, and filtered to remove low-expressing genes.
    ##            data_processing.1 data_processing.2 platform_id   contact_name
    ## GSM1550559 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550560 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550561 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550562 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550563 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ## GSM1550564 HuGene-2_0-st.pgf HuGene-2_0-st.mps    GPL16686 Karen,,Knudsen
    ##                                             contact_institute
    ## GSM1550559 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550560 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550561 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550562 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550563 Thomas Jefferson University - Kimmel Cancer Center
    ## GSM1550564 Thomas Jefferson University - Kimmel Cancer Center
    ##                     contact_address contact_city contact_state
    ## GSM1550559 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550560 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550561 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550562 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550563 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ## GSM1550564 233 S 10th St, BLSB 1008 Philadelphia  Pennsylvania
    ##            contact_zip/postal_code contact_country
    ## GSM1550559                   19107             USA
    ## GSM1550560                   19107             USA
    ## GSM1550561                   19107             USA
    ## GSM1550562                   19107             USA
    ## GSM1550563                   19107             USA
    ## GSM1550564                   19107             USA
    ##                                                                                               supplementary_file
    ## GSM1550559 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550559/suppl/GSM1550559_01_LN_CDT-CTRL_1.CEL.gz
    ## GSM1550560 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550560/suppl/GSM1550560_01_LN_CTS-CTRL_2.CEL.gz
    ## GSM1550561 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550561/suppl/GSM1550561_02_LN-CDT_CBTX_2.CEL.gz
    ## GSM1550562 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550562/suppl/GSM1550562_02_LN_CDT-CBTX_1.CEL.gz
    ## GSM1550563 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550563/suppl/GSM1550563_03_LN-CDT-DCTX_1.CEL.gz
    ## GSM1550564 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1550nnn/GSM1550564/suppl/GSM1550564_03_LN-CDT-DCTX_2.CEL.gz
    ##            data_row_count cell line:ch1
    ## GSM1550559          44629         LNCaP
    ## GSM1550560          44629         LNCaP
    ## GSM1550561          44629         LNCaP
    ## GSM1550562          44629         LNCaP
    ## GSM1550563          44629         LNCaP
    ## GSM1550564          44629         LNCaP
    ##                                                              treatment:ch1
    ## GSM1550559                  CTRL treated in charcoal dextran treated serum
    ## GSM1550560                  CTRL treated in charcoal dextran treated serum
    ## GSM1550561 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550562 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550563   16h docetaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550564   16h docetaxel (1nM) treated in charcoal dextran treated serum

``` r
table(sampleInfo$characteristics_ch1.1)
```

    ## 
    ## treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum 
    ##                                                                          2 
    ##                     treatment: 16h cabazitaxel (1nM) treated in full serum 
    ##                                                                          2 
    ##   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum 
    ##                                                                          2 
    ##                       treatment: 16h docetaxel (1nM) treated in full serum 
    ##                                                                          2 
    ##                  treatment: CTRL treated in charcoal dextran treated serum 
    ##                                                                          2 
    ##                                      treatment: CTRL treated in full serum 
    ##                                                                          2

``` r
#Let's pick just those columns that seem to contain factors we might 
#need for the analysis.
sampleInfo <- select(sampleInfo, characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo, sample = characteristics_ch1.1)

head(sampleInfo)
```

    ##                                                                                sample
    ## GSM1550559                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550560                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550561 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550562 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550563   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550564   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum

``` r
dim(sampleInfo)
```

    ## [1] 12  1

``` r
sampleInfo$sample
```

    ##  [1] "treatment: CTRL treated in charcoal dextran treated serum"                 
    ##  [2] "treatment: CTRL treated in charcoal dextran treated serum"                 
    ##  [3] "treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum"
    ##  [4] "treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum"
    ##  [5] "treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum"  
    ##  [6] "treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum"  
    ##  [7] "treatment: CTRL treated in full serum"                                     
    ##  [8] "treatment: CTRL treated in full serum"                                     
    ##  [9] "treatment: 16h cabazitaxel (1nM) treated in full serum"                    
    ## [10] "treatment: 16h cabazitaxel (1nM) treated in full serum"                    
    ## [11] "treatment: 16h docetaxel (1nM) treated in full serum"                      
    ## [12] "treatment: 16h docetaxel (1nM) treated in full serum"

``` r
library(stringr)
sampleInfo$group <- ""
for(i in 1:nrow(sampleInfo)){
  if(str_detect(sampleInfo$sample[i], "CTRL") && str_detect(sampleInfo$sample[i], "full"))
  {sampleInfo$group[i] <- "Conf"}
  
  if(str_detect(sampleInfo$sample[i], "CTRL") && str_detect(sampleInfo$sample[i], "dextran"))
  {sampleInfo$group[i] <- "Cond"}
  
  if(str_detect(sampleInfo$sample[i], "cabazitaxel") && str_detect(sampleInfo$sample[i], "full"))
  {sampleInfo$group[i] <- "cabazitaxelf"}
  
  if(str_detect(sampleInfo$sample[i], "cabazitaxel") && str_detect(sampleInfo$sample[i], "dextran"))
  {sampleInfo$group[i] <- "cabazitaxeld"}
  
  if(str_detect(sampleInfo$sample[i], "docetaxel") && str_detect(sampleInfo$sample[i], "full"))
  {sampleInfo$group[i] <- "docetaxelf"}
  
  if(str_detect(sampleInfo$sample[i], "docetaxel") && str_detect(sampleInfo$sample[i], "dextran"))
  {sampleInfo$group[i] <- "docetaxeld"}
}

sampleInfo 
```

    ##                                                                                sample
    ## GSM1550559                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550560                  treatment: CTRL treated in charcoal dextran treated serum
    ## GSM1550561 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550562 treatment: 16h cabazitaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550563   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550564   treatment: 16h docetaxel (1nM) treated in charcoal dextran treated serum
    ## GSM1550565                                      treatment: CTRL treated in full serum
    ## GSM1550566                                      treatment: CTRL treated in full serum
    ## GSM1550567                     treatment: 16h cabazitaxel (1nM) treated in full serum
    ## GSM1550568                     treatment: 16h cabazitaxel (1nM) treated in full serum
    ## GSM1550569                       treatment: 16h docetaxel (1nM) treated in full serum
    ## GSM1550570                       treatment: 16h docetaxel (1nM) treated in full serum
    ##                   group
    ## GSM1550559         Cond
    ## GSM1550560         Cond
    ## GSM1550561 cabazitaxeld
    ## GSM1550562 cabazitaxeld
    ## GSM1550563   docetaxeld
    ## GSM1550564   docetaxeld
    ## GSM1550565         Conf
    ## GSM1550566         Conf
    ## GSM1550567 cabazitaxelf
    ## GSM1550568 cabazitaxelf
    ## GSM1550569   docetaxelf
    ## GSM1550570   docetaxelf

``` r
sampleInfo$serum <- ""
for(i in 1:nrow(sampleInfo)){
  if(str_detect(sampleInfo$sample[i], "dextran"))
  {sampleInfo$serum[i] <- "dextran"}
  
  if(str_detect(sampleInfo$sample[i], "full"))
  {sampleInfo$serum[i] <- "full_serum"}
 
}

sampleInfo <- sampleInfo[,-1]
sampleInfo
```

    ##                   group      serum
    ## GSM1550559         Cond    dextran
    ## GSM1550560         Cond    dextran
    ## GSM1550561 cabazitaxeld    dextran
    ## GSM1550562 cabazitaxeld    dextran
    ## GSM1550563   docetaxeld    dextran
    ## GSM1550564   docetaxeld    dextran
    ## GSM1550565         Conf full_serum
    ## GSM1550566         Conf full_serum
    ## GSM1550567 cabazitaxelf full_serum
    ## GSM1550568 cabazitaxelf full_serum
    ## GSM1550569   docetaxelf full_serum
    ## GSM1550570   docetaxelf full_serum

# Sample clustering and Principal Components Analaysis

We can visualize the correlations between the samples by hierarchical
clustering.

The function ‘cor’ can calculate the correlation on the scale of 0 to 1,
in a pairwise fashion between all samples, then visualise on a heatmap.
There are many ways to create heatmaps in R, here I use ‘pheatmap’, the
only argument it requires is a matrix of numeric values.

We can add more sample info onto the plot to get a better pic of the
group and clustering. Here, we make use of the ‘sampleInfo’ file that
was created earlier, to match with the columns of the correlation
matrix.

``` r
library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)   
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
## Print the rownames of the sample information and check it matches the correlation matrix

rownames(sampleInfo)
```

    ##  [1] "GSM1550559" "GSM1550560" "GSM1550561" "GSM1550562" "GSM1550563"
    ##  [6] "GSM1550564" "GSM1550565" "GSM1550566" "GSM1550567" "GSM1550568"
    ## [11] "GSM1550569" "GSM1550570"

``` r
colnames(corMatrix)
```

    ##  [1] "GSM1550559" "GSM1550560" "GSM1550561" "GSM1550562" "GSM1550563"
    ##  [6] "GSM1550564" "GSM1550565" "GSM1550566" "GSM1550567" "GSM1550568"
    ## [11] "GSM1550569" "GSM1550570"

``` r
## If not, force the rownames to match the columns
#rownames(sampleInfo) <- colnames(corMatrix)

pheatmap(corMatrix, annotation_col= sampleInfo)
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

Another way is to use Principal component analysis (PCA). It has to note
that the data has to be transposed, so that the genelist is in the
column, while rownames are the samples, so the PCA process will not run
out of the memory in the oher way round.

Let’s add labels to plot the results, here, we use the ‘ggplots2’
package, while the ‘ggrepel’ package is used to position the text labels
more cleverly so they can be read. Here we can see that the samples are
divided into two groups based on the serum treatment types.

``` r
#make PCA
library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group, label=paste("",group))) + geom_point() + geom_text_repel()
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Differential expression analysis

In this section, we use the limma package to perform differential
expressions. Limma stands for “Linear models for microarray”. Here, we
need to tell limma what sample groups we want to compare. I choose
sampleInfo$group. A design matrix will be created, this is a matrix of 0
and 1, one row for each sample and one column for each sample group.

We can rename the column names so that it is easier to see.

Now, let’s check if the expression data contain any lowly-expressed
genes, this will affect the quality of DE analysis. A big problem in
doing statistical analysis like limma is the inference of type 1
statistical errors, also called false positive. One simple way to reduce
the possibility for type 1 errors is to do fewer comparisons, by
filtering the data. For example, we know that not all genes are
expressed in all tissues and many genes will not be expressed in any
sample. As a result, in DGE analysis, it makes sense to remove the genes
that are likely not expressed at all.

It is quite subjective how one defines a gene being expressed, here, I
follow the tutorial, to make the cut off at the median of the expression
values, which means to consider around 50% of the genes will not be
expressed. Keep those expressed genes if they are present in more than 2
samples.

We can see that around half of the genes are not qualified as an
“expressed” gene here, which makes sense, bcoz our cut-off is the
median value.

``` r
library(limma)
```

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
design <- model.matrix(~0 + sampleInfo$group)
design
```

    ##    sampleInfo$groupcabazitaxeld sampleInfo$groupcabazitaxelf
    ## 1                             0                            0
    ## 2                             0                            0
    ## 3                             1                            0
    ## 4                             1                            0
    ## 5                             0                            0
    ## 6                             0                            0
    ## 7                             0                            0
    ## 8                             0                            0
    ## 9                             0                            1
    ## 10                            0                            1
    ## 11                            0                            0
    ## 12                            0                            0
    ##    sampleInfo$groupCond sampleInfo$groupConf sampleInfo$groupdocetaxeld
    ## 1                     1                    0                          0
    ## 2                     1                    0                          0
    ## 3                     0                    0                          0
    ## 4                     0                    0                          0
    ## 5                     0                    0                          1
    ## 6                     0                    0                          1
    ## 7                     0                    1                          0
    ## 8                     0                    1                          0
    ## 9                     0                    0                          0
    ## 10                    0                    0                          0
    ## 11                    0                    0                          0
    ## 12                    0                    0                          0
    ##    sampleInfo$groupdocetaxelf
    ## 1                           0
    ## 2                           0
    ## 3                           0
    ## 4                           0
    ## 5                           0
    ## 6                           0
    ## 7                           0
    ## 8                           0
    ## 9                           0
    ## 10                          0
    ## 11                          1
    ## 12                          1
    ## attr(,"assign")
    ## [1] 1 1 1 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$`sampleInfo$group`
    ## [1] "contr.treatment"

``` r
## the column names are a bit ugly, so we will rename
colnames(design) <- c("Cabazitaxeld","Cabazitaxelf","Cond","Conf","Docetaxeld","Docetaxelf")

design
```

    ##    Cabazitaxeld Cabazitaxelf Cond Conf Docetaxeld Docetaxelf
    ## 1             0            0    1    0          0          0
    ## 2             0            0    1    0          0          0
    ## 3             1            0    0    0          0          0
    ## 4             1            0    0    0          0          0
    ## 5             0            0    0    0          1          0
    ## 6             0            0    0    0          1          0
    ## 7             0            0    0    1          0          0
    ## 8             0            0    0    1          0          0
    ## 9             0            1    0    0          0          0
    ## 10            0            1    0    0          0          0
    ## 11            0            0    0    0          0          1
    ## 12            0            0    0    0          0          1
    ## attr(,"assign")
    ## [1] 1 1 1 1 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$`sampleInfo$group`
    ## [1] "contr.treatment"

``` r
## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 3

## check how many genes are removed / retained.
table(keep)
```

    ## keep
    ## FALSE  TRUE 
    ## 20965 23664

``` r
## subset to just those expressed genes
gse <- gse[keep,]
```

Here there is a little extra step to find out the outliers. This has to
be done carefully so the filtered data won’t be too biased. We calculate
‘weights’ to define the reliability of each sample. The ‘arrayweights’
function will assign a score to each sample, with a value of 1 implying
equal weight. Samples with score less than 1 are down-weighed, or else
up-weighed.

``` r
# coping with outliers
## calculate relative array weights
aw <- arrayWeights(exprs(gse),design)
aw
```

    ##         1         2         3         4         5         6         7         8 
    ## 0.9704842 0.9704842 0.8790788 0.8790788 0.9799805 0.9799805 0.8296341 0.8296341 
    ##         9        10        11        12 
    ## 1.1632620 1.1632620 1.2393733 1.2393733

Now we have a design matrix, we need to estimate the coefficients. For
this design, we will essentially average the replicate arrays for each
sample level. In addition, we will calculate standard deviations for
each gene, and the average intensity for the genes across all
microarrays.

We are ready to tell limma which pairwise contrasts that we want to
make. For this experiment, we are going to contrast treatment (there are
two types of texane drugs) and control in each serum type. So there are
4 contrasts to specify.

To do the statistical comparisons, Limma uses Bayesian statistics to
minimize type 1 error. The eBayes function performs the tests. To
summarize the results of the statistical test, ‘topTable’ will adjust
the p-values and return the top genes that meet the cutoffs that you
supply as arguments; while ‘decideTests’ will make calls for DEGs by
adjusting the p-values and applying a logFC cutoff similar to topTable.

``` r
## Fitting the coefficients
fit <- lmFit(exprs(gse), design,
             weights = aw)

head(fit$coefficients)
```

    ##          Cabazitaxeld Cabazitaxelf     Cond     Conf Docetaxeld Docetaxelf
    ## 16657436     4.604278     4.590902 4.539703 4.606668   4.569910   4.504447
    ## 16657440     5.073400     5.194166 5.274991 5.001295   5.230574   5.016408
    ## 16657450     6.606782     6.532054 6.571707 6.754964   6.854709   6.495389
    ## 16657469     5.318536     5.298899 5.358629 5.191262   5.398528   5.356415
    ## 16657476     5.410581     5.443410 5.358263 5.343908   5.516173   5.335881
    ## 16657480     4.481099     4.454008 4.487220 4.453781   4.418727   4.379055

``` r
## Making comparisons between samples, can define multiple contrasts
contrasts <- makeContrasts(Docetaxeld - Cond, Cabazitaxeld - Cond, Docetaxelf - Conf, Cabazitaxelf - Conf, levels = design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)


topTable(fit2)
```

    ##          Docetaxeld...Cond Cabazitaxeld...Cond Docetaxelf...Conf
    ## 16681891      -0.427866483          0.55288766       -0.02322651
    ## 16840609       0.573090433         -0.20083754       -0.21477778
    ## 17017165       0.116986540         -0.32582451       -0.44959204
    ## 16782010       0.178206104         -0.01891042       -0.29766696
    ## 17099705      -0.437302094          0.02759174       -0.16952565
    ## 16691877      -0.149430506         -0.18664917       -0.90866405
    ## 16959582      -0.006190599         -0.17371327        0.40620159
    ## 16970902      -0.416408859          0.01679171       -0.02017417
    ## 16936214      -0.112315413          0.12707501       -0.28387277
    ## 17009126      -0.645309985         -0.59167505       -0.39348344
    ##          Cabazitaxelf...Conf  AveExpr        F      P.Value  adj.P.Val
    ## 16681891          -0.1414628 5.625598 43.01843 2.992480e-06 0.07081404
    ## 16840609          -0.3278230 5.769203 18.87958 1.217387e-04 0.99909331
    ## 17017165          -0.3710638 5.244400 14.18263 4.059196e-04 0.99909331
    ## 16782010           0.2850105 4.440952 14.12466 4.128114e-04 0.99909331
    ## 17099705           0.2776297 6.750407 13.62541 4.783683e-04 0.99909331
    ## 16691877          -0.5849898 4.667038 11.52373 9.376105e-04 0.99909331
    ## 16959582           0.3098904 5.531461 11.26699 1.024647e-03 0.99909331
    ## 16970902          -0.1618155 5.237839 11.18970 1.052724e-03 0.99909331
    ## 16936214          -0.6271607 5.482977 11.17492 1.058196e-03 0.99909331
    ## 17009126          -0.3389764 4.940635 11.00507 1.123616e-03 0.99909331

``` r
topTable1 <- topTable(fit2, coef=1)
topTable2 <- topTable(fit2, coef=2)
topTable3 <- topTable(fit2, coef=3)
topTable4 <- topTable(fit2, coef=4)

#if we want to know how many genes are differentially expressed overall, we can use the decideTest function.
summary(decideTests(fit2))
```

    ##        Docetaxeld - Cond Cabazitaxeld - Cond Docetaxelf - Conf
    ## Down                   0                   0                 0
    ## NotSig             23664               23664             23664
    ## Up                     0                   0                 0
    ##        Cabazitaxelf - Conf
    ## Down                     0
    ## NotSig               23664
    ## Up                       0

``` r
table(decideTests(fit2))
```

    ## 
    ##     0 
    ## 94656

# Further visualization with gene annotation

Now we want to know the gene name associated with the gene ID. The
annotation data can be retrieved with the ‘fData’ function. Let’s select
the ID, GB\_ACC, this is genbank accession ID. Add into fit2 table.

The “Volcano Plot” function is a common way of visualising the results
of a DE analysis. The x axis shows the log-fold change and the y axis is
some measure of statistical significance, which in this case is the
log-odds, or “B” statistic. We can also change the color of those genes
with p value cutoff more than 0.05, and fold change cut off more than 1.

``` r
anno <- fData(gse)
head(anno)
```

    ##                ID RANGE_STRAND RANGE_START RANGE_END total_probes    GB_ACC
    ## 16657436 16657436            +       12190     13639           25 NR_046018
    ## 16657440 16657440            +       29554     31109           28          
    ## 16657450 16657450            +      317811    328581           36 NR_024368
    ## 16657469 16657469            +      329790    342507           27          
    ## 16657476 16657476            +      459656    461954           27 NR_029406
    ## 16657480 16657480            +      523009    532878           12          
    ##                     SPOT_ID     RANGE_GB
    ## 16657436   chr1:12190-13639 NC_000001.10
    ## 16657440   chr1:29554-31109 NC_000001.10
    ## 16657450 chr1:317811-328581 NC_000001.10
    ## 16657469 chr1:329790-342507 NC_000001.10
    ## 16657476 chr1:459656-461954 NC_000001.10
    ## 16657480 chr1:523009-532878 NC_000001.10

``` r
anno <- select(anno,ID,GB_ACC)
fit2$genes <- anno

topTable(fit2)
```

    ##                ID       GB_ACC Docetaxeld...Cond Cabazitaxeld...Cond
    ## 16681891 16681891 NM_001013692      -0.427866483          0.55288766
    ## 16840609 16840609                    0.573090433         -0.20083754
    ## 17017165 17017165                    0.116986540         -0.32582451
    ## 16782010 16782010                    0.178206104         -0.01891042
    ## 17099705 17099705    NR_039696      -0.437302094          0.02759174
    ## 16691877 16691877                   -0.149430506         -0.18664917
    ## 16959582 16959582                   -0.006190599         -0.17371327
    ## 16970902 16970902                   -0.416408859          0.01679171
    ## 16936214 16936214    NR_037440      -0.112315413          0.12707501
    ## 17009126 17009126                   -0.645309985         -0.59167505
    ##          Docetaxelf...Conf Cabazitaxelf...Conf  AveExpr        F      P.Value
    ## 16681891       -0.02322651          -0.1414628 5.625598 43.01843 2.992480e-06
    ## 16840609       -0.21477778          -0.3278230 5.769203 18.87958 1.217387e-04
    ## 17017165       -0.44959204          -0.3710638 5.244400 14.18263 4.059196e-04
    ## 16782010       -0.29766696           0.2850105 4.440952 14.12466 4.128114e-04
    ## 17099705       -0.16952565           0.2776297 6.750407 13.62541 4.783683e-04
    ## 16691877       -0.90866405          -0.5849898 4.667038 11.52373 9.376105e-04
    ## 16959582        0.40620159           0.3098904 5.531461 11.26699 1.024647e-03
    ## 16970902       -0.02017417          -0.1618155 5.237839 11.18970 1.052724e-03
    ## 16936214       -0.28387277          -0.6271607 5.482977 11.17492 1.058196e-03
    ## 17009126       -0.39348344          -0.3389764 4.940635 11.00507 1.123616e-03
    ##           adj.P.Val
    ## 16681891 0.07081404
    ## 16840609 0.99909331
    ## 17017165 0.99909331
    ## 16782010 0.99909331
    ## 17099705 0.99909331
    ## 16691877 0.99909331
    ## 16959582 0.99909331
    ## 16970902 0.99909331
    ## 16936214 0.99909331
    ## 17009126 0.99909331

``` r
## Create volcano plot
full_results1 <- topTable(fit2, coef=1, number=Inf)
library(ggplot2)
ggplot(full_results1,aes(x = logFC, y=B)) + geom_point()
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1


full_results1 %>% 
  mutate(Significant = P.Value < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
```

![](GSE_analysis_microarray_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

# Further visualization of selected gene

I think at this point, we are quite clear about data structure of GSE
data. It has an experiment data, pData; the expression data, exprs; and
also annotation data, fData. And we have learned how to check the
expression data, normalize them, and perform differential expression
analysis.

Now, with the differential expression gene tables, there are some
downstream analyses that we can continue, such as to export a full table
of DE genes, to generate a heatmap for your selected genes, get the gene
list for a particular pathway, or survival analysis (but this is only
for those clinical data).

Here, I just want to look into the fold change data of a selected gene,
whether it is significantly differential expressed or not.

``` r
## Get the results for particular gene of interest
#GB_ACC for Nkx3-1 is NM_001256339 or NM_006167
##no NM_001256339 in this data
full_results2 <- topTable(fit2, coef=2, number=Inf)
full_results3 <- topTable(fit2, coef=3, number=Inf)
full_results4 <- topTable(fit2, coef=4, number=Inf)
filter(full_results1, GB_ACC == "NM_006167")
```

    ##                ID    GB_ACC      logFC  AveExpr         t    P.Value adj.P.Val
    ## 17075536 17075536 NM_006167 -0.2379718 8.055042 -3.056052 0.01219317 0.9819604
    ##                  B
    ## 17075536 -3.374675

``` r
filter(full_results2, GB_ACC == "NM_006167")
```

    ##                ID    GB_ACC       logFC  AveExpr         t   P.Value adj.P.Val
    ## 17075536 17075536 NM_006167 -0.09964001 8.055042 -1.244539 0.2418169 0.9999406
    ##                  B
    ## 17075536 -4.561219

``` r
filter(full_results3, GB_ACC == "NM_006167")
```

    ##                ID    GB_ACC      logFC  AveExpr         t   P.Value adj.P.Val
    ## 17075536 17075536 NM_006167 0.02902121 8.055042 0.3762533 0.7146269 0.9999674
    ##                  B
    ## 17075536 -4.926393

``` r
filter(full_results4, GB_ACC == "NM_006167")
```

    ##                ID    GB_ACC      logFC  AveExpr         t   P.Value adj.P.Val
    ## 17075536 17075536 NM_006167 0.01019426 8.055042 0.1304658 0.8987979 0.9998834
    ##                 B
    ## 17075536 -4.95523

That’s all for the walk-through, thanks for reading, I hope you have
learned something new here.

# Acknowlegdement

Many thanks to the following tutorials made publicly available:

1.  Introduction to microarray analysis GSE15947, by Department of
    Statistics, Purdue Univrsity
    <https://www.stat.purdue.edu/bigtap/online/docs/Introduction_to_Microarray_Analysis_GSE15947.html>

2.  Mark Dunning, 2020, GEO tutorial, by Sheffield Bioinformatics Core
    <https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html#Further_processing_and_visualisation_of_DE_results>
