# QTLSorghum
QTLseqr is an R package for QTL mapping using NGS Bulk Segregant
Analysis.

QTLseqr is still under development and is offered with out any
guarantee.

### **For more detailed instructions please read the vignette [here](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf)**

### For updates read the [NEWS.md](https://github.com/bmansfeld/QTLseqr/blob/master/NEWS.md)

# Installation

<!-- You can install and update QTLseqr by using our [drat](http://dirk.eddelbuettel.com/code/drat.html) repository hosted on our github page: -->

<!-- ```{r drat-install, eval = FALSE} -->

<!-- install.packages("QTLseqr", repos = "http://bmansfeld.github.io/drat") -->

<!-- ``` -->

<!-- OR You can install QTLseqr from github with: -->

You can install QTLseqr from github with:

``` r
# install devtools first to download packages from github
install.packages("devtools")

# use devtools to install QTLseqr
devtools::install_github("PBGLMichaelHall/QTLseqr")
```


**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).

**If you use QTLseqr in published research, please cite:**

> Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
> analysis with next-generation sequencing *The Plant Genome*
> [doi:10.3835/plantgenome2018.01.0006](https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006)

We also recommend citing the paper for the corresponding method you work
with.

QTL-seq method:

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
> C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
> L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
> quantitative trait loci in rice by whole genome resequencing of DNA
> from two bulked populations. *Plant J*, 74: 174–183.
> [doi:10.1111/tpj.12105](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105)

G prime method:

> Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
> Segregant Analysis Using Next Generation Sequencing. *PLOS
> Computational Biology* 7(11): e1002255.
> [doi.org/10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)

## Abstract

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is
efficient in detecting quantitative trait loci (QTL). Despite the
popularity of NGS-BSA and the R statistical platform, no R packages are
currently available for NGS-BSA. We present QTLseqr, an R package for
NGS-BSA that identifies QTL using two statistical approaches: QTL-seq
and G’. These approaches use a simulation method and a tricube smoothed
G statistic, respectively, to identify and assess statistical
significance of QTL. QTLseqr, can import and filter SNP data, calculate
SNP distributions, relative allele frequencies, G’ values, and
log10(p-values), enabling identification and plotting of QTL.

# Examples:

# Load/install libraries


``` r {r libraries}
# install.packages("tinytex")
# install.packages("vcfR")
# install.packages("tidyr")
# install.packages("ggplot2")
devtools::install_github("PBGLMichaelHall/QTLseqr",force = TRUE)
library(QTLseqr)
library(tinytex)
library(vcfR)
library(tidyr)
library(ggplot2)

```
# Set the Working Directory to where VCF file is stored in file system

``` r 
setwd("/home/michael/Desktop/QTLseqr/extdata")
```
# Vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTLseqR functions will only take SNPS, ie, length of REF and ALT== 1
```r
vcf <- read.vcfR(file = "freebayes_D2.filtered.vcf")
```


![Screenshot from 2022-03-30 15-12-19](https://user-images.githubusercontent.com/93121277/160842876-d35bcbdf-c487-42ad-ac92-01f98f436eea.png)

```r

#Convert to tidy data frame
VCF_TIDY <- vcfR2tidy(vcf)
```

![Screenshot from 2022-03-30 15-18-37](https://user-images.githubusercontent.com/93121277/160843944-e37e77e9-acb4-401f-95eb-f75f894e953f.png)

# Call the Parser
```r
QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "D2_F2_tt",LowBulk = "D2_F2_TT", filename = Hall)
```
![Screenshot from 2022-03-30 15-10-38](https://user-images.githubusercontent.com/93121277/160842389-5d1d3915-c652-4fa4-a0a4-2829d04ad9a0.png)


# Preview the CSV file
![mcsv](https://user-images.githubusercontent.com/93121277/158783968-db377510-7852-4359-a48f-afb34b8efb5a.png)

# Invoke unique command to extract Sample names reverse comapatible to the VCF

```r
unique(VCF_TIDY$gt$Indiv)
```
![Screenshot from 2022-03-30 15-33-29](https://user-images.githubusercontent.com/93121277/160846874-b44284ed-9eb5-44a7-9ef0-65a5819e182e.png)



```r
#Set High bulk and Low bulk sample names and parser generated file name
#The file name is generated from the QTLParser_1_MH function in line 119

HighBulk <- "D2_F2_tt"
LowBulk <- "D2_F2_TT"
file <- "Hall.csv"

#Choose which chromosomes/contigs will be included in the analysis,

Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")

df <-
  importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  ) 

```
![Screenshot from 2022-03-30 15-38-14](https://user-images.githubusercontent.com/93121277/160847939-5d0feb2d-8d54-4d59-949d-be1939df8de2.png)

# Inspect the head of the df object
![Screenshot from 2022-03-30 15-41-45](https://user-images.githubusercontent.com/93121277/160848653-41ab84f1-1d37-49aa-9bd4-ea1319704c8f.png)


```r
#plot histograms associated with filtering arguments such as mamximum and minumum Total Depths and reference Allele Frequency to determine cut off values 

ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)

ggsave(filename = "Depth_Histogram.png",plot=last_plot())

```

![hist34](https://user-images.githubusercontent.com/93121277/156784666-7abcc556-e0ef-4b4b-b981-c97074fccb0c.png)

```r

ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))

ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

```

![hist4](https://user-images.githubusercontent.com/93121277/156784371-a96fcac0-fe21-43de-ad6c-568c69e3d25a.png)



``` r {r Filtering, warning = FALSE}

#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    #    minGQ = 0
  )
```
![Screenshot from 2022-04-01 08-29-37](https://user-images.githubusercontent.com/93121277/161207708-a50c7061-94ce-4729-b764-417bdb488d32.png)


```r

#Run G' analysis
df_filt<-runGprimeAnalysis_MH(
  SNPset = df_filt,
  windowSize = 5000000,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)
```
![Screenshot from 2022-04-01 08-31-03](https://user-images.githubusercontent.com/93121277/161207885-45119458-fa34-4259-80be-a8cd2f30fb3c.png)

# G' Distribution Plot
```r
#The plot reveals a skewed G Prime statistic with a really small variance. Perhaps it is due to the small number of variants called.
#In addition, Hampels outlier filter in the second argument, can also be changed to "deltaSNP"
plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
```
![Screenshot from 2022-04-01 08-53-51](https://user-images.githubusercontent.com/93121277/161211346-780ed554-5fa6-4d22-bada-d446714c06aa.png)

```r
#We can see raw data before and after our filtering step
plotGprimeDist_MH(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)
```
![Screenshot from 2022-04-01 08-59-04](https://user-images.githubusercontent.com/93121277/161212125-9daf552a-ddb1-4a50-8603-08adbf21a869.png)


```




``` r {r QTLSEQ, warning = FALSE}

#Run QTLseq analysis
df_filt2 <- runQTLseqAnalysis_MH(
  SNPset = df_filt,
  windowSize = 5000000,
  popStruc = "F2",
  bulkSize = c(45, 38),
  replications = 10000,
  intervals = c(95, 99)
)

```
![Screenshot from 2022-04-01 08-39-28](https://user-images.githubusercontent.com/93121277/161209279-3bd76a3c-97d4-44f0-801d-c5c3e7ba1942.png)


# Plot G Statistic Distribution

```r
hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")
```

![Screenshot from 2022-04-01 08-44-38](https://user-images.githubusercontent.com/93121277/161209927-59cd31d5-5663-4d0b-9d7d-ae6ee64a180d.png)



```r
#Plot Snps as a function of chromosome and position values
plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggsave(filename = "nSNPs.png",plot = last_plot())
```
![snps](https://user-images.githubusercontent.com/93121277/158374828-5878811f-20c0-449c-8ac3-10d195bd3375.png)

```r
#Using QTLStats funciton plot Gprime Statistic with False Discovery Rate Threhshold as a third argument boolean operator as TRUE. The q value is used as FDR threshold null value is 0.05%.
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggsave(filename = "GPrime.png",plot = last_plot())
```
![Gprime](https://user-images.githubusercontent.com/93121277/156785067-1253709f-b19d-4326-a4c3-8db8abf33052.png)


```r
#Again using plotQTLStats change second argument varaible to deltaSNP and plot.
plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())
```

![deltasnp](https://user-images.githubusercontent.com/93121277/158375244-b509026f-8f57-4368-a18e-93711afe8661.png)

```r
#Finally with plotQTLStats plot negLog10Pval
plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01)
ggsave(filename = "negLog10Pval.png",plot = last_plot())
```

![neglog](https://user-images.githubusercontent.com/93121277/158375824-99eefc2d-1d09-4a25-926d-7a1d4a0e1340.png)

```r
#Add subset argument to focus on particular chromosomes one, three, four, and six.
#The reason is due to signficant QTL regions
plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("1","3","4","6"))

```

![Gprime2](https://user-images.githubusercontent.com/93121277/158375088-768456af-870e-4ec3-9427-900ce6aa8ea7.png)


# Export summary CSV



```r
QTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
# Preview the Summary QTL

![Screenshot from 2022-04-01 09-18-54](https://user-images.githubusercontent.com/93121277/161214771-b0ded302-846c-43c9-a8a4-a87b97423f63.png)

``` r {r AlleleFreq,warning=FALSE}
#Use the function to plot allele frequencies per chromosome
#Second argument size specifes size of scalar factor on nSNPs
Obs_Allele_Freq(SNPSet = df_filt, size = .001)
```

![Screenshot from 2022-04-01 12-38-33](https://user-images.githubusercontent.com/93121277/161247921-2fab6d12-6b11-433a-af3b-a5c969aabc8a.png)



``` r
##Use the function to investigate chromosomal region of interest
Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = 4, threshold = .90)
```
![lb2](https://user-images.githubusercontent.com/93121277/158378100-43ec8daa-ce2a-475b-8515-4a6f0195924e.png)



![freq2](https://user-images.githubusercontent.com/93121277/158391011-7fcc46ae-bc6f-4f48-a11d-24ff59b8bdbe.png)

################################################################################################################################################

# Rice QTL Analysis: High Bulk sample size of 430 Tolerant to cold environments and Low Bulk sample size of 385 Suseptilble to cold environments
```r


#Set Working Directory
setwd("/home/michael/Desktop/RiceCold2")

#vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTL-Rice-Cold functions will only take SNPS, ie, length of REF and ALT== 1
vcf <- read.vcfR(file = "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf.gz")

```
![Screenshot from 2022-04-01 09-27-24](https://user-images.githubusercontent.com/93121277/161215916-c3328eaf-721b-4778-940a-e356fe60e9ca.png)

```r
#Convert to tidy data frame
VCF_TIDY <- vcfR2tidy(vcf)
```
![Screenshot from 2022-04-01 09-37-35](https://user-images.githubusercontent.com/93121277/161217479-f8a99317-1dec-4b7b-a7de-04d488abbdf0.png)

```r
#Call the Parser
QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "ET-pool-385",LowBulk = "ES-pool-430", filename = "Hall")
```
# Standard Parser Output makes a list of Chromosome as named in VCF file
![Screenshot from 2022-04-01 09-41-00](https://user-images.githubusercontent.com/93121277/161218164-169207c9-7039-4857-9066-f59f07c206b7.png)
# Also find unique sample names reverse compatible with VCF file
```r
unique(VCF_TIDY$gt$Indiv)
```
![Screenshot from 2022-04-01 09-42-58](https://user-images.githubusercontent.com/93121277/161218336-79ca6524-e2ba-47fd-bbd8-68125bdd966b.png)

```r
#Set High bulk and Low bulk sample names and parser generated file name

HighBulk <- "ET-pool-385"
LowBulk <- "ES-pool-430"
file <- "Hall.csv"
```
# Preview the CSV file

![Screenshot from 2022-04-01 10-33-41](https://user-images.githubusercontent.com/93121277/161226740-5b648928-97bd-40e4-8622-e8f845bbe271.png)


```r
#Input chromosomes values which will be included in the analysis,
Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")



df <-
  importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  ) 

```

![Screenshot from 2022-04-01 10-35-10](https://user-images.githubusercontent.com/93121277/161226938-8f50aa1b-27b4-4fe1-9b77-b06d14ec8b31.png)



# Plot histograms associated with filtering arguments to determine if cut off values are appropriate




```r

ggplot(data =df) +geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)
ggsave(filename = "Depth_Histogram.png",plot=last_plot())

```
![dplowhigh](https://user-images.githubusercontent.com/93121277/158780363-0939b60a-4a19-4104-b435-e352d715f5df.png)

```r

ggplot(data = df) +geom_histogram(aes(x = REF_FRQ))
ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

```

![reffreq](https://user-images.githubusercontent.com/93121277/158780481-24a6992e-e239-41d4-9868-3b6da831e757.png)



```r


#Filter SNPs based on some criteria
df_filt <-
  filterSNPs(
    SNPset = df,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    #    minGQ = 0
  )
  
  ```
  
  ![Screenshot from 2022-04-01 09-54-07](https://user-images.githubusercontent.com/93121277/161220139-0c079197-99db-4292-949a-cd23d0e2e7c7.png)


```r
#Run G' analysis
df_filt<-runGprimeAnalysis_MH(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)


```
![Screenshot from 2022-04-01 09-54-50](https://user-images.githubusercontent.com/93121277/161220373-00be3d7e-b67d-44ca-a1af-5c4670751f39.png)


```r


#Run QTLseq analysis
df_filt2 <- runQTLseqAnalysis_MH(
  SNPset = df_filt,
  windowSize = 1e6,
  popStruc = "F2",
  bulkSize = c(430, 385),
  replications = 10000,
  intervals = c(95, 99)
)

```
![Screenshot from 2022-04-01 09-56-09](https://user-images.githubusercontent.com/93121277/161220504-ba5f8e90-126f-4ee7-9a1d-06d5e09e6a94.png)

# Plot G Statistic Distribution
```r
hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")
```
![gstat](https://user-images.githubusercontent.com/93121277/158780626-0dd9efaa-8c2b-448e-8e22-ce94a1ff6fbf.png)


```r

# G' Distribution Plot
plotGprimeDist_MH(SNPset = df_filt2, outlierFilter = "Hampel")
ggsave(filename = "Hampel_GPrime.png",plot = last_plot())

```

![gprime](https://user-images.githubusercontent.com/93121277/158780745-ce684a8b-5267-42f2-aab1-9038b66490a5.png)


```r


plotGprimeDist_MH(SNPset = df_filt2, outlierFilter = "deltaSNP",filterThreshold = 0.1)
ggsave(filename = "DeltaSNP.png",plot = last_plot())

```

![deltaSNP](https://user-images.githubusercontent.com/93121277/158780846-58095997-5814-4440-8594-e2c560412eee.png)


```r


#make the Plot
plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggsave(filename = "nSNPs.png",plot = last_plot())


```
![Screenshot from 2022-04-01 10-36-37](https://user-images.githubusercontent.com/93121277/161227187-9a3d3138-e996-4965-9659-6a879d0be70a.png)



```


```r
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggsave(filename = "GPrime.png",plot = last_plot())

```

![Screenshot from 2022-04-01 10-37-50](https://user-images.githubusercontent.com/93121277/161227403-aa546c4a-c6f0-4c3c-8c25-485bae654739.png)


```r
plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())

```
![Screenshot from 2022-04-01 10-39-24](https://user-images.githubusercontent.com/93121277/161227673-447d11f9-1245-4042-affd-56ca43c39799.png)





```r
plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))
ggsave(filename = "negLog10Pval.png",plot = last_plot())

```
![Screenshot from 2022-04-01 10-05-35](https://user-images.githubusercontent.com/93121277/161222325-69f4ad58-8fb8-4eb7-bc91-fc900be47416.png)


```r

plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("NC_029256.1","NC_029257.1","NC_029263.1","NC_029265.1"))

```
![Screenshot from 2022-04-01 10-06-19](https://user-images.githubusercontent.com/93121277/161222339-dbbafc12-0631-4d22-b477-5f038256ee6e.png)


# Export summary CSV
```r

getQTLTable(SNPset = df_filt2, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
# Preview the QTL Summary
![Screenshot from 2022-04-01 09-58-43](https://user-images.githubusercontent.com/93121277/161220947-979d5bbf-8438-4110-a950-a33224878a01.png)

```r

#Use the function to plot allele frequencies per chromosome
Obs_Allele_Freq(SNPSet = df_filt)
```

# Looks Dense
![LowB](https://user-images.githubusercontent.com/93121277/158788612-ebe92c64-ed0b-48f8-9ef6-f7e8d6f7bd2d.png)
![HighB](https://user-images.githubusercontent.com/93121277/158788656-ee153e0a-1a6a-4b28-86b8-f7f218bb1287.png)


# Filter Low Allelic Depth Frequencies
```r
##Use the function to investigate chromosomal region of interest
Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = 8, threshold = .75)

```

# Preview the plot with idenitfied SNP positions
![LB3](https://user-images.githubusercontent.com/93121277/158788921-b35622eb-9926-4c0b-9b6b-fd06c71dc0c2.png)
![HB3](https://user-images.githubusercontent.com/93121277/158789001-69b463a5-56df-44f6-94a7-7d016c795833.png)

# Investigate SNP@POS 24525659
![header](https://user-images.githubusercontent.com/93121277/158789386-a163d764-607b-4a4c-8884-11ce4a9debb0.png)
![values](https://user-images.githubusercontent.com/93121277/158789420-806f5cfe-ed44-48b7-80fa-f738eaf7a5e8.png)



