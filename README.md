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
#plot histograms associated with filtering arguments to determine if cut off values are appropriate




ggplot(data =df) + geom_histogram(aes(x = DP.LOW + DP.HIGH)) + xlim(0,400)

ggsave(filename = "Depth_Histogram.png",plot=last_plot())

ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))

ggsave(filename = "Ref_Freq_Histogram.png",plot = last_plot())

```


![hist34](https://user-images.githubusercontent.com/93121277/156784666-7abcc556-e0ef-4b4b-b981-c97074fccb0c.png)
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


#Run G' analysis
df_filt<-runGprimeAnalysis_MH(
  SNPset = df_filt,
  windowSize = 5000000,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)


#Plot G Statistic Distribution
hist(df_filt2$G,breaks = 950,xlim = c(0,10),xlab = "G Distribution",main = "Histogram of G Values")


# G' Distribution Plot
plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
ggsave(filename = "Hampel_GPrime.png",plot = last_plot())

setwd("/home/michael/Desktop/SorghumQTL/DeltaSNP/")
plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP",filterThreshold = 0.1)
ggsave(filename = "DeltaSNP.png",plot = last_plot())

```
![Gstat](https://user-images.githubusercontent.com/93121277/156784882-f70f2144-88dd-4a87-ab6e-8a04a31c7fe1.png)



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



#make the Plot
snpnumber <- plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggsave(filename = "nSNPs.png",plot = last_plot())

setwd("/home/michael/Desktop/SorghumQTL/GPrimeDistributionPlots/")
Gprime<-plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggsave(filename = "GPrime.png",plot = last_plot())
setwd("/home/michael/Desktop/SorghumQTL/DeltaSNP/")
deltaSNP<-plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())
setwd("/home/michael/Desktop/SorghumQTL/negLog10Pval/")
neglog<-plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01)
ggsave(filename = "negLog10Pval.png",plot = last_plot())
Gprime2<-plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("1","3","4","6"))
#plotiing snps per Chromosome
snpnumber

```
![snps](https://user-images.githubusercontent.com/93121277/158374828-5878811f-20c0-449c-8ac3-10d195bd3375.png)



Gprime

![Gprime](https://user-images.githubusercontent.com/93121277/156785067-1253709f-b19d-4326-a4c3-8db8abf33052.png)

Gprime2

![Gprime2](https://user-images.githubusercontent.com/93121277/158375088-768456af-870e-4ec3-9427-900ce6aa8ea7.png)

deltaSNP

![deltasnp](https://user-images.githubusercontent.com/93121277/158375244-b509026f-8f57-4368-a18e-93711afe8661.png)

neglog

![neglog](https://user-images.githubusercontent.com/93121277/158375824-99eefc2d-1d09-4a25-926d-7a1d4a0e1340.png)



# Export summary CSV



```r
QTLTable <- getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
write.csv(QTLTable, file = "QTLTablePeaks.csv", row.names = FALSE, col.names = TRUE)
Table4 <- read.table(file = "QTLTablePeaks.csv",header = TRUE, sep = ",", fill=TRUE)
```
# Preview the Summary QTL


``` r {r AlleleFreq,warning=FALSE}
#Use the function to plot allele frequencies per chromosome
Obs_Allele_Freq(SNPSet = df_filt)
```

![lb](https://user-images.githubusercontent.com/93121277/158377569-1b9183b6-2419-4eba-811d-7d95f9cc36ce.png)


![FREQ1](https://user-images.githubusercontent.com/93121277/156783048-88f98c3c-4c89-4519-bf1b-7ed15b6b5dea.png)




``` r
##Use the function to investigate chromosomal region of interest
Obs_Allele_Freq2(SNPSet = df_filt, ChromosomeValue = 4, threshold = .90)
```
![lb2](https://user-images.githubusercontent.com/93121277/158378100-43ec8daa-ce2a-475b-8515-4a6f0195924e.png)



![freq2](https://user-images.githubusercontent.com/93121277/158391011-7fcc46ae-bc6f-4f48-a11d-24ff59b8bdbe.png)


# Rice QTL Analysis: High Bulk sample size of 430 Tolerant to cold environments and Low Bulk sample size of 385 Suseptilble to cold environments
```r


#Set Working Directory
setwd("/home/michael/Desktop/RiceCold2")

#vcf file must only contain bialleleic variants. (filter upstream, e.g., with bcftools view -m2 -M2), also the QTL-Rice-Cold functions will only take SNPS, ie, length of REF and ALT== 1
vcf <- read.vcfR(file = "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf.gz")

#Convert to tidy data frame
VCF_TIDY <- vcfR2tidy(vcf)

#Call the Parser
QTLParser_1_MH(vcf = VCF_TIDY, HighBulk = "ET-pool-385",LowBulk = "ES-pool-430")


#Set High bulk and Low bulk sample names and parser generated file name

HighBulk <- "ET-pool-385"
LowBulk <- "ES-pool-430"
file <- "Hall.csv"
```
# Preview the CSV file
![csv](https://user-images.githubusercontent.com/93121277/158783309-96b56b84-18f7-4221-b3b2-58ffee1be281.png)


```r
#Choose which chromosomes will be included in the analysis,
#the tidy data frame makes a CHROMKEY so no need to change chromosome names
Chroms <- 1:12


df <-
  importFromTable(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  ) 


#plot histograms associated with filtering arguments to determine if cut off values are appropriate



ggplot(data = df) +geom_histogram(aes(x = AD_ALT.LOW + AD_ALT.HIGH)) + xlim(0,400)
ggsave(filename = "AD_Histogram.png",plot = last_plot())
```

![altlowhigh](https://user-images.githubusercontent.com/93121277/158780030-405d4137-d84e-40e5-a994-c2afd82844ef.png)


```r

ggplot(data = df) +geom_histogram(aes(x = AD_REF.LOW + AD_REF.HIGH)) + xlim(0,400)
ggsave(filename = "AD_Ref_Histogram.png",plot = last_plot())

```

![reflow](https://user-images.githubusercontent.com/93121277/158780242-92fa67d2-1929-450a-ae57-5868938c8349.png)

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


#Run G' analysis
df_filt<-runGprimeAnalysis_MH(
  SNPset = df_filt,
  windowSize = 1e6,
  outlierFilter = "deltaSNP",
  filterThreshold = 0.1)


```


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



#Plot G Statistic Distribution
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
snpnumber <- plotQTLStats(SNPset = df_filt2, var = "nSNPs")
ggsave(filename = "nSNPs.png",plot = last_plot())
snpnumber

```
![snps9](https://user-images.githubusercontent.com/93121277/158781092-45cd72ca-ca9d-4b52-8a94-d0202ba7d6a9.png)

```


```r
Gprime<-plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
ggsave(filename = "GPrime.png",plot = last_plot())

```

![GPrime](https://user-images.githubusercontent.com/93121277/158781260-207f2f26-85b7-424e-9c54-df7b65784283.png)

```r
deltaSNP<-plotQTLStats(SNPset = df_filt2, var = "deltaSNP", plotIntervals  = TRUE)
ggsave(filename = "DeltaSNPInterval.png",plot = last_plot())
deltaSNP
```
![deltaSNP8](https://user-images.githubusercontent.com/93121277/158781414-bf2d0243-01bf-4cfb-a053-7fac67d895b4.png)



```r
neglog<-plotQTLStats(SNPset = df_filt2, var = "negLog10Pval",plotThreshold = TRUE,q=0.01,subset = c("1","2","8","10"))
ggsave(filename = "negLog10Pval.png",plot = last_plot())
neglog
```
![neglog'](https://user-images.githubusercontent.com/93121277/158781581-badbc13f-405e-4ff7-8fe8-68b91f421811.png)

```r

Gprime2<-plotQTLStats(SNPset = df_filt2, var = "Gprime",plotThreshold = TRUE,q=0.01,subset = c("1","2","8","10"))
Gprime2
```

![Gprime22](https://user-images.githubusercontent.com/93121277/158781815-8700bad5-fc58-4ce6-a2fe-72ecba9d66cd.png)


# Export summary CSV
```r

QTLTable <- getQTLTable(SNPset = df_filt2, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
write.csv(QTLTable, file = "QTLTablePeaks.csv", row.names = FALSE, col.names = TRUE)
Table4 <- read.table(file = "QTLTablePeaks.csv",header = TRUE, sep = ",", fill=TRUE)
```
# Preview the QTL Summary



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


```r


obs_MH(SNPSet = df_filt2, ChromosomeValue1 = 1,ChromosomeValue2 = 2,ChromosomeValue3 = 8,ChromosomeValue4 = 10, threshold = .01)
for(i in 1:12){
  obs_MH(SNPSet = df_filt2, ChromosomeValue1 = i,ChromosomeValue2 = i,ChromosomeValue3 = i,ChromosomeValue4 = i, threshold = .01)
}

```

# From the function I choose A few plots to show an interesting relationship between High Bulk allelic depth frequencies and Gprime Values

![Screenshot from 2022-03-17 11-26-49](https://user-images.githubusercontent.com/93121277/158790139-43e3e639-08d6-4c6a-8cc5-b4b0f88b94a9.png)


