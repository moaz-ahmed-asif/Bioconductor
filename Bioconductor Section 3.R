BiocManager::install(c("Biobase",
                       "GEOquery",
                       "genomicsclass/GSE5859Subset",
                       "affy",
                       "hgu95acdf",
                       "genefilter",
                       "parathyroidSE",
                       "airway",
                       "pasillaBamSubset",
                       "Rsamtools",
                       "GenomicAlignments",
                       "ArrayExpress",
                       "NGScopyData",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "AnnotationDbi"))

library(Biobase)
library(GEOquery)

geoq <- getGEO("GSE9514")    # download a microarray dataset from GEO
names(geoq)    
e <- geoq[[1]]    # extract ExpressionSet
e

dim(e)    # number of features and samples in ExpressionSet
ncol(e)
nrow(e)

exprs(e)[1:3,1:3] # exprs gives matrix of microarray values
head(exprs(e))[,1]    # first column
exprs(e)[1,]    # first row
exprs(e)["10000_at",]    # can also index by name
rownames(e)[1]    # row names are probe sets
dim(exprs(e))    # rows are features, columns are samples
summary(exprs(e))
boxplot(log2(exprs(e)),outline=FALSE)

pData(e)[1:3,1:6] # pData gives phenotype data (sample information)
names(pData(e))
pData(e)$characteristics_ch1 # column in GEO to describe experimental state/condition
as.numeric(as.factor(pData(e)$characteristics_ch1)) # help see replicates of each state
dim(pData(e))    # rows of pData correspond to columns of exprs
dim(e)

library(pheatmap)
pheatmap(cor(exprs(e), use = "c")) # argument use = "c" stops an error if there are any missing data points

# fData gives feature data (probe information)
fData(e)[1:3,1:3]
dim(fData(e))    # rows of fData correspond to rows of exprs
names(fData(e))
head(as.factor(fData(e)$"Gene Symbol"))
head(rownames(e))

# additional annotation tied to ExpressionSet
experimentData(e)
annotation(e)

library(GSE5859Subset)
data(GSE5859Subset)
dim(geneExpression)
dim(sampleInfo)
dim(geneAnnotation)

pd = AnnotatedDataFrame(sampleInfo)
rownames(pd) = colnames(geneExpression)
pData(pd)["GSM136530.CEL.gz","date"]
varLabels(pd)[1]
pData(pd)

fd = AnnotatedDataFrame(geneAnnotation)
rownames(fd) = rownames(geneExpression)
pData(fd)["204810_s_at", "CHR"] # What chromosome is this gene from?

eset = ExpressionSet(geneExpression,phenoData = pd,featureData = fd)
ind1 <- which(featureData(eset)$CHR == "chrY")
ind2 <- pData(eset)$group == 1
femaleY <- colMeans(exprs(eset)[ind1, ind2]) 
maleY <- colMeans(exprs(eset)[ind1, !ind2]) 
boxplot(maleY, femaleY)
median(maleY) - median(femaleY)
eset[1:10,1:5] # first 10 features, first 5 samples



wd <- getwd()
datadir <- paste0(wd, "/rawdata-master")
basedir <- paste0(datadir, "/celfiles")
setwd(basedir)
library(affy)
tab <- read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)
rownames(tab) <- tab$filenames
tab
fns <- list.celfiles(basedir)
fns
fns %in% tab[,1] # check
ab <- ReadAffy(phenoData=tab)
dim(pm(ab))
dim(pData(ab))
annotation(ab)
rownames(pData(ab))
colnames(pm(ab))
pData(ab[,1:4])
e <- rma(ab) # preprocess probe-level data into gene-level data
dim(exprs(e))
ejust <- justRMA(phenoData=tab) # read and process data to gene-level in one command
dim(ejust)

library(limma)
library(rafalib)
basedir <- paste0(datadir, "/agilent")
setwd(basedir)
targets <- readTargets("TargetBeta7.txt")
targets
RG <- read.maimages(targets$FileName, source="genepix")
MA <- MA.RG(RG,bc.method="none")
dim(RG$R)
dim(RG$G)
dim(MA$M)
dim(MA$A)
plot(MA$A[,1], MA$M[,1]) # MA plot for first sample
mypar(1,1)
imageplot(MA$M[,2], RG$printer, zlim=c(-3,3))

detach("package:affy")
library(oligo)
basedir <- paste0(datadir,"/celfiles")
setwd(basedir)
tab <- read.delim("sampleinfo.txt",check.names=FALSE,as.is=TRUE)
fns <- list.celfiles(listGzipped=TRUE)
fns %in% tab[,1] ##check
pd <- as(tab, "AnnotatedDataFrame")
efs <- read.celfiles(filenames=tab[,1],phenoData=pd,sampleNames=sampleNames(pd))
e <- rma(efs)
dim(exprs(e))
setwd(wd)



library(parathyroidSE)
data(parathyroidGenesSE)
se <- parathyroidGenesSE
se
# assay contains results of the assay 
dim(se)
assay(se)[1:3,1:3]
dim(assay(se))    # rows = features (ranges), columns = samples
# colData contains sample information
colData(se)[1:3,1:6]
dim(colData(se))
names(colData(se))
colData(se)$treatment
as.numeric(colData(se)$treatment)
# rowRanges contains feature information
rowRanges(se)[1]
class(rowRanges(se))
length(rowRanges(se))
length(rowRanges(se)[[1]])
head(rownames(se))
metadata(rowRanges(se))
# additional metadata, including sample information
metadata(se)$MIAME
abstract(metadata(se)$MIAME)



library(airway)
data(airway)
airway
as.numeric(metadata(airway)[[1]]@pubMedIds) # What is the PubMed ID (PMID) for the original paper supplying these data?
as.numeric(airway@metadata[[1]]@pubMedIds)
dim(airway)
rowRanges(airway)
length(rowRanges(airway)) # How many features are in the dataset?
colData(airway)@nrows # How many samples are in the dataset?
names(colData(airway))
which(colnames(airway)=="SRR1039509")
as.character(colData(airway)[which(airway$Run=="SRR1039509"),]$cell) # Which cell line is associated with sample SRR1039509?
colData(airway) # What is the name of the metadata column specifying whether the sample was treated with dexamethasone (a steroid)?
"trt" %in% colData(airway)$dex # TRUE
as.character(unique(colData(airway)$dex)) 
min((colData(airway))$avgLength)
row.names(colData(airway)[which((colData(airway))$avgLength==min(colData(airway)$avgLength)),]) # Which sample (column namme has the shortest average read length)
as.character(colData(airway)[which((colData(airway))$avgLength==min(colData(airway)$avgLength)),]$Run)
length(rowRanges(airway)) # Each feature (row) in the dataset corresponds to a single gene. How many genes are in the dataset?
ir <- rowRanges(airway)[[100]]
length(ir) # The row metadata for the airway object include GRanges specifying the exon boundaries for each gene. How many exons are in the 100th gene?
as.numeric(unique(seqnames(ir))) # On which chromosome is the 100th gene?
sum(width(ranges(ir)))
width(ranges(ir))
intersect(ranges(ir), ranges(ir))
intersect(ir, ir)
ir[ir %over% ir]
sum(seqlengths(ir))
sum(width(reduce(ranges(ir))))
width(range(reduce(ranges(ir)))) # How many bases long is the 100th gene (including introns)?
width(range(ir)) 
start(resize(range(ir),1)) # What is the transcription start site (TSS) of the 100th gene?
end(range(ir)) # for negative strands, end point is the start point
mean(assay(airway)["ENSG00000103196",]) # What is the mean expression level of the gene with ENSEMBL ID "ENSG00000103196" across all samples?
trt <- mean(assay(airway)["ENSG00000103196",
                          which(colData(airway)$dex=="trt")]) # What is the mean expression level of this gene in samples treated with dexamethasone?
untrt <- mean(assay(airway)["ENSG00000103196",
                          which(colData(airway)$dex=="untrt")]) # What is the mean expression level of this gene in untreated control samples?
log2(trt / untrt) # What is the log (base 2) ratio of mean expression of this gene between treated and untreated samples?

# Next Generation Sequencing (NGS)

library(pasillaBamSubset)
library(Rsamtools)
filename <- untreated1_chr4()
(bf <- BamFile(filename))
seqinfo(bf)
(sl <- seqlengths(bf))
quickBamFlagSummary(bf)
(gr <- GRanges("chr4",IRanges(1, sl["chr4"])))
countBam(bf, param=ScanBamParam(which = gr))
reads <- scanBam(BamFile(filename, yieldSize=5))
class(reads)
names(reads[[1]])
reads[[1]]$pos # the aligned start position
reads[[1]]$rname # the chromosome
reads[[1]]$strand # the strand
reads[[1]]$qwidth # the width of the read
reads[[1]]$seq # the sequence of the read
gr <- GRanges("chr4",IRanges(500000, 700000))
reads <- scanBam(bf, param=ScanBamParam(what=c("pos","strand"), which=gr))
hist(reads[[1]]$pos)
readsByStrand <- split(reads[[1]]$pos, reads[[1]]$strand)
myHist <- function(x) table(cut(x, 50:70 * 10000 ))
tab <- sapply(readsByStrand, myHist)
barplot( t(tab) )

gr <- GRanges("chr4",IRanges(440000, 470000))
reads <- scanBam(bf, param=ScanBamParam(what=c("seq"), which=gr))
# How are the start positions distributed?
reads[[1]]$seq
class(reads[[1]]$seq)
# GC content = proportion of the sequence with is 'G' or 'C', a number between 0 and 1
total_GC <- sum(letterFrequency(reads[[1]]$seq, "GC"))
total <- sum(width(reads[[1]]$seq))
GC_content <- total_GC / total
AT_content <- 1 - GC_content

library(GenomicAlignments)
(ga <- readGAlignments(bf))
length(ga)
granges(ga[1])
gr <- GRanges("chr4", IRanges(700000, 800000))
(fo <- findOverlaps(ga, gr)) # which reads over this range
countOverlaps(gr, ga) # count overlaps of range with the reads
table(ga %over% gr) # logical vector of read overlaps with the range
countOverlaps(ga, gr) # integer vector with the number of overlaps for each read with the range in `gr`

ga <- readGAlignments(BamFile(filename))
hist(start(ga), breaks=100)

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
g <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
length(g)
g2 <- g[g %over% GRanges("chr4", IRanges(200000, 300000))]
g2
class(g2['FBgn0039890'])
findOverlaps(g,g2["FBgn0039890"])
countOverlaps(g2["FBgn0039890"], ga)



library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

grl <- exonsBy(txdb, by="gene") # make GRangesList of exons for each gene
grl[100] # GRangesList of exons for 100th gene
grl[[100]] # GRanges with exons of 100th gene
grl[[100]][1] # first exon of 100th gene

# paths to BAM files
fl1 <- untreated1_chr4()
fl2 <- untreated3_chr4()

# libraries for importing BAM files
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)

# specify files with BamFileList
fls <- BamFileList(c(fl1, fl2))
names(fls) <- c("first","second")

# find reads that overlap exons
so1 <- summarizeOverlaps(features=grl,
                         reads=fls,
                         ignore.strand=TRUE)
so1

# examine count matrix
head(assay(so1))
colSums(assay(so1))

# examine rest of SummarizedExperiment components
rowRanges(so1) # feature or genes
colData(so1) # samples
colData(so1)$sample <- c("one","two")    # add sample information
colData(so1)
metadata(rowRanges(so1)) 

# exploratory data analysis of counts
x <- assay(so1)[,1]
hist(x[x > 0], col="grey")
hist(x[x > 0 & x < 40000], col="grey")
plot(assay(so1) + 1, log="xy")

# count second file as paired-end reads
# ?untreated3_chr4
# ?summarizeOverlaps
fls <- BamFileList(fl2)
so2 <- summarizeOverlaps(features=grl,
                         reads=fls,
                         ignore.strand=TRUE,
                         singleEnd=FALSE, 
                         fragments=TRUE)
colSums(assay(so2))
colSums(assay(so1))

# show there are half as many reads in so2 as so1
plot(assay(so1)[,2], assay(so2)[,1], xlim=c(0,5000), ylim=c(0,5000),
     xlab="single end counting", ylab="paired end counting")
abline(0,1)
abline(0,.5)

g <- genes(txdb)
seqnames(g)
g <- g[seqnames(g) == "chr4"]
seqnames(g)
grl <- exonsBy(txdb, by="gene")
grl <- grl[names(g)]
all.equal(names(g), names(grl))

library(Rsamtools)
bf <- BamFile(untreated1_chr4())
so1 <- summarizeOverlaps(features=g,
                         reads=bf,
                        ignore.strand=TRUE)
so2 <- summarizeOverlaps(features=grl,
                         reads=bf,
                       ignore.strand=TRUE)

x = assay(so1)[,1]
y = assay(so2)[,1]
mean ((y/x)[x>0]) # removing zero counts

plot(assay(so1) + 1, assay(so2) + 1, log="xy")
abline(0,1)

propy <- y/sum(y)
# Now multiply these proportions by 1 million. This operation scales each column of the count table such that we get the number of reads expected if the sample were sequenced to have 1 million reads mapping to genes on chromosome 4
propym <- propy*1e6
head(propym,1)
propym[["FBgn0002521"]]

ebp = sum(width(reduce(grl)))
mean(ebp)
summary(ebp)

head(propym)

count = assay(so2)
fpm = (count/sum(count)) * 1e6
head(fpm)
fpkm = (propym/ebp) * 1e3
head(fpkm)
fpkm[["FBgn0002521"]]



library(GEOquery)
glioMA = getGEO("GSE78703")[[1]]
glioMA

names(pData(glioMA))    # variable names
glioMA$molecule_ch1    # molecule being assayed (RNA)
table(glioMA$`treated with:ch1`, glioMA$`cell type:ch1`)    # experimental variables



library(hgu133plus2.db)
hgu133plus2.db
library(hgu133plus2probe)
head(hgu133plus2probe)
dim(hgu133plus2probe)
select(hgu133plus2.db, keytype="PROBEID", 
       columns=c("SYMBOL", "GENENAME", "PATH", "MAP"), keys="1007_s_at")
length(unique(select(hgu133plus2.db, keytype="SYMBOL", 
       columns=c("PROBEID", "GENENAME", "PATH", "MAP"), keys="EGFR")$PROBEID))

library(GO.db)
keytypes(GO.db)
columns(GO.db)
Gos <- select(GO.db, keytype="TERM", 
       columns=c("GOID"),keys ="glial cell proliferation")$GOID
Gos
keytypes(hgu133plus2.db)
columns(hgu133plus2.db)
length(unique(select(hgu133plus2.db, keytype="GO", 
                     columns=c("PROBEID"), keys=Gos)$PROBEID))



# Given an ExpressionSet x
exprs(x) # access experimental data
pData(x) # access sample information
fData(x) # access feature information

# Given a SummarizedExperiment x
assay(x) # access experimental data
colData(x) # access sample information
rowData(x) # access feature information

library(GEOquery)
library(NGScopyData)
library(Rsamtools)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
# mp = getGEO("GSE3725")[[1]]

# How many samples are in the ExpressionSet obtained above?
nrow(pData(mp))
# How many features?
nrow(fData(mp))

# remove 6 samples corresponding to experimental controls
pData(mp)$title
lstem <- mp[, !grepl("^GMP expressing", pData(mp)$title)]

pData(lstem)$title
titles <- as.character(pData(lstem)$title)
cell_type <- gsub(".*\\((.*?)( enriched)?\\).*", "\\1", titles)
# add cell_type column to pData
pData(lstem)$cell_type <- factor(cell_type)

# How many samples are from the cell type L-GMP?
lstem
head(pData(lstem))
length(which(pData(lstem)$cell_type=="L-GMP"))

# Inspect the feature data. What "Gene Symbol" field is associated with probe "1421579_at" ?
head(fData(lstem),1)
fData(lstem)$ID
fData(lstem)[which(fData(lstem)$ID == "1421579_at"),]$"Gene Symbol"

# What is the mean expression level of "1421579_at" across all samples?
mean(exprs(lstem)["1421579_at",])

# What is the mean expression level of "1421579_at" in samples with cell_type "L-GMP"?
rownames(pData(lstem)[which(pData(lstem)$cell_type == "L-GMP"),])
mean(exprs(lstem)["1421579_at",rownames(pData(lstem)[which(pData(lstem)$cell_type == "L-GMP"),])])



library(Rsamtools)
library(NGScopyData)
tps_27 <- tps_27.chr6()$bamFpath # path to BAM file
tps_27 <- BamFile(tps_27)
