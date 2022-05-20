library(dslabs)
data(heights)
heights
class(heights)
class(heights$sex)
class(heights$height)
class("Male")
class(75.00000)
nrow(heights)
heights$height[777]
heights[777,2]
heights$sex[777]
heights[1, 777]
heights[777,1]
max(heights$height)
which.min(heights$height)
mean(heights$height)
median(heights$height)
mean(heights$sex == "Male")
sum(heights$height >= 78)
sum(heights$sex == "Female" & heights$height > 78)



library(Biobase) # load one of the core Bioconductor packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)

# get help through the documentation
help.start()
?mean
help(mean)
help(package="genefilter")

# inspect objects, classes and methods

?ExpressionSet
?"ExpressionSet-class"
methods(class = ExpressionSet)

# inspect the source code for functions and methods
read.csv
plotMA
showMethods("plotMA")
getMethod("plotMA","data.frame")

# vignettes teach you how to use various functions in a package
vignette(package="Biobase")
vignette("ExpressionSetIntroduction")
browseVignettes(package="Biobase")

# report key details about your R session
sessionInfo()

BiocManager::version()

# section 1
# options(timeout = 3000)
BiocManager::install(c("genefu",
                       "COPDSexualDimorphism.data",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression",
                       "SNPlocs.Hsapiens.dbSNP144.GRCh37"))
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]
sum(is.na(sig.gene70$NCBI.gene.symbol)) # How many components of the signature have a missing value for the associated NCBI gene symbol?
index <- grep("kinase",sig.gene70$Description)
length(index) # How many of the members of the 70-gene signature are genes coding for kinases?
sig.gene70$Description[index]

library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)
expr.meta
sum(expr.meta$GENDER=='2-Female') # What is the number of female participants in this study?
median(expr.meta$pkyrs) # What is the median of the distribution of pack years smoked in this cohort (women and men)?

# True or False: The distribution of pack-years smoked is well-approximated by a Gaussian (Normal) probability distribution
hist(expr.meta$pkyrs)
qqnorm(expr.meta$pkyrs)
qqline(expr.meta$pkyrs)

# Which of the following is an aspect of the display that would suggest caution in using the t test in comparing males and females with respect to pack years smoked?
boxplot(pkyrs~gender, data=expr.meta)
boxplot(age~gender, data=expr.meta)

# Under this model we use a number denoted lambda that for our purposes is used as an exponent to transform the dependent variable pyp1 of the regression to have a distribution that is approximately Gaussian. Thus, if lambda is 1, we use pyp1 untransformed, if lambda is 0.5, we use sqrt(pyp1), and so on
# For what value of lambda does the likelihood reach its highest value for the model lm1?
expr.meta$pyp1 = expr.meta$pkyrs+1
library(MASS)
lm1 = lm(pyp1~gender, data=expr.meta)
boxcox(lm1) # ~ 0.5

lambda <- 0.43
boxplot(I(pyp1^lambda)~gender, data=expr.meta)
# for 0.34 ~ 0.43, no outlier

library(gwascat)
data(ebicat_2020_04_30)
ebicat_2020_04_30
sort(table(mcols(ebicat_2020_04_30)[,"DISEASE/TRAIT"]),decreasing=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19
chr11seq <- BSgenome.Hsapiens.UCSC.hg19[["chr11"]]
subseq(chr11seq,start=10^6,width=25)
# which of the following sequences is most common on chromosome 11: "ATG", "TGA", "TAA", and "TAG"
?countPattern
countPattern("ATG",chr11seq)
countPattern("TGA",chr11seq)
countPattern("TAA",chr11seq)
countPattern("TAG",chr11seq) # or
seqs <- c("ATG","TGA","TAA","TAG")
seqs[as.numeric(which.max(sapply(seqs, function(sq){countPattern(sq,chr11seq)})))]

chr7seq <- BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
?alphabetFrequency
# what percent of chromosome 7 is T,C,G, and A
alphabetFrequency(chr7seq,as.prob=TRUE)[c("A","C","G","T","N")]*100 # %
# What proportion are Cs (including counts of N in the total)
sum(alphabetFrequency(chr7seq,as.prob=TRUE)[c("C","N")])

library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snps144 = SNPlocs.Hsapiens.dbSNP144.GRCh37
class(snps144)
s17 = snpsBySeqname(snps144, "17")
head(s17)
class(s17)
head(s17$RefSNP_id)
snpsById(snps144, "rs73971683")
s17[which(s17$RefSNP_id=="rs73971683")]



library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)

library(SummarizedExperiment)
tissSE = SummarizedExperiment(list(rma=e))
colData(tissSE) = DataFrame(tab)

tissSE
assay(tissSE)
mean(assay(tissSE["209169_at",]))
boxplot(assay(tissSE)["209169_at",]~tissSE$Tissue, mar=c(6,4,2,2) )
# This gene is expressed in the brain but not the other tissues
IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")
# Which of the following ID(s) appear to represent a gene specific to placenta
tissSE[which(tissSE$Tissue=="placenta")]@NAMES
tissSE[which(tissSE$Tissue=="placenta" && tissSE@NAMES==IDs)]

library(rafalib)
mypar(3,2)
sapply(IDs,function(x){boxplot(e[x,]~tissue,las=2,main=x)})

library(hgu133a.db)
sym = mapIds(hgu133a.db, keys=rownames(tissSE), column="SYMBOL", keytype="PROBEID")
nm = mapIds(hgu133a.db, keys=rownames(tissSE), column="GENENAME", keytype="PROBEID")
rowData(tissSE) = DataFrame(symbol=sym, genename=nm)

grep("H2AX", rowData(tissSE)$symbol, value=TRUE)
length((grep("^H2A", rowData(tissSE)$symbol, value=TRUE)))
assay(tissSE["205436_s_at",])

mypar()
boxplot(as.numeric(assay(tissSE["205436_s_at",]))~tissSE$Tissue)
# Expression of H2AFX is greater in hippocampus than in cerebellum

emails <- c("john.doe@ivyleague.edu", "education@world.gov", "dalai.lama@peace.org",
            "invalid.edu", "quant@bigdatacollege.edu", "cookie.monster@sesame.tv")
hits <- grep('edu',emails)
# Subset emails using hits
emails[hits]
grep("GAPDH", rowData(tissSE)$symbol, value=TRUE)
sum(rowData(tissSE)$symbol=='H2AX',na.rm=TRUE) # How many features in this SummarizedExperiment measure expression of gene H2AFX?

nrow(tissSE[grep("kinase", rowData(tissSE)$genename),])

library(hgu133aprobe)
head(hgu133aprobe)
class(hgu133aprobe)
length(which(hgu133aprobe$Probe.Set.Name=="206269_at"))

