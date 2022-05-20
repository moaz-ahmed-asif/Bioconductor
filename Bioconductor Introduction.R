# options(timeout = 3000)
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
                       "BSgenome.Mmusculus.UCSC.mm10",
                       "RMariaDB",
                       "org.Hs.eg.db",
                       "TxDb.Mmusculus.UCSC.mm10.knownGene",
                       "ShortRead",
                       "Rfastp",
                       "Rsubread",
                       "BSgenome.Hsapiens.UCSC.hg38",
                       "Rbowtie2",
                       "QuasR"))

library(BSgenome.Mmusculus.UCSC.mm10)
class(BSgenome.Mmusculus.UCSC.mm10)
BSgenome.Mmusculus.UCSC.mm10
contigNames <- seqnames(BSgenome.Mmusculus.UCSC.mm10)
contigNames[1:22]
contigLengths <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
contigLengths[1:4]
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10$chr19
chr19_Seq
chr19_Seq <- BSgenome.Mmusculus.UCSC.mm10[["chr19"]]
chr19_Seq
class(chr19_Seq)
system.time(subseq(chr19_Seq,start=1,end=10000000))
system.time(chr19_Seq[1:10000000])
subseq(chr19_Seq,start=1,end=10000000)

library(Biostrings)
alphabetFrequency(chr19_Seq) 
chr19_SeqComp <- complement(chr19_Seq)
chr19_SeqRev <- reverse(chr19_Seq)
chr19_SeqRevComp <- reverseComplement(chr19_Seq[10000000:10000010])
chr19_Seq[10000000:10000010]
chr19_SeqRevComp
length(chr19_Seq[10000000:10000008])
chr19_SeqTranslation <- translate(chr19_Seq[10000000:10000008])
chr19_SeqTranslation
chr19_Count <- countPattern(pattern="ATCTGCAATG",
                            subject=chr19_Seq)
chr19_Count
chr19_Count <- countPattern(pattern="ATCTGCAATG",
                            subject=chr19_Seq,
                            max.mismatch = 2,
                            min.mismatch = 0)
chr19_Count

IUPAC_CODE_MAP
chr19_Count <- countPattern(pattern="RYKHBNKYSRR",
                            subject=chr19_Seq,
                            fixed=FALSE)
chr19_Count
chr19_Count_Watson <- countPattern(pattern="ATCTGCAATG", subject=chr19_Seq)
chr19_Count_Crick <- countPattern(pattern="ATCTGCAATG", subject=reverseComplement(chr19_Seq))
Total_chr19_Count <- chr19_Count_Watson+chr19_Count_Crick

chr19_SeqSet <- DNAStringSet(chr19_Seq[10000000:10000010])
names(chr19_SeqSet) <- "chr19"
writeXStringSet(chr19_SeqSet,filepath = "chr19_Seq.fa")

chr19_FromFASTA <- readDNAStringSet(filepath = "chr19_Seq.fa")
chr19_FromFASTA$chr19

library(BSgenome.Hsapiens.UCSC.hg19)
length(seqnames(BSgenome.Hsapiens.UCSC.hg19)) # Count the number of contigs

contigLengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
sum(contigLengths[order(contigLengths,
                        decreasing = FALSE)][1:3]) # Give the sum of lengths of the 3 smallest chromosomes

alpFreq <-alphabetFrequency(BSgenome.Hsapiens.UCSC.hg19$chr20)
alpFreq["N"] # How many unknown bases - base N - are in chromosome 20

library(ggplot2)
alpFreq <-alphabetFrequency(BSgenome.Hsapiens.UCSC.hg19$chr20)
atcgFreq <- alpFreq[c("A","T","C","G")]
atcgFreqDF <- data.frame(Bases=names(atcgFreq),Frequency=atcgFreq)
ggplot(atcgFreqDF, # Create a barchart of the total number of the A,T,C,G bases on chromosome 20
       aes(y=Frequency,x=Bases,fill=Bases)) +
  geom_bar(stat = "identity") +
  theme_minimal()

mySeq <- BSgenome.Hsapiens.UCSC.hg19$chr20[1000000:1000020]
complement(mySeq) # Extract the sequence from chromosome 20 at position 1,000,000 to 1,000,020 and retrieve the complement sequence

mySet <- DNAStringSet(complement(mySeq))
names(mySet) <- "testFasta"
writeXStringSet(mySet,filepath ="testFasta.fa") # Write this complement sequence to a FASTA file

# Look up the position of MYC in IGV (Human hg19) and find the genomic coordinates of its first exon
mySeq <- BSgenome.Hsapiens.UCSC.hg19$
  chr8[128748315:128748869] # Extract the sequence for the first exon
translate(mySeq[1:3]) # Compare the sequence to that found in IGV and identify start of translated region in gene
countPattern("ATG",mySeq) # Count the number of classical start codons (ATG) in the first exon

gapdhSeq <- BSgenome.Hsapiens.UCSC.hg19$chr12[6643976:6644027]
translate(gapdhSeq[1:3])
countPattern("ATG",gapdhSeq) # count occurrence of ATG in exon 2 of NM_001289745 transcript. (chr12:6643976-6644027)



library(GenomicRanges)

myIntervals <- IRanges(start=c(1,10,20),end=c(2,11,30))
class(myIntervals)
myIntervals
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals)
class(myGenomicIntervals)
myGenomicIntervals
myGenomicIntervals[c(2,3),]

start(myGenomicIntervals)
start(myGenomicIntervals) <- c(1,5,15)
start(myGenomicIntervals)
start(myGenomicIntervals) <- c(10,50,150)
contigNames <- seqnames(myGenomicIntervals)
contigNames
as.character(contigNames)
seqnames(myGenomicIntervals) <- c("chr2","chr2","chr1")
seqnames(myGenomicIntervals)
seqnames(myGenomicIntervals) <- c("chr3","chr2","chr1")

seqlevels(myGenomicIntervals) <- c("chr1","chr2","chr3")
seqlevels(myGenomicIntervals)
seqnames(myGenomicIntervals) <- c("chr1","chr2","chr3")
myGenomicIntervals

myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"), 
                              myIntervals,
                              strand=c("+","-","*"))
myGenomicIntervals
strand(myGenomicIntervals) <- c("+","+","-")
strand(myGenomicIntervals)
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals,
                              strand=c("+","+","+"))
strand(myGenomicIntervals)

names(myGenomicIntervals) <- c("peak1","peak2","peak3")
myGenomicIntervals
myGenomicIntervals["peak2"]
myIntervals <- IRanges(start=c(1,2,2),end=c(2,11,30))
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              myIntervals,strand=c("+","+","+"),
                              Score=c(10,20,40),
                              Comment=c("GoodQC","GoodQC","BadQC"))
myGenomicIntervals$Score
myGenomicIntervals$Comment
mcols(myGenomicIntervals)
as.data.frame(myGenomicIntervals)
GRanges("chr1:110-120")

library(stringi)
myRange <- "chr1:110:120"
newRange <- stri_replace_last_fixed(myRange, ':', '-')
newRange

myRange <- "MyID:chr1:110:120"
newRange <- stri_replace_last_fixed(myRange, ':', '-')
newerRange <- stri_replace_first_regex(newRange, '\\w*:', '')
newerRange
GRanges(newerRange)

myGenomicIntervals[1]
shift(myGenomicIntervals[1],shift = 10)
myGenomicIntervals[3]
resize(myGenomicIntervals[3], width=20)
myGenomicIntervals[3]
resize(myGenomicIntervals[3], width=20, fix = "end")
resize(myGenomicIntervals[3], width=20, fix = "start")
strand(myGenomicIntervals)[3] <- "-"
resize(myGenomicIntervals[3],width=20, fix="start")
myGenomicIntervals <- GRanges(seqnames=c("chr1","chr1","chr2"),
                              ranges=IRanges(start=c(1,2,2),end=c(2,11,30)),
                              strand=c("+","+","+"))
myGenomicIntervals
reduce(myGenomicIntervals)
strand(myGenomicIntervals) <- c("+","-","+")
reduce(myGenomicIntervals)
strand(myGenomicIntervals) <- c("+","-","+")
reduce(myGenomicIntervals, ignore.strand=TRUE)
myGenomicIntervals1 <- GRanges(seqnames=c("chr1","chr1"),
                               ranges=IRanges(start=c(1,25),
                                              end=c(20,30)),
                               strand=c("+","+"))

myGenomicIntervals2 <- GRanges(seqnames=c("chr1","chr1"),
                               ranges=IRanges(start=c(22,100),
                                              end=c(27,130)),
                               strand=c("+","+"))
myGenomicIntervals1 %over% myGenomicIntervals2
myGenomicIntervals1[.Last.value]

myOverlaps <- findOverlaps(myGenomicIntervals1,myGenomicIntervals2)
class(myOverlaps)
myOverlaps
queryHits(myOverlaps)
myGenomicIntervals1[.Last.value]
subjectHits(myOverlaps)
myGenomicIntervals2[.Last.value]

myGenomicIntervals1 <- GRanges(seqnames=c("chr1","chr1"),
                               ranges=IRanges(start=c(10,20),
                                              end=c(25,30)),
                               strand=c("+","+"))

myGenomicIntervals2 <- GRanges(seqnames=c("chr1","chr1"),
                               ranges=IRanges(start=c(1,10000),
                                              end=c(2,10002)),
                               strand=c("+","+"))
indexOfNearest <- nearest(myGenomicIntervals1,myGenomicIntervals2)
indexOfNearest
myGenomicIntervals2[indexOfNearest]
precedeIndex <- precede(myGenomicIntervals1,myGenomicIntervals2)
followIndex <- follow(myGenomicIntervals1,myGenomicIntervals2)
myGenomicIntervals2[precedeIndex]
distances <- distanceToNearest(myGenomicIntervals1,myGenomicIntervals2)
distances

library(rtracklayer)
mySicerPeaks <- import.bed("SicerPeaks.bed")
mySicerPeaks
export.bed(mySicerPeaks, "moreSicerPeaks.bed")

library(BSgenome.Mmusculus.UCSC.mm10)
subseq(BSgenome.Mmusculus.UCSC.mm10$chr10,1,100)
sicerSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, names=mySicerPeaks)
sicerSeq
writeXStringSet(sicerSeq,"sicerSeq.fa")

# Read in the file data/Myc_Ch12_1_withInput_Input_Ch12_peaks.xls and create a GRanges object which includes values of fold_enrichment and - log10 pvalue provided in file.

MYCpeaks <- read.delim("Myc_Ch12_1_withInput_Input_Ch12_peaks.xls",sep="\t",comment.char = "#")
MYCgranges <- GRanges(MYCpeaks$chr,IRanges(MYCpeaks$start,MYCpeaks$end),
                      FE=MYCpeaks$fold_enrichment,
                      minuslog10Pval=MYCpeaks$X.log10.pvalue.)
# Create a boxplot of the fold enrichments in genomic intervals over every chromosome
MYCDF <- as.data.frame(MYCgranges)
library(ggplot2)
ggplot(MYCDF, # Create a boxplot of the fold enrichments in genomic intervals over every chromosome
       aes(x=seqnames,y=FE,color=seqnames)) +
  geom_boxplot() +
  theme_bw()

# Create a GRanges object of genomic intervals on chromosome 1 with scores greater than 10 and pvalue less than 0.0001 and export as BED file filteredMyc.bed
MYCfilteredGRanges <- MYCgranges[MYCgranges$FE > 10 & MYCgranges$minuslog10Pval > -log10(0.0001)]
MYCfilteredGRangesChr1 <- MYCfilteredGRanges[seqnames(MYCfilteredGRanges) %in% "chr1"]
export.bed(MYCfilteredGRangesChr1,"filteredMyc.bed")

# Read in TXT file of gene positions (containing contig - named as seqnames column-, genomic start and end) from the file data/mm10_GenePosForIGV.txt and export to a BED file Genes.bed
IGVGenePositions <- read.delim("mm10_GenePosForIGV.txt",sep="\t")
genePos <- GRanges(IGVGenePositions$seqnames,
                   IRanges(IGVGenePositions$start,IGVGenePositions$end),
                   strand = IGVGenePositions$strand)
names(genePos) <- IGVGenePositions$gene_id
export.bed(genePos,con="Genes.bed")

# Create a GRanges of the transcriptional start site positions of every gene (1bp exact TSS).
tssPos <- resize(genePos,width = 1,fix = "start")
# Extend this GRanges to be +/- 500 bp around the transciptional start sites.
tssPosExt <- resize(tssPos,width = 1000,fix = "center")
# Create a BED file (called filteredMycOnTSS.bed) of our Myc peaks from question 4 which overlap our new GRanges of +/- 500 bp around transciptional start sites.
MYCforBED <- MYCfilteredGRangesChr1[MYCfilteredGRangesChr1 %over% tssPosExt]
export.bed(MYCforBED,con="filteredMycOnTSS.bed")

# Import the Myc_Ch12_1_withInput_Input_Ch12_summits.bed BED file to a GRanges. The 5th column in file represents summit height
mySummits <- import.bed("Myc_Ch12_1_withInput_Input_Ch12_summits.bed")
mySummits$Overlap <- ifelse(mySummits %over% tssPosExt,"TSS","Not_TSS")
mySummitsDF <- as.data.frame(mySummits)
ggplot(mySummitsDF,aes(x=score,fill=Overlap)) +
  geom_density(alpha=0.5) +
  facet_wrap(~seqnames) +
  theme_bw()
# Filter the Summits GRanges to the top 500 ranked by the GRanges score column.
orderSummits500 <- mySummits[order(mySummits$score,decreasing = TRUE)][1:500,]
# Extend these top 500 Summits GRanges to 50bps around the centre of the GRanges intervals. Extract the sequences under the peaks and write to a file.
centredSummits <- resize(orderSummits500,width = 50,fix="center")
library(BSgenome.Mmusculus.UCSC.mm10)
centredSummitsSeq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,centredSummits)
writeXStringSet(centredSummitsSeq,filepath = "top500bySummit.fa")





library(rtracklayer)
myBedG <- import.bedGraph("TSS_ENCFF940MBK.bedGraph")
class(myBedG)
myBedG[1:3]
strand(myBedG)
myBigWig <- import.bw("TSS_ENCFF940MBK.bw") # ENCFF940MBK.bigWig
class(myBigWig)
myBigWig[1:3]

myGRanges <- GRanges("chr1", IRanges(72823698,72824485))
filteredBigWig <- myBigWig[myBigWig %over% myGRanges]
filteredBigWig[1:3]
myBigWig <- import.bw("TSS_ENCFF940MBK.bw", # ENCFF940MBK.bigWig
                      as = "RleList")
class(myBigWig)
myBigWig[1:2]
chr1_rle <- myBigWig$chr1
chr1_rle # or
chr1_rle <- myBigWig[["chr1"]]
chr1_rle
chr1_rle[1:10]
chr1_rle[1:10] <- 100
chr1_rle
rleAsVector <- as.vector(chr1_rle[1:10])
rleAsVector
rleAsDF <- as.data.frame(chr1_rle[1:10])
rleAsDF
chr1_rle+1000
(chr1_rle+1000)*10
chr1_rle < 10
chr1_rle[chr1_rle < 10] <- 0
chr1_rle
mean(chr1_rle)
max(chr1_rle)
sum(chr1_rle)

myBigWig <- import.bw("TSS_ENCFF940MBK.bw", # ENCFF940MBK.bigWig
                      as = "RleList")
myBigWig + 10

chromosomeMax <- max(myBigWig)
chromosomeMax[1:4]

myRanges <- GRanges("chr1", IRanges(72811055,72856974))
mycPeaks <- import.bed("Myc_Ch12_1_withInput_Input_Ch12_summits.bed")
mycPeaks <- resize(mycPeaks,50,fix="center")
newMycPeaks <- mycPeaks[mycPeaks %over% myRanges]
newMycPeaks

rleOverGranges <- myBigWig[newMycPeaks] # for RLE
rleOverGranges # for RLE
sum(rleOverGranges)

myRleList <- RleList(chr1=chr1_rle)
myRleList
export.bw(myRleList,con="chr1_Myc.bigWig")

newMycPeaks

mySelection <- BigWigSelection(newMycPeaks)
import.bw("TSS_ENCFF940MBK.bw", 
          selection=mySelection, 
          as="RleList")
# Read in the bigWig signal for the region (chr1:133040000-133149400) as an RleList
regionToImport <- GRanges(seqnames="chr1",
                          ranges=IRanges(133040000,133149400))
mySelection <- BigWigSelection(regionToImport)
bigWigRegion <- import.bw("ENCFF940MBK.bigWig",
                          selection=mySelection,
                          as="RleList")
# Find the maximum and minimum of value of bigWig signal in the region
myMax <- max(bigWigRegion$chr1[133040000:133149400])
myMax <- max(bigWigRegion)["chr1"]

# Produce an area plot (using geom_area) of the signal over the transcriptional start site (+/- 500bp of gene start) for the Ppp1r15b gene using ggplot2
TSSpos <- seq(133131166-500,133131166+500)
tssDF <- as.data.frame(bigWigRegion$chr1[TSSpos])
tssDF$position <- TSSpos
library(ggplot2)
ggplot(tssDF,aes(y=value,x=position)) +
  geom_area()+theme_minimal()+xlab("Genome")+ylab("Score")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),axis.title.y = element_text(angle = 0))
# Apply a log2 transformation to the genomic scores and save as a new variable
newBigWigRegion <- log2(bigWigRegion)
export.bw(newBigWigRegion,"log2.bigWig")
# Read in the file Myc_Ch12_1_withInput_Input_Ch12_peaks.xls and create a GRanges object of all peaks on chromosome 1
MYCpeaks <- read.delim("Myc_Ch12_1_withInput_Input_Ch12_peaks.xls",sep="\t",comment.char = "#")
MYCgranges <- GRanges(MYCpeaks$chr,IRanges(MYCpeaks$start,MYCpeaks$end),
                      FE=MYCpeaks$fold_enrichment,
                      minuslog10Pval=MYCpeaks$X.log10.pvalue.)
MYCgrangesChr1 <- MYCgranges[seqnames(MYCgranges) %in% "chr1"]
# Using this GRanges of Myc peaks, import in the Myc ChIPseq bigwig signal over the peaks on chromosome 1
mySelection <- BigWigSelection(MYCgrangesChr1)
bigWigPeaks <- import.bw("ENCFF940MBK.bigWig",
                         selection=mySelection,
                         as="RleList")
# Import the gene positions from mm10_GenePosForIGV.txt as a GRanges. Filter to genes on chromosome 1 and 
# identify peaks overlapping the gene’s promoter (here defined as a region ranging between 2000bp upstream from 
# a gene’s transcriptional start site to the transcriptional start site) and those not overlapping a promoter
IGVGenePositions <- read.delim("mm10_GenePosForIGV.txt",sep="\t")
genePos <- GRanges(IGVGenePositions$seqnames,
                   IRanges(IGVGenePositions$start,IGVGenePositions$end),
                   strand = IGVGenePositions$strand)
names(genePos) <- IGVGenePositions$gene_id
tssPos <- resize(genePos,width = 1,fix = "start")
promoterPos <- resize(tssPos,width = 2000,fix = "end")

promoterPeaks <- MYCgrangesChr1[MYCgrangesChr1 %over% promoterPos]
nonPromoterPeaks <- MYCgrangesChr1[!MYCgrangesChr1 %over% promoterPos]
# Create a violin plot of sum of bigWig signal within peaks at promoters and those not at promoters using ggplot2
PeakCoveragePromoter <- bigWigPeaks[promoterPeaks]
PeakMaxs_Promoter <- sum(PeakCoveragePromoter)

PeakCoverageNonPromoter <- bigWigPeaks[nonPromoterPeaks]
PeakMaxs_NonPromoter <- sum(PeakCoverageNonPromoter)

df_promoter <- data.frame(Sum=PeakMaxs_Promoter,PromoterOrNot="Promoter")
df_notPromoter <- data.frame(Sum=PeakMaxs_NonPromoter,PromoterOrNot="NotPromoter")
toPlot_SumSignal <- rbind(df_promoter,df_notPromoter)

ggplot(toPlot_SumSignal,aes(x=PromoterOrNot,y=Sum,fill=PromoterOrNot))+
  geom_violin()+scale_y_log10()+scale_fill_brewer(palette = "Pastel2")+theme_minimal()
# Export a bed file of TSS peaks whose sum of ChIPseq signal is 2 fold greater than the median of all peaks’ summed signals
toExport <- bigWigPeaks
peakSumMedian <- median(c(sum(PeakCoveragePromoter),sum(PeakCoverageNonPromoter)))
filteredPeaks <- promoterPeaks[PeakMaxs_Promoter > peakSumMedian*2]
export.bed(filteredPeaks,"filteredPeaks.bedGraph")



library(TxDb.Hsapiens.UCSC.hg19.knownGene)
class(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb.Hsapiens.UCSC.hg19.knownGene
myGenes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
myGenes
myExons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
myExons
myTranscripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
myTranscripts
myCDS <- cds(TxDb.Hsapiens.UCSC.hg19.knownGene)
myCDS
myTranscripts <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             columns=c("gene_id","tx_id"))
myTranscripts[1:2]
myPromoters <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                         upstream=2000,downstream=50)
myPromoters[1:2]
transcriptByGenes <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   by="gene")
transcriptByGenes[1:2]
transcriptByGenes$`1` # or
transcriptByGenes[[1]]
exonsByTranscript <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                             by="tx")
exonsByTranscript[1:2]
myCustomTxDb <- makeTxDbFromGFF("Xkr4.gtf")
class(myCustomTxDb)
myCustomTxDb
genes(myCustomTxDb)
exonsBy(myCustomTxDb,by="gene")
availableGenomes <- ucscGenomes()
availableGenomes[1:4,]

library(RMariaDB)
# hg18TxDb <- makeTxDbFromUCSC(genome="hg18")
hg18TxDb
class(hg18TxDb)
hg18Promoters <- promoters(hg18TxDb,upstream=2000,downstream=50)
export.gff(myCustomTxDb,"customTxDbb.gff",format="gff")
export.gff(myCustomTxDb,"customTxDbb.gtf",format="gtf")

length(exonsBy(hg18TxDb,by="gene"))
lengths(exonsBy(hg18TxDb,by="gene"))[1:5]
transcriptLengths(hg18TxDb)[1:5,]

library(BSgenome.Hsapiens.UCSC.hg19)
hg19TransSeq <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg19, 
                                      transcripts=TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19TransSeq
writeXStringSet(hg19TransSeq,"myTranscriptSequences.fa")

library(GenomeInfoDb)
allMappings <- genomeStyles()
names(allMappings)
allMappings[["Homo_sapiens"]]
# myGenes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
seqlevelsStyle(myGenes)
myGenes[1:2,]
seqlevelsStyle(myGenes) <- "Ensembl"
myGenes[1:2,]

library(org.Hs.eg.db)
class(org.Hs.eg.db)
columns(org.Hs.eg.db)
help(GENENAME)
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db, keytype="SYMBOL")[1:10]
select(org.Hs.eg.db, keys = "A1BG", keytype = "SYMBOL", 
       columns = c("SYMBOL", "GENENAME", "CHR"))

geneLocations <- genes(hg18TxDb)
geneLocations

IDs <- geneLocations$gene_id
myTable <- select(org.Hs.eg.db, keys = IDs, keytype = "ENTREZID",
                  columns = c("SYMBOL", "GENENAME", "ENTREZID") )
myTable[1:2,]

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# Count the number of genes and transcripts
length(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))
length(transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene))
# Plot the distribution of log10 gene lengths (as on genome including introns) across selected chromosomes (chr1, chr2, chr3) as density plots using ggplot2
geneTable <- as.data.frame(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))
library(ggplot2)
ggplot(geneTable[geneTable$seqnames %in% c("chr1","chr2","chr3"),],
       aes(x=width,fill=seqnames))+geom_density()+scale_x_log10()+
  facet_grid(seqnames~.)+theme_minimal()+scale_fill_discrete(name="Chromosome")
# Retrieve the log10 transcript lengths (sum of exons in transcripts) and plot against the total log10 transcript lengths (as on genome including introns) for transcripts on chromosome 1 using ggplot2. Scale the size of points by the number of exons in a transcript
allTX <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
transcriptTable <- as.data.frame(allTX)
filtTranscriptTable <- transcriptTable[transcriptTable$seqnames %in% c("chr1"),]
txLengths <- transcriptLengths(TxDb.Mmusculus.UCSC.mm10.knownGene)
toPlotTxLen <- merge(filtTranscriptTable,txLengths,all=FALSE,by="tx_name")
ggplot(toPlotTxLen,aes(x=width,y=tx_len,size=nexon))+
  geom_point()+scale_x_log10()+scale_y_log10()+theme_minimal()
# Plot the nucleotide content for the longest transcript (by genomic size including introns) overlapping the region (chr1:3,083,473-3,876,510.) using ggplot2
filtTX <- allTX[allTX %over% GRanges("chr1",IRanges(3083473,3876510))]
longestTX <- filtTX[order(width(filtTX),decreasing = TRUE),][1]
allSeqs <- extractTranscriptSeqs(BSgenome.Mmusculus.UCSC.mm10,TxDb.Mmusculus.UCSC.mm10.knownGene)
alpFreq <-alphabetFrequency(allSeqs[longestTX$tx_id])
atcgFreq <- alpFreq[,c("A","T","C","G")]
atcgFreqDF <- data.frame(Bases=names(atcgFreq),Frequency=atcgFreq)
ggplot(atcgFreqDF,aes(y=Frequency,x=Bases,fill=Bases))+
  geom_bar(stat = "identity")+theme_minimal()
# In the GM12878_minus_HeLa.csv file, Column 1 contains Entrez IDs and column 2 contains symbols. Add a column of gene names to the table
myDETable <- read.delim("GM12878_minus_HeLa.csv",sep=",",stringsAsFactors = FALSE)
myAnnotationTable <- select(org.Hs.eg.db,keys=as.character(myDETable$ID),keytype = "ENTREZID",columns=c("ENTREZID","GENENAME"))
newTable <- merge(myAnnotationTable,myDETable,by=1,all.x=FALSE,all.y=TRUE)
newTable <- newTable[order(newTable$padj),]
# Identify genes with a padj < 0.05 and logFC greater than 1. Export promoter positions (upstream 2000bp, downstream 50bp) of these genes’s transcripts, merge overlapping promoters and export to a bed file and visualise in IGV.
genesToSelect <- newTable$ENTREZID[newTable$padj < 0.05 & !is.na(newTable$padj) & newTable$log2FoldChange > 1]
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
myPromoters <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene,upstream = 2000,downstream = 50,columns=c("gene_id"))
selectedPromoters <- myPromoters[as.vector(myPromoters$gene_id) %in% genesToSelect]
reducedselectedPromoters <- reduce(selectedPromoters)
export.bed(reducedselectedPromoters,con="DE_Promoters.bed")

library(ShortRead)
fastQ <- readFastq("sampled_ENCFF000CXH.fastq.gz")
class(fastQ)
fastQ
length(fastQ)
readLengths <- width(fastQ)
readLengths[1:10]
fastQ[1:10]
sequenceOfReads <- sread(fastQ)
class(sequenceOfReads)
sequenceOfReads
alpFreq <- alphabetFrequency(sequenceOfReads)
alpFreq[1:2,]
idsOfReads <- id(fastQ)
class(idsOfReads)
idsOfReads[1:2]
Ids <- as.character(idsOfReads)
Ids[1:4]
quals <- quality(fastQ)
class(quals)
quals
qualityEncoding <- encoding(quals)
qualityEncoding
quals[1]
quals[[1]]
toTranslateList <- strsplit(as.character(quals[[1]]),"")
toTranslate <- unlist(toTranslateList)
toTranslate
qualityEncoding[toTranslate]
readScores<- alphabetScore(quals)
readScores[1]
sum(qualityEncoding[toTranslate])
matrixOfQualities <- as(quals,"matrix")
rowSums(matrixOfQualities)[1]
alpByCyc <- alphabetByCycle(sequenceOfReads)
alpByCyc[1:4,1:15]
qualsByCyc <- alphabetByCycle(quals)
qualsByCyc[1:4,1:15]
readOccurence <- table(sequenceOfReads)
sort(readOccurence,decreasing = TRUE)[1:2]
duplicates <- srduplicated(fastQ)
duplicates[1:3]
table(duplicates)
my_QA <- qa("sampled_ENCFF000CXH.fastq.gz")
my_QA
myReport <- report(my_QA)
myReport
browseURL(myReport)
myReport <- report(my_QA, dest="QC_report")
browseURL(myReport)
TrimmedFastq <- trimTails(fastQ,20,"5")
TrimmedFastq
writeFastq(TrimmedFastq,"myTrimmed_Fastq.fastq.gz")

library(Rfastp)
rfastp_report <- rfastp(read1 = "sampled_ENCFF000CXH.fastq.gz", outputFastq ="Rfastp_ENCFF000CXH")
dfsummary <- qcSummary(rfastp_report)
dfsummary
curvePlot(rfastp_report, curve="quality_curves")
curvePlot(rfastp_report, curve="content_curves")
sampleToRead <- FastqSampler("sampled_ENCFF000CXH.fastq.gz",
                             n=100)
yield(sampleToRead)
sampleToRead <- FastqStreamer("sampled_ENCFF000CXH.fastq.gz",
                              n=100)
first100Reads <- yield(sampleToRead)
second100Reads <- yield(sampleToRead)
fq <- FastqStreamer("sampled_ENCFF000CXH.fastq.gz",
                    n=25000)
while (length(fq_stream <- yield(fq)) > 0) {
  print(length(fq_stream ))
}

# Plot the occurrence of A, G, T, C and N in the read with ID “HWI-BRUNOP16X_0001:7:66:3246:157000#0/1” using ggplot2
library(ShortRead)
my_fastq <- readFastq("heart.bodyMap.fq")
fastqRead <- my_fastq[as.character(id(my_fastq)) == "HWI-BRUNOP16X_0001:7:66:3246:157000#0/1"]
fastqReadSeq <- sread(fastqRead)
alpFreq <- alphabetFrequency(fastqReadSeq)
alpFreq5 <- alpFreq[,c("A","C","G","T","N")]
toPlot <- data.frame(bases=names(alpFreq5),Freq=alpFreq5)
library(ggplot2)
ggplot(toPlot,aes(x=bases,y=Freq,fill=bases))+geom_bar(stat="identity") + theme_bw()

# Create a boxplot of quality scores over cycles for first 10,000 reads using base graphics
allQuals <- quality(my_fastq[1:10000,])
forBox <- as(allQuals,"matrix")
colnames(forBox) <- paste0("Cycle",1:ncol(forBox))
boxplot(forBox)

# Plot the frequency of A, G, C, T and N bases over cycles using ggplot2
alpByCyle <- alphabetByCycle(sread(my_fastq))
alpByCyleFilt <-  alpByCyle[c("A","G","C","T","N"),]
AbyCycFrame <- data.frame(Base="A",Freq=alpByCyleFilt["A",],Cycle=1:75)
CbyCycFrame <- data.frame(Base="C",Freq=alpByCyleFilt["C",],Cycle=1:75)
TbyCycFrame <- data.frame(Base="T",Freq=alpByCyleFilt["T",],Cycle=1:75)
GbyCycFrame <- data.frame(Base="G",Freq=alpByCyleFilt["G",],Cycle=1:75)
NbyCycFrame <- data.frame(Base="N",Freq=alpByCyleFilt["N",],Cycle=1:75)
myFrame <- rbind(AbyCycFrame,CbyCycFrame,TbyCycFrame,GbyCycFrame,NbyCycFrame)
ggplot(myFrame,aes(x=Cycle,y=Freq,colour=Base))+geom_line()+theme_bw()

# How many reads are in the original fastq? Read in a random 10000 bases from the original fastq
mySampledFile <- FastqSampler("heart.bodyMap.fq",n=10000)
my_fastq_10000 <- yield(mySampledFile)
my_fastq_10000

# Count the number of duplicate sequences and filter file to only the non-duplicated sequences
dupLogical <- srduplicated(my_fastq_10000)
table(dupLogical)
nonDups <- my_fastq_10000[dupLogical]

# Trim the resulting non-duplicated reads to remove reads when a succession of 10 bases fall below quality of 5
trimmedNonDups <- trimTails(nonDups,10,"5")

# Plot a histogram of the resulting read lengths using base graphics
hist(width(trimmedNonDups),xlab="Trimmed Length",main="Histogram of trimmed lengths")

# Filter out reads less than 32 bases and write to file
filteredFastq <- trimmedNonDups[width(trimmedNonDups) >= 32]
writeFastq(filteredFastq,file="filteredFastq.fq.gz")
unlink("filteredFastq.fq.gz")

library(Rsubread)
library(BSgenome.Hsapiens.UCSC.hg38)
chr7hg38 <- BSgenome.Hsapiens.UCSC.hg38[["chr7"]]
chr7hg38Set <- DNAStringSet(list(chr7=chr7hg38))
writeXStringSet(chr7hg38Set,file="chr7.fa")
buildindex("chr7","chr7.fa", memory=8000)
align("chr7","sampledActin.fq.gz",
      output_format = "BAM",
      output_file = "Rsubread_NoSplicing_sampledActin.bam")
library(Rsamtools)
sortBam("Rsubread_NoSplicing_sampledActin.bam","SortedActB")
indexBam("SortedActB.bam")
quickBamFlagSummary("SortedActB.bam")
subjunc("chr7","sampledActin.fq.gz",
        output_format = "BAM",
        output_file = "RsubreadsampledActin.bam")
sortBam("RsubreadsampledActin.bam",
        "SortedActBSpliced")
indexBam("SortedActBSpliced.bam")
quickBamFlagSummary("SortedActBSpliced.bam")

library(Rbowtie2)
bowtie2_build(references="chr7.fa",
              bt2Index=file.path("chr7hg38"))
# R.utils::gunzip("sampledActin.fq.gz")
bowtie2(bt2Index = "chr7hg38",
        samOutput = "sampledActin.sam",
        seq1 = "sampledActin.fq")
bamFile_Bowtie2 <- asBam("sampledActin.sam")
bamFile_Bowtie2
sortBam(bamFile_Bowtie2,"SortedActBSpliced_bowtie")
indexBam("SortedActBSpliced_bowtie.bam")
library(QuasR)
FileName <- "sampledActin.fq.gz"
SampleName <- "sampledActin"
sampleTable <- data.frame(FileName,SampleName)
write.table(sampleTable,file="sampleTable.txt",sep="\t",quote=FALSE,row.names = FALSE)
sampleTable
qAlign("sampleTable.txt", "BSgenome.Hsapiens.UCSC.hg38")
qAlign("sampleTable.txt","chr7.fa")
qAlign("sampleTable.txt","chr7.fa", aligner="Rhisat2")

library(Rsamtools)
coordSorted <- sortBam("liver.bodyMap.bam",
                       "Sorted_liver")
coordSorted
readnameSorted <- sortBam("liver.bodyMap.bam",
                          "SortedByName_liver",
                          byQname=TRUE)
readnameSorted
coordSorted <- sortBam("liver.bodyMap.bam",
                       "Sorted_liver",
                       maxMemory=1)
coordSorted
indexBam("Sorted_liver.bam")
quickBamFlagSummary("Sorted_liver.bam")
idxstatsBam("Sorted_liver.bam")
library(GenomicAlignments)
myHeader <- scanBamHeader("Sorted_liver.bam")
str(myHeader)
names(myHeader)
names(myHeader$Sorted_liver.bam)
myHeader$Sorted_liver.bam$targets
myHeader$Sorted_liver.bam$text
myHeader$Sorted_liver.bam$text["@HD"]
myHeader <- scanBamHeader("SortedByName_liver.bam")
myHeader$SortedByName_liver.bam$text["@HD"]
myHeader <- scanBamHeader("Sorted_liver.bam")
myHeader$Sorted_liver.bam$text["@PG"]
myReads <- readGAlignments("Sorted_liver.bam")
class(myReads)
myReads[1:2,]
seqnames(myReads)
cigar(myReads)[1:2]
njunc(myReads)[1:2]
myReads[strand(myReads) == "+"]
my5primeReads <- narrow(myReads, start=1, width = 1)
my5primeReads[1:2]
myReadsPos <- narrow(myReads[strand(myReads) == "+"],
                     start=1, width = 1)
myReadsNeg <- narrow(myReads[strand(myReads) == "-"],
                     end=-1, width = 1)
my5primeReads <- c(myReadsPos,myReadsNeg)
my5primeReads[1:2]
myReadAsGRanges <- granges(myReads,use.mcols = TRUE)
myReadAsGRanges
myReadAsGRangesList <- grglist(myReads,use.mcols = TRUE)
myReadAsGRangesList[njunc(myReads) == 1]
myReadAsGRanges <- granges(myReads, use.mcols = TRUE)
myReadsAgain <- as(myReadAsGRanges, "GAlignments")
myReadsAgain[1:2]
myReadAsGRanges <- granges(myReads, use.mcols = TRUE)
my5Prime <- resize(myReadAsGRanges, fix = "start", width = 1)
my5PrimeAsReads <- as(my5Prime, "GAlignments")
my5PrimeAsReads
library(rtracklayer)
export(my5PrimeAsReads, con="myModifiedReads.bam")
myRanges <- GRanges("chr12", IRanges(98591400,98608400))
myParam <- ScanBamParam(which=myRanges)
myParam
filteredReads <- readGAlignments("Sorted_liver.bam", param = myParam)
filteredReads
myParam <- ScanBamParam(what=c("qname", "seq", "qual"))
infoInReads <- readGAlignments("Sorted_liver.bam", param = myParam)
infoInReads[1]
mcols(infoInReads)
bamHeader <- scanBamHeader("Sorted_liver.bam")
myChromosomes <- bamHeader$Sorted_liver.bam$targets
for(i in 1:length(myChromosomes)){
  grangesForImport <- GRanges(names(myChromosomes)[i],
                              IRanges(1,myChromosomes)[i])
  myParam <- ScanBamParam(which = grangesForImport)
  myReads <- readGAlignments("Sorted_liver.bam", 
                             param=myParam)
  print(length(myReads))
}
# Coordinate sort and index the aligned reads in BAM files Heart.bam and Liver.bam
sortedHeart <- sortBam("heart.bodyMap.bam","SortedHeart")
indexBam(sortedHeart)
sortedLiver <- sortBam("liver.bodyMap.bam","SortedLiver")
indexBam(sortedLiver)
# Plot the number of mapped reads on every chromsome in the Heart and Liver BAM files using ggplot2
library(ggplot2)
idxHeart <- idxstatsBam(sortedHeart)
idxLiver <- idxstatsBam(sortedLiver)
idxHeartDF <- data.frame(Sample="Heart",idxHeart)
idxLiverDF <- data.frame(Sample="Liver",idxLiver)
toPlot <- rbind(idxHeartDF,idxLiverDF)
ggplot(toPlot,aes(x=seqnames,y=mapped,fill=seqnames))+
  geom_bar(stat="identity")+
  facet_wrap(~Sample)+
  coord_flip()+
  theme_bw()+xlab("Chromosome")
# Using the qwidth() and the width() functions,
# plot the length of reads vs the length of their alignment for the Heart bam file using ggplot2.
# Facet the plot by the number of junctions a read spans.
myReads <- readGAlignments("heart.bodyMap.bam")
toPlot <- data.frame(readLength=qwidth(myReads),alignmentLength=width(myReads),junctions=factor(njunc(myReads)))
ggplot(toPlot,aes(x=readLength,y=alignmentLength))+
  geom_point()+facet_grid(~junctions)+
  theme_minimal()+xlab("Read Length")+ylab("Alignment Length")
# Plot this again, but add a limit to the alignment of 2.5Kb in the plotting
ggplot(toPlot,aes(x=readLength,y=alignmentLength))+
  geom_point()+facet_grid(~junctions)+
  theme_minimal()+xlab("Read Length")+ylab("Alignment Length")+ylim(0,2500)
# Export any aligned reads spanning more than 40000 bp on the genome to a BAM file and review in IGV
export(myReads[width(myReads) > 40000],"longAlignments.bam")
# Import the read IDs, sequence and qualities from the Heart BAM file
myParam <- ScanBamParam(what=c("qname","seq","qual"))
infoInReads <- readGAlignments("SortedHeart.bam",param = myParam)
# Find the number of unique read IDs and compare to total reads in file
readnames <- mcols(infoInReads)$qname
uniqueIDs <- length(unique(readnames))
uniqueIDs
totalReads <- length(infoInReads)
totalReads
# Plot the A,G,C,T,N content of 75bp reads in file
uniqueReads <- infoInReads[!duplicated(readnames) & qwidth(infoInReads) == 75]
seqOfReads <- mcols(uniqueReads)$seq
alpFreq <- alphabetFrequency(seqOfReads)
sumedAlpFreq <- colSums(alpFreq)
mainBases <- sumedAlpFreq[c("A","C","G","T","N")]
toPlot <- data.frame(Base=names(mainBases),Freq=mainBases)
ggplot(toPlot,aes(x=Base,y=Freq,fill=Base))+geom_bar(stat="identity")+theme_bw()
# Using a loop and ScanBamParam, count the number number reads in Heart file 
# overlapping each exon for the SLC25A3 gene.
# Remember we can use TxDb objects to extract GRanges of exon positions for genes
library(org.Hs.eg.db)
myIds <- AnnotationDbi::select(org.Hs.eg.db, keys = "SLC25A3", keytype = "SYMBOL", 
                               columns = c("SYMBOL", "GENENAME","ENTREZID"))
entrezIDforSLC25A3 <- myIds[,"ENTREZID"]
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
allExons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
exonsforSLC25A3 <- allExons[[entrezIDforSLC25A3]]
seqlevels(exonsforSLC25A3) <- "chr12"
for(i in 1:length(exonsforSLC25A3)){
  myRegionOfInterest <- exonsforSLC25A3[i]
  myParam <- ScanBamParam(which = myRegionOfInterest)
  ReadsInExons <- readGAlignments("SortedHeart.bam",param = myParam)
  print(length(ReadsInExons))
}

library(GenomicAlignments)
sortedHeart <- sortBam("heart.bodyMap.bam","Heart")
indexBam(sortedHeart)
heartCoverage <- coverage("Heart.bam")
class(heartCoverage)
heartCoverage
chr12Cov <- heartCoverage[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()
heartAln <- readGAlignments("Heart.bam")
heartCov1 <- coverage(heartAln)
chr12Cov <- heartCov1[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()
heartGR <- granges(heartAln)
heartCov2 <- coverage(heartGR)
chr12Cov <- heartCov2[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal() # slightly different image
heartGRL <- grglist(heartAln)
heartCov3 <- coverage(heartGRL)
chr12Cov <- heartCov3[["chr12"]]
signalDepth <- chr12Cov[98591400:98608400]
signalDepthScaled <- data.frame(Position=98591400:98608400,
                                Signal=signalDepth*1000)
library(ggplot2)
ggplot(signalDepthScaled,aes(x=Position,y=Signal))+
  geom_line()+theme_minimal()
heartAlnPos <- heartAln[strand(heartAln) == "+"]
heartAlnPos <- coverage(heartAlnPos)
heartAlnPos["chr12"]
export.bw(heartAlnPos,con="heartPos.bw")
heartCoverageX10 <- coverage("Heart.bam",
                             weight = 10)
heartCoverageX10["chr12"]
heartCoverage["chr12"]
allChromosomeStats <- idxstatsBam("Heart.bam")
allChromosomeStats
mapped <- sum(allChromosomeStats[,"mapped"])
heartCoverageNorm <- coverage("Heart.bam",
                              weight = (10^6)/mapped)
heartCoverageNorm["chr12"] # export to .bw and view in IGV
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
exonsOfGenes <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,
                        by="gene")
slc25A3 <- exonsOfGenes[["5250"]]
slc25A3
heartCoverageNorm <- coverage("Heart.bam")
myMeanCovOverExons <- mean(heartCoverageNorm[slc25A3])
myMeanCovOverExons
geneBody <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSS <- promoters(geneBody,500,500)
myTSScounts <- summarizeOverlaps(TSS,"Heart.bam")
class(myTSScounts)
myTSScounts
countMatrix <- assay(myTSScounts)
countMatrix["5250",]
Granges <- rowRanges(myTSScounts)
Granges
geneExons <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
geneExons["5250"]
myGeneCounts <- summarizeOverlaps(geneExons,"Heart.bam")
myGeneCounts
countMatrix <- assay(myGeneCounts)
countMatrix["5250",]
grgList <- rowRanges(myGeneCounts)
grgList
allGeneCounts <- summarizeOverlaps(geneExons,
                                   c("Heart.bam","Liver.bam"))
countMatrix <- assay(allGeneCounts)
countMatrix["5250",]
myBam <- BamFile("Heart.bam")
class(myBam)
myBam <- BamFile("Heart.bam", yieldSize = 1000)
heartGeneCounts <- summarizeOverlaps(geneExons,myBam)
heartGeneCounts
myBam <- BamFileList(c("Heart.bam","Liver.bam"),
                     yieldSize = 1000)
allGeneCounts <- summarizeOverlaps(geneExons,myBam)
allGeneCounts
# Sort and index both the heart.bodyMap.bam and liver.bodyMap.bam files
library(GenomicAlignments)
library(Rsamtools)
sortBam("liver.bodyMap.bam","Liver")
indexBam("Liver.bam")
sortBam("heart.bodyMap.bam","Heart")
indexBam("Heart.bam")
# Calculate the coverage from the sorted, indexed BAM files normalised to the total mapped reads per sample
liverMapped <- idxstatsBam("Liver.bam")
totalLiver <- sum(liverMapped[,"mapped"])
liverCoverageNorm <- coverage("Liver.bam",weight = 1/totalLiver)
heartMapped <- idxstatsBam("Heart.bam")
totalHeart <- sum(heartMapped[,"mapped"])
heartCoverageNorm <- coverage("Heart.bam",weight = 1/totalHeart)
# Plot the coverage for Heart and Liver samples over region Chr12 98,986,183-98,998,558 using ggplot2
heartCoverageNorm12 <- heartCoverageNorm[["chr12"]]
signalDepthHeart <- heartCoverageNorm12[98592587:98607803]
signalDepthScaled <- data.frame(Sample="Heart",
                                Position=98592587:98607803,
                                Signal=signalDepthHeart)
liverCoverageNorm12 <- liverCoverageNorm[["chr12"]]
signalDepthLiver <- liverCoverageNorm12[98592587:98607803]
signalDepthScaled2 <- data.frame(Sample="Liver",
                                 Position=98592587:98607803,
                                 Signal=signalDepthLiver)
toPlot <- rbind(signalDepthScaled,signalDepthScaled2)
library(ggplot2)
ggplot(toPlot,aes(x=Position,y=Signal,colour=Sample))+
  geom_line()+theme_minimal()+ggtitle("Coverage over Chr12 98,986,183-98,998,558 for liver and heart samples")
# Using the normalised coverage for heart and liver samples,
# calculate the sum coverage within TSS (+500bp/-500bp) of the TMPO, SLC25A3 and IKBIP genes.
# Create a barplot of the sum coverage for each gene using ggplot2
library(org.Hs.eg.db)
myIds <- AnnotationDbi::select(org.Hs.eg.db, keys = c("SLC25A3","TMPO","IKBIP"), keytype = "SYMBOL", 
            columns = c("SYMBOL", "GENENAME","ENTREZID"))
entrezIDforGenes <- myIds[,"ENTREZID"]
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38Genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
myPromoters <- promoters(hg38Genes,500,500)
prms <- myPromoters[myPromoters$gene_id %in% entrezIDforGenes,]
liverSums <- sum(liverCoverageNorm[prms])
heartSums <- sum(heartCoverageNorm[prms])
lFfra <- data.frame(Sample="Liver",ENTREZID=prms$gene_id,Coverage=liverSums)
hFfra <- data.frame(Sample="Heart",ENTREZID=prms$gene_id,Coverage=heartSums)
covFrame <- rbind(lFfra,hFfra)
toPlot <- merge(myIds,covFrame,by="ENTREZID")
ggplot(toPlot,aes(x=SYMBOL,y=Coverage,fill=Sample))+geom_bar(stat="identity",position = "dodge")+theme_bw()
# Using the summarizeOverlaps() function, count reads from Heart and liver samples in the TMPO, SLC25A3 and IKBIP genes.
# Normalise reads to total number of mapped reads in the sample per million mapped reads and plot the mean counts per gene versus the log2 fold difference between heart and liver samples
library(org.Hs.eg.db)
myIds <- AnnotationDbi::select(org.Hs.eg.db, keys = c("SLC25A3","TMPO","IKBIP"), keytype = "SYMBOL", 
            columns = c("SYMBOL", "GENENAME","ENTREZID"))
entrezIDforGenes <- myIds[,"ENTREZID"]
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38Genes <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,by="gene")
exByGenes <- hg38Genes[names(hg38Genes) %in% entrezIDforGenes,]
sortBam("liver.bodyMap.bam","Liver")
indexBam("Liver.bam")
liverMapped <- idxstatsBam("Liver.bam")
totalLiver <- sum(liverMapped[,"mapped"])
sortBam("heart.bodyMap.bam","Heart")
indexBam("Heart.bam")
heartMapped <- idxstatsBam("Heart.bam")
totalHeart <- sum(heartMapped[,"mapped"])
newCounts <- summarizeOverlaps(exByGenes,c("Liver.bam","Heart.bam"))
cMatrix <- assay(newCounts)
cMatrix[,"Liver.bam"] <- (cMatrix[,"Liver.bam"]/totalLiver)*(10^6)
cMatrix[,"Heart.bam"] <- (cMatrix[,"Heart.bam"]/totalHeart)*(10^6)
toPlot <- data.frame(ENTREZID=rownames(cMatrix),mean=rowMeans(cMatrix),log2FC=log2(cMatrix[,"Liver.bam"]/cMatrix[,"Heart.bam"]))
toPlot <- merge(myIds,toPlot,by="ENTREZID")
ggplot(toPlot,aes(x=mean,y=log2FC))+geom_point()+theme_bw()+ylim(-3,3)+ggtitle("Liver_Minus_Heart")+geom_text(label=toPlot$SYMBOL)
