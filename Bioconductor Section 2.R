BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))

library(Biobase)
library(GenomicRanges)
library(ERBS)
data(HepG2)  # cell line of liver origin
data(GM12878)  # immortalized B cell
HepG2  # locations of ER binding peaks
class(HepG2) # What is the class of HepG2 ?
values(HepG2) # or, mcols(HepG2)
length(HepG2) # How many regions are represented?
median(HepG2$signalValue) # What is the median of the signalValue column for the HepG2 data?
as.character(HepG2[which.max((HepG2$signalValue))]@seqnames@values) # In what chromosome is the region with the highest signalValue
length(which(HepG2@seqnames@values=="chr16")) # How many regions are from chromosome 16
median(width(HepG2)) # What is the median width?
hist(width(HepG2)) # heavy right tail

chr <- seqnames(HepG2)
chr
as.character(chr)
table(chr)
table(chr)[1:24]
HepG2[chr=="chr20",] # subset of chr20
granges(HepG2[chr=="chr20",])
x = HepG2[order(HepG2),]
seqnames(x)     # demonstrate usefulness of Rle type
as.character(seqnames(x))

library(IRanges)

ir <- IRanges(5,10)
ir
ir <- IRanges(5, width = 6)
ir
start(ir)
end(ir)
width(ir)
IRanges(start=c(3,5,17), end=c(10,8,20))
shift(ir, -2)
narrow(ir, start=2)
narrow(ir, end=5)
flank(ir, width=3, start=TRUE, both=FALSE)
flank(ir, width=3, start=FALSE, both=FALSE)
flank(ir, width=3, start=TRUE, both=TRUE)
flank(ir, width=3, start=FALSE, both=TRUE)
ir * 2
ir * -2
ir + 2
ir - 2
resize(ir, 1)

plot(0,0,xlim=c(0,23),ylim=c(0,14),type="n",xlab="",ylab="",xaxt="n")
axis(1,0:15)
abline(v=0:14 + .5,col=rgb(0,0,0,.5))
# plot the original IRange
plotir <- function(ir,i) { arrows(start(ir)-.5,i,end(ir)+.5,i,code=3,angle=90,lwd=3) }
plotir(ir,1)
# draw a red shadow for the original IRange
polygon(c(start(ir)-.5,start(ir)-.5,end(ir)+.5,end(ir)+.5),c(-1,15,15,-1),col=rgb(1,0,0,.2),border=NA)
# draw the different ranges
plotir(shift(ir,-2), 2)
plotir(narrow(ir, start=2), 3)
plotir(narrow(ir, end=5), 4)
plotir(flank(ir, width=3, start=TRUE, both=FALSE), 5)
plotir(flank(ir, width=3, start=FALSE, both=FALSE), 6)
plotir(flank(ir, width=3, start=TRUE, both=TRUE), 7)
plotir(flank(ir, width=3, start=FALSE, both=TRUE), 8)
plotir(ir * 2, 9)
plotir(ir * -2, 10)
plotir(ir + 2, 11)
plotir(ir - 2, 12)
plotir(resize(ir, 1), 13)
text(rep(15,12), 1:12, c("ir","shift(ir,-2)","narrow(ir,start=2)",
                         "narrow(ir,end=5)",
                         "flank(ir, start=T, both=F)",
                         "flank(ir, start=F, both=F)",
                         "flank(ir, start=T, both=T)",
                         "flank(ir, start=F, both=T)",
                         "ir * 2","ir * -2","ir + 2","ir - 2",
                         "resize(ir, 1)"), pos=4)

ir <- IRanges(start=c(3,5,17), end=c(10,8,20))
range(ir)
reduce(ir)
gaps(ir)
disjoin(ir)

ir <- IRanges(start=101,end=200)
ir
ir*2
start(ir*2) # What is the starting point of the resulting range
start(narrow(ir,start = 20)) # what is the new starting point of the range
width(ir+25) # what is the width of the resulting range
sum(width(IRanges(start = c(1,11,21),end=c(3,15,27)))) # What is the sum of the widths of all the ranges

ir1 <- IRanges(start = c(101,106,201,211,221,301,306,311,351,361,401,411,501),end=c(150,160,210,270,225,310,310,330,390,380,415,470,510))

library("ph525x")
plotRanges(ir1)
sum(width(gaps(ir1))) # What is the total width from 101 to 510 which is not covered by ranges in x?
length(disjoin(ir1)) # How many disjoint ranges are contained within the ranges in ir1 from the previous question?

library(rafalib)
mypar(2,1)
plotRanges(ir1)
plotRanges(resize(ir1,1))
all.equal(start(resize(ir1,1)), end(resize(ir1,1)))
# it gives you just the starting point of each range

library(GenomicRanges)
gr <- GRanges("chrZ", IRanges(start=c(5,10),end=c(35,45)),
              strand="+", seqlengths=c(chrZ=100L))
gr
genome(gr) <- "hg19"
gr
seqnames(gr)
seqlengths(gr)
shift(gr, 10)
shift(gr, 80)
trim(shift(gr, 80))

mcols(gr)
mcols(gr)$value <- c(-1,4)
gr
mcols(gr)$value <- NULL
gr

gr2 <- GRanges("chrZ",IRanges(11:13,51:53))
grl <- GRangesList(gr, gr2)
grl
length(grl)
elementNROWS(grl)
grl[[1]]
width(grl)
sum(width(grl))
mcols(grl)$value <- c(5,7)
grl
mcols(grl)

(gr1 <- GRanges("chrZ",IRanges(c(1,11,21,31,41),width=5),strand="*"))
(gr2 <- GRanges("chrZ",IRanges(c(19,33),c(38,35)),strand="*"))
fo <- findOverlaps(gr1, gr2)
fo
queryHits(fo)
subjectHits(fo)
gr1 %over% gr2
gr1[gr1 %over% gr2]

gr1 <- GRanges("chrZ",IRanges(1,10),strand="+")
gr2 <- GRanges("chrZ",IRanges(1,10),strand="-")
gr1 %over% gr2

(r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20))))
str(r)
as.numeric(r)
(v <- Views(r, start=c(4,2), end=c(7,6)))
str(v)

ir <- IRanges(c(3, 8, 14, 15, 19, 34, 40),
              width = c(12, 6, 6, 15, 6, 2, 7))

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}

par(mfrow=c(4,1), mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))

gir = GRanges(seqnames="chr1", ir, strand=c(rep("+", 4), rep("-",3)))

plotGRanges = function (x, xlim = x, col = "black", sep = 0.5, xlimits = c(100,1100), ...) {
  main = deparse(substitute(x))
  ch = as.character(seqnames(x)[1])
  x = ranges(x)
  height <- 1
  if (is(xlim, "Ranges")) 
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim = xlimits, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom, end(x) + 0.5, ybottom + height, 
       col = col, ...)
  title(main, xlab = ch)
  axis(1)
}

par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlimits = c(0,50))
plotGRanges(resize(gir,1), xlimits = c(0,50), col="green")
plotGRanges(flank(gir,3), xlimits = c(0,50), col="purple")
plotGRanges(flank(gir,2,start=FALSE), xlimits = c(0,50), col="brown")

x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
# x = GRanges("chr1", IRanges(c(101,201,401,501), width = c(150,250,450,550)), strand = "+")
# y = GRanges("chr1", IRanges(c(101,221,301,401,541), width = c(150,250,350,470,550)), strand = "+")
par(mfrow=c(2,1))
plotGRanges(x, xlimits = c(0,1200), col="green")
plotGRanges(y, xlimits = c(0,1200), col="red")

GRangesList(x,y) # keep the information about which set the ranges belong to
subsetByOverlaps(x, y)
c(x,y) # combine them into a single GRanges
par(mfrow=c(2,1))
plotGRanges(subsetByOverlaps(x, y), xlimits = c(0,1200), col="purple")
plotGRanges(c(x,y), xlimits = c(0,1200), col="brown")

abs(start(ranges(disjoin(c(x,y))))-start(ranges(c(x,y))))
sum(width(x))+sum(width(y))-sum(width(ranges(disjoin(c(x,y)))))
sum(width(disjoin(c(x,y)))) # total width which is in x or y but not in both?
sum(width(x))+sum(width(y))-sum(width(ranges(intersect(x,y))))
ranges(x)

dis = disjoin(c(x,y))
both = dis %over% x & dis %over% y
sum(width(dis[both]))

not =! (dis %over% x & dis %over% y)
sum(width(dis[not]))

z = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="-")
sum(x %over% z)



library(ArrayExpress)
if (!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")
library(minfi)
pref = unique(substr(dir(patt="idat"),1,17)) # find the prefix strings
raw = read.metharray(pref)
glioMeth = preprocessQuantile(raw) # generate SummarizedExperiment

MbyGene = function(mset, symbol="TP53", rad=5000) {
  # phase 1: annotated GRanges for the gene
  require(erma)
  require(Gviz)
  gmod = suppressMessages(genemodel(symbol))     # erma utility
  gseq = as.character(seqnames(gmod)[1])
  gmod$transcript = symbol
  # phase 2: filter down to the region of interest
  mlim = mset[which(seqnames(mset)==gseq),] # restrict to chromosome
  # now focus the methylation data to vicinity of gene
  d1 = subsetByOverlaps(GRanges(rowRanges(mlim),,, getM(mlim)),
                        range(gmod)+rad)
  # phase 3: use Gviz
  plotTracks(list(DataTrack(d1), 
                  GeneRegionTrack(gmod, 
                                  transcriptAnnotation="transcript", name=gseq), 
                  GenomeAxisTrack(name=gseq, showTitle=TRUE)))
}

MbyGene(glioMeth, symbol="TERT")




library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(ERBS)
data(HepG2)
data(GM12878)
res = findOverlaps(HepG2, GM12878)
class(res)
res
index = queryHits(res)
erbs = HepG2[index,]
erbs
granges(erbs)
erbs

HepG2[17,]
x <- distanceToNearest(GM12878,HepG2[17,])
x
m <- which(mcols(x)$distance==min(mcols(x)$distance))
m
GM12878[queryHits(x[m,]),]
min(mcols(x)$distance) # What is the distance between the 17th region of HepG2 and its closest region in GM12878?
length(HepG2)
d = distanceToNearest(HepG2,GM12878)
mean(mcols(d)$distance < 2000) # What proportion of these distances are smaller than 2000 base pairs?
order(erbs)
erbs = granges(erbs)
erbs[order(erbs),]
erbs2 = intersect(HepG2,GM12878)
erbs2[order(erbs2),]



library(Homo.sapiens)
ghs = genes(Homo.sapiens)
ghs
length(ghs)
dim(ghs)
which.max(table(seqnames(ghs)))
hist(width(granges(ghs)))
median(width(granges(ghs)))

tssgr = resize(ghs,1)
granges(ghs)[1:3,]
granges(tssgr)[1:3,]

(index <- precede(erbs,ghs))
ghs[index[1:3]]
erbs[1:3]
# HepG2[1:3,]
# granges(HepG2)[1:3,]
# granges(ghs)[index[1:3],]
(d = distance(erbs,ghs[index]))
tssgr
(d = distanceToNearest(erbs, tssgr))
try(d[3])
queryLength(d)
dists = values(d)$distance
hist(dists, nc = 1000, xlim = c(0,100000))
(index = subjectHits(d))
(index = subjectHits(d)[dists < 1000])
sdists = ifelse(end(HepG2) < start(tssgr[index]), dists, -dists)
hist(sdists, xlim=c(-100000,100000), main="Density of d(ER binding peak, nearest TSS)" ,breaks=seq(min(sdists),max(sdists),len=1000))
abline(v=-c(10000,10000), lty=2)

columns(Homo.sapiens)
keytypes(Homo.sapiens)
geneids <- as.character(values(tssgr[index])$GENEID)
result = AnnotationDbi::select(Homo.sapiens,keys=geneids,columns=c("SYMBOL","GENENAME"),keytype="GENEID")
result[1:2,]

tss <- resize(ghs,1)
start(tss["100113402"]) # What is the TSS (Transcription Start Site) of the gene with ID: 100113402?
d <- subjectHits(distanceToNearest(erbs[4,], tss)) # or
d <- nearest(erbs[4], tss)
geneids <- as.character(values(tss[d,])$GENEID) # or mcols
AnnotationDbi::select(Homo.sapiens,keys=geneids,columns=c("SYMBOL","ENTREZID"),keytype="GENEID")
AnnotationDbi::select(Homo.sapiens,keys=geneids,columns=c("SYMBOL","GENENAME"),keytype="GENEID")



erbs2 = GenomicRanges::intersect(HepG2,GM12878)
(index2 <- precede(erbs2,ghs))
(d2 = distance(erbs,ghs[index2]))
tssgr
(d2 = distanceToNearest(erbs, tssgr))
queryLength(d2)
dists2 = values(d2)$distance
hist(dists2, nc = 1000, xlim = c(0,100000))
(index2 = subjectHits(d2))
(index2 = subjectHits(d2)[dists2 < 1000])
geneids2 <- as.character(values(tssgr[index2])$GENEID)
result2 = AnnotationDbi::select(Homo.sapiens,keys=geneids2,columns="GENENAME",keytype="GENEID")
result2[1:2,]



library(Biostrings)
(dna <- DNAString("TCGAGCAAT"))
length(dna)
DNAString("JQX")
DNAString("NNNACGCGC-TTA-CGGGCTANN")
dna[4:6]
as.character(dna)
(set1 <- DNAStringSet(c("TCA", "AAATCG", "ACGTGCCTA", "CGCGCA", "GTT", "TCA")))
set1[2:3]
set1[[4]] # preferable if extract one string
set1[4]
length(set1)
width(set1)
duplicated(set1)    # detect which sequences are duplicated
unique(set1)    # keep only unique sequences
sort(set1)

dna_seq <- DNAString("ATCGCGCGCGGCTCTTTTAAAAAAACGCTACTACCATGTGTGTCTATC")
letterFrequency(dna_seq, "A")
letterFrequency(dna_seq, "GC") # count G or C in sequence
dinucleotideFrequency(dna_seq)
trinucleotideFrequency(dna_seq)
reverseComplement(dna_seq)
translate(dna_seq)
countPattern("CG", dna_seq)    # pattern "CG" occurs 5 times
matchPattern("CG", dna_seq)
start(matchPattern("CG", dna_seq))
matchPattern("CTCTTTTAAAAAAACGCTACTACCATGTGT", dna_seq)
countPattern("TAG", dna_seq)
countPattern(reverseComplement(DNAString("TAG")), dna_seq)

set2 <- DNAStringSet(c("AACCGGTTTCGA", "CATGCTGCTACA", "CGATCGCGCCGG", "TACAACCGTACA"))
vcountPattern("CG", set2)
vmatchPattern("CG", set2)
vmatchPattern("CG", set2)[[1]]

eco <- DNAString("GGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGA")
eco
length(eco)
translate(eco)
countPattern("ATG",eco) # How many potential start codons are in the eco sequence?
matchPattern("ATG",eco) # Find the locations of these ATG trinucleotides.
start(matchPattern("ATG",eco)[1]) # What is the start location of the first ATG?
AA <- translate(eco[start(matchPattern("ATG",eco))[1]:length(eco)])
length(AA)
matchPattern("*",AA)
AA[1:start(matchPattern("*",AA))-1]
length(AA[1:start(matchPattern("*",AA))-1])
n <- as.numeric(letterFrequency(AA,"DE")) # number of negatively charged AA
p <- as.numeric(letterFrequency(AA,"HKR")) # number of positively charged AA
p-n # net charge of AAstring



library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
c17 = Hsapiens$chr17
c17

?getSeq
class(Hsapiens)
showMethods("getSeq")

# collection of DNA strings with ChIP-seq binding peaks
data(HepG2)
hepseq = getSeq(Hsapiens, HepG2)
length(HepG2)    # same number of sequences
width(HepG2)[1:5]    # widths match

# collection of shifted DNA strings with no relationship to binding sequences - essentially random
rhepseq = getSeq(Hsapiens, shift(HepG2, 2500))

# count occurrences of a motif in DNA sequences
mot = "TCAAGGTCA"
?vmatchPattern
vcountPattern(mot, hepseq)

# consider both forward matches and reverse complement matches 
sum(vcountPattern(mot, hepseq))    # forward pattern match
sum(vcountPattern(mot, reverseComplement(hepseq)))    # reverse pattern match

## compare motif occurrence in binding peak to random upstream sequences
# count of motifs in binding peaks
sum(vcountPattern(mot, hepseq)) +
  sum(vcountPattern(mot, reverseComplement(hepseq)))
# count of motifs in randomly selected regions of equal length
sum(vcountPattern(mot, rhepseq)) +
  sum(vcountPattern(mot, reverseComplement(rhepseq)))

# for real analysis, use MotifDb package, probabilistic binding packages like MEME and FIMO

library(stringr)
erbseq = getSeq(Hsapiens,erbs)
seqs= 1:75
props <- function (x) {
  (str_count(as.character(erbseq[[x]]),"C") +
     str_count(as.character(erbseq[[x]]),"G"))/length(erbseq[[x]])
  }
meds <- sapply(seqs, props)
median(meds)

control = shift(erbs,10000)
controlseqs = getSeq(Hsapiens,control)
gc = alphabetFrequency(controlseqs)[,2:3]
n = width(control)
controlgccontent = rowSums(gc)/n
median(controlgccontent)
