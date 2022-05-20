BiocManager::install(c("BSgenome",
                       "BSgenome.Hsapiens.UCSC.hg19.masked",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene",
                       "org.Hs.eg.db",
                       "ensembldb",
                       "EnsDb.Hsapiens.v75",
                       "AnnotationHub",
                       "rtracklayer",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "KEGGREST",
                       "KEGGgraph",
                       "RBGL",
                       "rols",
                       "GSEABase"))
# install.packages(c("R.utils", "png", "DT"))
library(BSgenome)
library(Biostrings)
ag = available.genomes()
length(ag)
grep("Scerev", ag, value=TRUE)
grep("Hsap", ag, value=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
head(genome(Hsapiens))
seqnames(Hsapiens)
length(Hsapiens)
class(Hsapiens)
methods(class="BSgenome")
Hsapiens$chrX
substr(Hsapiens$chrX, 5e6, 5.1e6)
nchar(Hsapiens$chrY)
nchar(Hsapiens[[24]])
sum(unlist(lapply(18:24,function(x) nchar(Hsapiens[[x]]))))

library(parallel)
system.time(sum(unlist(lapply(18:24,function(x) nchar(Hsapiens[[24]])))))
detectCores()
options(mc.cores=4)
system.time(sum(unlist(mclapply(18:24,function(x) nchar(Hsapiens[[24]])))))

# How many Bioconductor packages provide reference genomic sequence for zebrafish (Danio rerio). Exclude .masked suffix
grep("Drerio", ag, value = TRUE)

library(BSgenome.Hsapiens.UCSC.hg19.masked)
class(BSgenome.Hsapiens.UCSC.hg19.masked$chr17) # class
genome(BSgenome.Hsapiens.UCSC.hg19.masked)

BSgenome.Hsapiens.UCSC.hg19.masked$chr22 # In build hg19, what percentage of the length of chromosome 22 is occupied by "assembly gaps"?

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene # abbreviate
txdb
class(txdb)
methods(class="TxDb")
ghs = genes(txdb)
ghs
table(strand(ghs))
summary(width(ghs))
exons(txdb, columns=c("EXONID", "TXNAME", "GENEID"),
      filter=list(gene_id=c(100, 101)))
ghs[which.max(width(ghs))] # inspect largest gene in genome

# compare total size of exons to total size of genes
ex = exons(txdb)
rex = reduce(ex)
ex_width = sum(width(rex))    # bases in exons
gene_width = sum(width(genes(txdb)))    # bases in genes
ex_width / gene_width
transcripts(txdb) # ensembldb gives higher results

library(ensembldb)
library(EnsDb.Hsapiens.v75)
edb = EnsDb.Hsapiens.v75  # abbreviate
names(listTables(edb))
txs <- transcripts(edb, filter = GeneNameFilter("ZBTB16"),
                   columns = c("protein_id", "uniprot_id", "tx_biotype"))
txs
transcripts(edb)
table(transcripts(edb)$tx_biotype)



library(AnnotationHub)
# ah = AnnotationHub(localHub = FALSE)
ah
query(ah, "HepG2")
display(query(ah, "HepG2")) # can search all HepG2 data
query(ah, c("HepG2", "H4K5"))
query(ah, c("HepG2", "H3K4me3"))
query(query(ah, "HepG2"), "H3K4me3")
query(query(ah, "HepG2"), "H3K4me3")$tags
length(unique(ah$species))
subset(ah, species == "Homo sapiens")
query(query(query(ah, "HepG2"), "H3K4me3"), c("E118", "broadPeak"))$ah_id
AH29728 <- ah[[query(query(query(ah, "HepG2"), "H3K4me3"), c("E118", "broadPeak"))$ah_id]]

# cheatsheet http://genomicsclass.github.io/book/pages/bioc1_annoCheat.html

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db) # columns() gives same answer
head(select(org.Hs.eg.db, keys="ORMDL3", keytype="SYMBOL", 
            columns="PMID"))
select(org.Hs.eg.db, keys="286297", keytype = "ENTREZID", 
       columns = c("SYMBOL", "GENENAME"))

library(GO.db)
GO.db # metadata
select(GO.db, keys = keys(GO.db)[1:5], columns = columns(GO.db)[1:3])
keys(GO.db, keytype = "TERM")[1:5]

# find GOID (gene ontology tag) for ribosome biogenesis
select(GO.db, keys = "ribosome biogenesis", keytype="TERM", columns="GOID")

# find symbols for genes involved in ribosome biogenesis 
select(org.Hs.eg.db, keys = "GO:0042254", keytype="GO", columns="SYMBOL")

# you can pull out multiple columns at once
select(org.Hs.eg.db, keys = "GO:0042254", keytype="GO", columns=c("SYMBOL", "ENTREZID"))

# find gene ontology tags for related to ZNF658, which has the specified ENTREZID
select(org.Hs.eg.db, keys = "26149", keytype ="ENTREZID", columns ="GO")

# save GO tags to a character vector
select(org.Hs.eg.db, keys = "26149", keytype ="ENTREZID", columns ="GO")$"GO"
myk = .Last.value

# identify biological processes ZNF658 is involved in
select(GO.db, keys = myk, columns = "TERM")



library(DBI)
dbListTables(GO_dbconn())
dbGetQuery(GO_dbconn(), "select _id, go_id, term from go_term limit 5")
dbGetQuery(GO_dbconn(), "select * from go_bp_parents where _id=30")
dbGetQuery(GO_dbconn(), "select _id, go_id, term from go_term where _id=26616")
dbGetQuery(GO_dbconn(), "select * from go_bp_parents where _id=26616")
dbGetQuery(GO_dbconn(), "select _id, go_id, term from go_term where _id=5932")

library(KEGGREST)
brca2K = keggGet("hsa:675") # reference to a specific gene
names(brca2K[[1]])
names(keggGet("path:hsa05212")[[1]]) # info on a pathway
keggGet("path:hsa05212")[[1]]$GENE[seq(1,132,2)] # entrez gene ids for pathway

# inspect some entrez ids
select(org.Hs.eg.db, keys="5888", keytype="ENTREZID", columns ="SYMBOL")
select(org.Hs.eg.db, keys="675", keytype="ENTREZID", columns ="SYMBOL")

library(png)
library(grid)
brpng = keggGet("hsa05212", "image")
grid.raster(brpng)

library(KEGGgraph)
toyKGML <- system.file("extdata/kgml-ed-toy.xml", package="KEGGgraph")
toyGraph <- parseKGML2Graph(toyKGML, genesOnly=FALSE)
toyGraph
nodes(toyGraph)

library(Rgraphviz)
nodeInfo <- getKEGGnodeData(toyGraph)
nodeType <- sapply(nodeInfo, getType)
makeNodeRenderAttrs <- function(g, label=nodes(g),
                                shape="ellipse", fill="#e0e0e0",...) {
  rv <- list(label=label, shape=shape, fill=fill, ...)
  nA <- nodeRenderInfo(g)
  for(i in seq(along=rv)) {
    if (length(rv[[i]]) == 1) {
      rv[[i]] <- rep(rv[[i]], numNodes(g))
    } else {
      if (length(rv[[i]]) != numNodes(g))
        stop("Attribute vector must have as many elements as 'g' has nodes.")
    }
    names(rv[[i]]) <- nodes(g)
    nA[[ names(rv)[[i]] ]] <- rv[[i]]
  }
  nodeRenderInfo(g) <- nA
  return(g)
}
toyDrawn <- plotKEGGgraph(toyGraph)
toyDrawnRefine <- makeNodeRenderAttrs(toyDrawn, fill=ifelse(nodeType=="gene", "lightblue", "orange"),
                                      shape=ifelse(nodeType=="gene", "ellipse","rectangle"))
renderGraph(toyDrawnRefine)

makeAttr <- function(graph, default, valNodeList) {
  tmp <- nodes(graph)
  x <- rep(default, length(tmp)); names(x) <- tmp
  if(!missing(valNodeList)) {
    stopifnot(is.list(valNodeList))
    allnodes <- unlist(valNodeList)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList[[i]]] <- names(valNodeList)[i]
    }
  }
  return(x)
}

mapfile <- system.file("extdata/map00260.xml",package="KEGGgraph")
map <- parseKGML(mapfile)
map
reactions <- getReactions(map)
reactions[[1]]
chemicalGraph <- KEGGpathway2reactionGraph(map)
outDegrees <- sapply(edges(chemicalGraph), length)
maxout <- names(sort(outDegrees,decreasing=TRUE))[1:3]
nAttrs <- list()
maxoutlabel <- as.list(maxout); names(maxoutlabel) <- maxout
nAttrs$label <- makeAttr(chemicalGraph, "", maxoutlabel)
nAttrs$fillcolor <- makeAttr(chemicalGraph, "lightblue", list(orange=maxout))
nAttrs$width <- makeAttr(chemicalGraph,"0.8", list("1.8"=maxout))
plot(chemicalGraph, nodeAttrs=nAttrs)

if(require(SPIA)) {
  data(colorectalcancer,package="SPIA")
} else {
  load(system.file("extdata/colorectalcancerSPIA.RData", package="KEGGgraph"))
}
head(top)

library(hgu133plus2.db)
x <- hgu133plus2ENTREZID
top$ENTREZ <- unlist(as.list(x[top$ID]))
top <- top[!is.na(top$ENTREZ),]
top <- top[!duplicated(top$ENTREZ),]
tg1 <- top[top$adj.P.Val < 0.05,]
DE_Colorectal <- tg1$logFC
names(DE_Colorectal) <- as.vector(tg1$ENTREZ)
ALL_Colorectal <- top$ENTREZ
tmp <- "hsa05210.xml"
retrieveKGML(pathwayid="05210", organism="hsa", destfile=tmp)
colFile <- system.file("extdata/hsa05210.xml",
                       package="KEGGgraph")
g <- parseKGML2Graph(colFile)
deKID <- translateGeneID2KEGGID(names(DE_Colorectal))
allKID <- translateGeneID2KEGGID(ALL_Colorectal)
isDiffExp <- nodes(g) %in% deKID
sprintf("%2.2f%% genes differentially-expressed", mean(isDiffExp)*100)

library(RColorBrewer)
library(org.Hs.eg.db)
library(RBGL)
library(grid)
ar <- 20
cols <- rev(colorRampPalette(brewer.pal(6, "RdBu"))(ar))
logfcs <- DE_Colorectal[match(nodes(g), deKID)]
names(logfcs) <- nodes(g)
logfcs[is.na(logfcs)] <- 0
incol <- round((logfcs+2)*5); incol[incol>ar] <- ar
undetected <- !nodes(g) %in% allKID
logcol <- cols[incol]; logcol[logfcs==0] <- "darkgrey"; logcol[undetected] <- "yellow"
names(logcol) <- names(logfcs)
nA <- makeNodeAttrs(g, fillcolor=logcol, label="", width=10, height=1.2)
par(mar=c(3,5,0,5), mgp=c(0,0,0))
layout(mat=matrix(c(rep(1,8),2), ncol=1, byrow=TRUE))
plot(g, "dot", nodeAttrs=nA)
image(as.matrix(seq(1,ar)), col=cols, yaxt="n", xaxt="n")
mtext("down-regulation", side=1,  at=0, line=1)
mtext("up-regulation", side=1,  at=1, line=1)

gDeg <- degree(g)
gIsSingle <- gDeg[[1]] + gDeg[[2]] == 0
options(digits=3)
gGeneID <- translateKEGGID2GeneID(nodes(g))
gSymbol <-  sapply(gGeneID, function(x) mget(x, org.Hs.egSYMBOL, ifnotfound=NA)[[1]])
isUp <- logfcs > 0
isDown <- logfcs < 0
singleUp <- isUp & gIsSingle
singleDown <- isDown & gIsSingle

ups <- nodes(g)[logfcs > 0]
upNs <- unique(unlist(neighborhood(g, ups, return.self=TRUE)))
upSub <- subKEGGgraph(upNs, g)
upNeighbor <- nodes(upSub)[sapply(neighborhood(upSub, nodes(upSub)), length)>0]
upNeighbor <- setdiff(upNeighbor, nodes(g)[undetected])
upSub <- subKEGGgraph(upNeighbor, upSub)
upSubGID <- translateKEGGID2GeneID(nodes(upSub))
upSymbol <- gSymbol[upSubGID]
upnA <- makeNodeAttrs(upSub, fillcolor=logcol[nodes(upSub)], label=upSymbol, fixedsize=TRUE, width=10, height=10, font=20)
plot(upSub, "dot", nodeAttrs=upnA)



library(rols)
# oo = Ontologies()
oo
oo[[1]]

glis = OlsSearch("glioblastoma")
glis
res = olsSearch(glis)
dim(res)
resdf = as(res, "data.frame") # get content
resdf[1:4,1:4]
resdf[1,5] # label 
resdf[1,6] # full description for one instance

library(GSEABase)
glioG = getGmt(system.file("gmt/glioSets.gmt", package="ph525x"))
glioG
head(geneIds(glioG[[1]]))

library(Homo.sapiens)
class(Homo.sapiens)
Homo.sapiens
tx = transcripts(Homo.sapiens)
keytypes(Homo.sapiens)
columns(Homo.sapiens)

library(hgu133plus2.db)
hgu133plus2.db

library(hgu133plus2probe)
head(hgu133plus2probe)
dim(hgu133plus2probe)

select(hgu133plus2.db, keytype="PROBEID", 
       columns=c("SYMBOL", "GENENAME", "PATH", "MAP"), keys="1007_s_at")

library(rtracklayer)
?liftOver

# chromosome 1 gene locations in hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(tx38, pruning.mode="coarse") = "chr1"
g1_38 <- genes(tx38)

# download the hg38 to hg19 chain file
library(AnnotationHub)
ah <- AnnotationHub()
ah.chain <- subset(ah, rdataclass == "ChainFile" & species == "Homo sapiens")
AH14108 <- ah[[query(ah.chain, c("hg19", "hg38"))$ah_id[1]]]

# perform the liftOver
g1_19L <- liftOver(g1_38, AH14108)
g1_19L

f1 = dir(system.file("extdata",package="ERBS"), full=TRUE)[1] # access data
readLines(f1, 4) # look at a few lines

library(rtracklayer)
imp = import(f1, format="bedGraph")
imp
genome(imp)
genome(imp) <- "hg19" # set genome
imp
genome(imp)
export(imp, "demoex.bed")  # implicit format choice, export as BED format
cat(readLines("demoex.bed", n=5), sep="\n") # check output file

