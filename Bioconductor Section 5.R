BiocManager::install(c("SpikeInSubset",
                       "genefilter",
                       "qvalue",
                       "limma",
                       "genomicsclass/maPooling",
                       "leukemiasEset",
                       "preprocessCore"))
library(Biobase)
library(maPooling)
data(maPooling)
head(pData(maPooling))

library(rafalib)
mypar()
flipt <- function(m) t(m[nrow(m):1,])
myimage <- function(m,...) {
  image(flipt(m),xaxt="n",yaxt="n",...)
  }

myimage(as.matrix(pData(maPooling)),col=c("white","black"),
        xlab="experiments",
        ylab="individuals",
        main="phenoData")

data(maPooling)
pd=pData(maPooling)
pooled=which(rowSums(pd)==12)
factor(as.numeric(grepl("b",names(pooled))))

i=11425;j=11878
pooled_y=exprs(maPooling[,pooled])
pooled_g=factor(as.numeric(grepl("b",names(pooled))))
mypar(1,2)
stripchart(split(pooled_y[i,],pooled_g),vertical=TRUE,method="jitter",col=c(1,2),
           main="Gene 1",xlab="Group",pch=15)
stripchart(split(pooled_y[j,],pooled_g),vertical=TRUE,method="jitter",col=c(1,2),
           main="Gene 2",xlab="Group",pch=15)

library(genefilter)
pooled_tt=rowttests(pooled_y,pooled_g)
pooled_tt$p.value[i]
pooled_tt$p.value[j]

individuals=which(rowSums(pd)==1)
individuals=individuals[-grep("tr",names(individuals))]
y=exprs(maPooling)[,individuals]
g=factor(as.numeric(grepl("b",names(individuals))))

technicalsd <- rowSds(pooled_y[,pooled_g==0])
biologicalsd <- rowSds(y[,g==0])
LIM=range(c(technicalsd,biologicalsd))
mypar(1,1)
boxplot(technicalsd,biologicalsd,names=c("technical","biological"),ylab="standard deviation")

mypar(1,2)
stripchart(split(y[i,],g),vertical=TRUE,method="jitter",col=c(1,2),xlab="Gene 1",pch=15)
points(c(1,2),tapply(y[i,],g,mean),pch=4,cex=1.5)
stripchart(split(y[j,],g),vertical=TRUE,method="jitter",col=c(1,2),xlab="Gene 2",pch=15)
points(c(1,2),tapply(y[j,],g,mean),pch=4,cex=1.5)

library(genefilter)
tt=rowttests(y,g)
tt$p.value[i]
tt$p.value[j]

NsigAt01 = sum(tt$p.value<0.01)
NsigAt01
NsigAt05 = sum(tt$p.value<0.05)
NsigAt05

set.seed(0)
shuffledIndex <- factor(sample(c(0,1),sum(g==0),replace=TRUE ))
nulltt <- rowttests(y[,g==0],shuffledIndex)
NfalselySigAt01 = sum(nulltt$p.value<0.01)
NfalselySigAt01
NfalselySigAt05 = sum(nulltt$p.value<0.05)
NfalselySigAt05

library(qvalue)
qvals = qvalue(tt$p.value)$qvalue
sum(qvals<0.05)
sum(qvals<0.01)

nullqvals = qvalue(nulltt$p.value)$qvalue
sum(nullqvals<0.05)
sum(nullqvals<0.01)

pool = exprs(maPooling)[,pooled] 
indiv = exprs(maPooling)[,individuals]
strain = ifelse(grepl("a",rownames(pData(maPooling))),0,1)
g_pool = strain[pooled]
g_indiv = strain[individuals]

# multiple testing

set.seed(1)
population = unlist( read.csv("femaleControlsPopulation.csv") )
alpha <- 0.05
N <- 12
m <- 10000
pvals <- replicate(m,{
  control = sample(population,N)
  treatment = sample(population,N)
  t.test(treatment,control)$p.value
})
pvals
sum(pvals < 0.05)

alpha <- 0.05
N <- 12
p0 <- 0.90 # 10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

set.seed(1)
calls <- sapply(1:m, function(i){
  control <- sample(population,N)
  treatment <- sample(population,N)
  if(!nullHypothesis[i]) treatment <- treatment + delta
  ifelse( t.test(treatment,control)$p.value < alpha, 
          "Called Significant", # if condition
          "Not Called Significant") # else condition
})

null_hypothesis <- factor( nullHypothesis, levels=c("TRUE","FALSE"))
table(null_hypothesis,calls)

B <- 10 # number of simulations
VandS <- replicate(B,{
  calls <- sapply(1:m, function(i){
    control <- sample(population,N)
    treatment <- sample(population,N)
    if(!nullHypothesis[i]) treatment <- treatment + delta
    t.test(treatment,control)$p.val < alpha
  })
  cat("V =",sum(nullHypothesis & calls), "S =",sum(!nullHypothesis & calls),"\n")
  c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
  })

B <-10000
minpval <- replicate(B, min(runif(10000,0,1))<0.01)
mean(minpval>=1)

set.seed(1)
pvals <- sapply(1:m, function(i){
  control <- sample(population,N)
  treatment <- sample(population,N)
  if(!nullHypothesis[i]) treatment <- treatment + delta
  t.test(treatment,control)$p.value
})
sum(pvals < 0.05/10000)

set.seed(1)
pvals.6 <- sapply(1:m, function(i){
  control <- sample(population,6)
  treatment <- sample(population,6)
  if(!nullHypothesis[i]) treatment <- treatment + delta
  t.test(treatment,control)$p.value
})
sum(pvals.6 < 0.05/10000)

library(genefilter) # rowttests is here
set.seed(1)
g <- factor(c(rep(0,N),rep(1,N))) # Define groups to be used with rowttests
B <- 1000 # number of simulations
Qs <- replicate(B,{
  controls <- matrix(sample(population, N*m, replace=TRUE),
                     nrow=m) # matrix with control data (rows are tests, columns are mice)
  treatments <-  matrix(sample(population, N*m, replace=TRUE),
                        nrow=m)  # matrix with control data (rows are tests, columns are mice)
  treatments[which(!nullHypothesis),] <-
    treatments[which(!nullHypothesis),] + delta # add effect to 10% of them
  dat <- cbind(controls,treatments) # combine to form one matrix
 calls <- rowttests(dat,g)$p.value < alpha
 R = sum(calls)
 Q = ifelse(R>0,sum(nullHypothesis & calls)/R,0)
 return(Q)
})
Qs

library(rafalib)
mypar(1,1)
hist(Qs) # Q is a random variable, this is its distribution
FDR = mean(Qs)
FDR

set.seed(1)
controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
dat <- cbind(controls,treatments)
pvals <- rowttests(dat,g)$p.value 

h <- hist(pvals,breaks=seq(0,1,0.05))
polygon(c(0,0.05,0.05,0),c(0,0,h$counts[1],h$counts[1]),col="grey")
abline(h=m0/20)

h <- hist(pvals,breaks=seq(0,1,0.01))
polygon(c(0,0.01,0.01,0),c(0,0,h$counts[1],h$counts[1]),col="grey")
abline(h=m0/100)

alpha <- 0.05
i = seq(along=pvals)
mypar(1,2)
plot(i,sort(pvals))
abline(0,i*alpha/m)
plot(i[1:30],sort(pvals)[1:30],main="Close-up") # close-up
abline(0,i*alpha/m)

# Identify outliers
# outliers <- boxplot(df$mpg, plot = FALSE)$out
# Remove outliers
# df[!(df$mpg %in% outliers),]

k <- max(which(sort(pvals) < i*alpha/m))
cutoff <- sort(pvals)[k]
cat("k =",k,"p-value cutoff=",cutoff)

fdr <- p.adjust(pvals, method="fdr")
mypar(1,1)
plot(pvals,fdr,log="xy")
abline(h=alpha,v=cutoff) # cutoff was computed above

alpha <- 0.05
B <- 1000 # number of simulations. We should increase for more precision
res <- replicate(B,{
  controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
  treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
  treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
  dat <- cbind(controls,treatments)
  pvals <- rowttests(dat,g)$p.value 
  # then the FDR
  calls <- p.adjust(pvals,method="fdr") < alpha
  R=sum(calls)
  Q=ifelse(R>0,sum(nullHypothesis & calls)/R,0)
  return(c(R,Q))
})
Qs <- res[2,]
mypar(1,1)
hist(Qs) # Q is a random variable, this is its distribution
FDR = mean(Qs)
FDR

Rs <- res[1,]
mean(Rs==0)*100 # 1.3 % of the simulations, we did not call any genes significant

hist(pvals,breaks=seq(0,1,0.05),freq=FALSE)
lambda = 0.1
pi0 = sum(pvals > lambda)/((1-lambda)*m)
abline(h = pi0)
print(pi0) # this is close to the true pi0 = 0.9

library(qvalue)
res <- qvalue(pvals)
qvals <- res$qvalues
plot(pvals,qvals)
res$pi0



library(rafalib)
library(SpikeInSubset)
data(rma95)
fac <- factor(rep(1:2,each=3))
pData(rma95)
par(mfrow=c(2,2))
for (i in 1:4) {
  spg = names(pData(rma95))
  plot(1:6, exprs(rma95)[spg[i+6],], main=spg[i+6], ylab="RMA",
    xlab="nominal", axes=FALSE)
  axis(2)
  axis(1, at=1:6, labels=pData(rma95)[[spg[i+6]]])
}

library(genefilter)
rtt <- rowttests(exprs(rma95),fac)
mask <- with(rtt, abs(dm) < .2 & p.value < .01)
spike <- rownames(rma95) %in% colnames(pData(rma95))
cols <- ifelse(mask,"red",ifelse(spike,"dodgerblue","black"))
mypar()
with(rtt, plot(-dm, -log10(p.value), cex=.8, pch=16,
     xlim=c(-1,1), ylim=c(0,5),
     xlab="difference in means",
     col=cols))
abline(h=2,v=c(-.2,.2), lty=2)
sum(p.adjust(rtt$p.value,method = "BH")[spike] < 0.05)
table(top50 = rank(rtt$p.value) <= 10, spike) # t-stat and p-val rank is the same

rtt$s <- apply(exprs(rma95), 1, function(row) sqrt(.5 * (var(row[1:3]) + var(row[4:6]))))
with(rtt, plot(s, -log10(p.value), cex=.8, pch=16,
              log="x",xlab="estimate of standard deviation",
              col=cols))

library(limma)
options(digits=3)
fit <- lmFit(rma95, design=model.matrix(~ fac))  # step 1 least squares estimates
colnames(coef(fit))
fit <- eBayes(fit) # step 2 moderate the t statistics
tt <- topTable(fit, coef=2)  # step 3 report
tt
dim(tt)

dim(topTable(fit, coef=2, number=Inf, sort.by="none"))
limmares <- data.frame(dm=coef(fit)[,"fac2"], p.value=fit$p.value[,"fac2"])
with(limmares, plot(dm, -log10(p.value),cex=.8, pch=16,
     col=cols,xlab="difference in means",
     xlim=c(-1,1), ylim=c(0,5)))
abline(h=2,v=c(-.2,.2), lty=2)
table(top50 = rank(limmares$p.value) <= 10, spike) 

n <- 40
qs <- seq(from=0,to=.2,length=n)
idx <- sapply(seq_len(n),function(i) which(as.integer(cut(rtt$s^2,qs)) == i)[1])
idx <- idx[!is.na(idx)]

par(mar=c(5,5,2,2))
plot(1,1,xlim=c(0,.21),ylim=c(0,1),type="n",
     xlab="variance estimates",ylab="",yaxt="n")
axis(2,at=c(.1,.9),c("before","after"),las=2)
segments((rtt$s^2)[idx],rep(.1,n),
         fit$s2.post[idx],rep(.9,n))

# Gene set analysis

library(GSE5859Subset)
data(GSE5859Subset)
library(sva)
library(limma)
library(matrixStats)
library(rafalib)
X = sampleInfo$group
mod <- model.matrix(~X)
svafit <- sva(geneExpression,mod)

svaX<-model.matrix(~X+svafit$sv)
lmfit <- lmFit(geneExpression,svaX)
tt<-lmfit$coef[,2]*sqrt(lmfit$df.residual)/(2*lmfit$sigma)
pval <- 2*(1-pt(abs(tt),lmfit$df.residual[1]))
qval <- p.adjust(pval,"BH")

gsets <- getGmt("c1.all.v7.4.entrez.gmt") # must provide correct path to file
length(gsets)
head(names(gsets))
gsets[["chrYq11"]]
head(geneIds(gsets[["chrYq11"]]))

mapGMT2Affy <- function(object,gsets){
  ann<-annotation(object)
  dbname<-paste(ann,"db",sep=".")
  require(dbname,character.only=TRUE)
  gns<-featureNames(object)
  # This call may generate warnings
  map<-select(get(dbname), keys=gns,columns=c("ENTREZID", "PROBEID"))
  map<-split(map[,1],map[,2])
  indexes<-sapply(gsets,function(ids){
    gns2<-unlist(map[geneIds(ids)])
    match(gns2,gns)
    })
  names(indexes)<-names(gsets)
  return(indexes)
}
# create an Expression Set
rownames(sampleInfo) <- colnames(geneExpression)
e = ExpressionSet(assay=geneExpression,
                  phenoData=AnnotatedDataFrame(sampleInfo),
                  annotation="hgfocus")
# can safely ignore the warning
gsids <- mapGMT2Affy(e,gsets)

# Approaches based on association tests
tab <- table(ingenset = 1:nrow(e) %in% gsids[["chrYq11"]], signif = qval < 0.05)
tab
chisq.test(tab)$p.val

library(rafalib)
mypar(1,1)
qs <- seq(0,1,len = length(tt)+1)-1/(2*length(tt));qs<-qs[-1]
qqplot(qt(qs,lmfit$df.resid),tt,ylim=c(-10,10),xlab="Theoretical quantiles",ylab="Observed") # note we leave some of the obvious ones out
abline(0,1)

# Gene set summary statistics

ind <- gsids[["chrXp11"]]
mypar(1,1)
plot(density(tt[-ind]),xlim=c(-7,7),main="",xlab="t-stat",sub="",lwd=4)
lines(density(tt[ind],bw=.7),col=2,lty=2,lwd=4)
rug(tt[ind],col=2)
legend(-6.5, .3, legend=c("present on xp11", "not on xp11"),
       col=c(2,1), lty=c(2,1), lwd=4)

es <- lmfit$coef[,2]
wilcox <- t(sapply(gsids,function(i){
  if(length(i)>2){
    tmp<- wilcox.test(es[i],es[-i])  
    n1<-length(i);n2<-length(es)-n1
    z <- (tmp$stat -n1*n2/2) / sqrt(n1*n2*(n1+n2+1)/12)
    return(c(z,tmp$p.value))
  } else return(rep(NA,2))  
}))
mypar(1,1)
cols <- rep(1,nrow(wilcox))
cols[grep("chrX",rownames(wilcox))]<-2
cols[grep("chrY",rownames(wilcox))]<-3
qqnorm(wilcox[,1],col=cols,cex=ifelse(cols==1,1,2),pch=16)
qqline(wilcox[,1])
legend("topleft",c("Autosome","chrX","chrY"),pch=16,col=1:3,box.lwd=0)

avgt <- sapply(gsids,function(i) sqrt(length(i))*mean(tt[i]))
qqnorm(avgt,col=cols,cex=ifelse(cols==1,1,2),pch=16)
qqline(avgt)
legend("topleft",c("Autosome","chrX","chrY"),pch=16,col=1:3,box.lwd=0)

avgt[order(-abs(avgt))[1:10]]

# Hypothesis testing for gene sets
N <- sapply(gsids,length)
ind1<-which(X==0)
ind2<-which(X==1)
corrs <- t(sapply(gsids,function(ind){
  if(length(ind)>=2){
    cc1<-cor(t(geneExpression[ind,ind1]))
    cc2<-cor(t(geneExpression[ind,ind2]))
  return(c(median(cc1[lower.tri(cc1)]),
    median(cc2[lower.tri(cc2)])))
  } else return(c(NA,NA))
}))
mypar(1,1)
plot(corrs[N>10,],xlim=c(-1,1),ylim=c(-1,1),xlab="Correlation within females",ylab="Correlation within males")
abline(h=0,v=0,lty=2)

avgcorrs <- rowMeans(corrs)
cf <- (1+(N-1)*avgcorrs)
cf[cf<1] <- 1 # we ignore negative correlations
correctedavgt <- avgt/sqrt(cf)
parampvaliid <- 2*pnorm(-abs(avgt))
parampval<- 2*pnorm(-abs(correctedavgt))
plot(avgt,correctedavgt,bg=cols,pch=21,xlab="Original",ylab="With correction factor",xlim=c(-7,20),ylim=c(-7,20),cex=1.5)
abline(0,1)
abline(h=0,v=0,lty=2)
thirdhighest <- order(-avgt)[3]
arrows(avgt[thirdhighest]+3,correctedavgt[thirdhighest],x1=avgt[thirdhighest]+0.5,lwd=2)

avgcorrs[thirdhighest]

length(gsids[[thirdhighest]])

cf[thirdhighest]

# Permutations
library(matrixStats)

set.seed(1)
B <- 400 # takes a few minutes
null <- sapply(1:B,function(b){
 nullX<- sample(X)
 nullsvaX<-model.matrix(~nullX+svafit$sv) ##note that we are not recomupting the surrogate values. 
 nulllmfit <- lmFit(geneExpression,nullsvaX)
 nulltt<-nulllmfit$coef[,2]*sqrt(nulllmfit$df.residual)/(2*nulllmfit$sigma)
 nullavgt <- sapply(gsids,function(i) sqrt(length(i))*mean(nulltt[i]))
 return(nullavgt)
})
permavgt <- avgt/rowSds(null)
permpval<- rowMeans(abs(avgt) < abs(null))

plot(correctedavgt,permavgt,bg=cols,pch=21,xlab="Parametric z-score (with correction)",ylab="Permutation z-score",cex=1.5,ylim=c(-5,15),xlim=c(-5,15))
abline(0,1)
abline(h=0,v=0,lty=2)

tab <- data.frame(avgt=p.adjust(signif(2*pnorm(1-abs(avgt)),2),method="BH"),
                  correctedavgt=p.adjust(signif(2*pnorm(1-abs(correctedavgt)),2),method="BH"),
                  permutations=p.adjust(permpval,method="BH"))

# include only gene sets with 10 or more genes in comparison
tab <- tab[N>=10,]
tab <- cbind(signif(tab,2),apply(tab,2,rank))
tab<-tab[order(tab[,1]),]
tab <- tab[tab[,1]< 0.25,]
tab

# Gene set testing

library(GEOquery)
# g <- getGEO("GSE34313")
e <- g[[1]]
e$condition <- e$characteristics_ch1.2
levels(e$condition) <- c("treatment: dexamethasone for 24 hr","treatment: dexamethasone for 4 hr","treatment: none")
table(e$condition)
boxplot(exprs(e), range=0)
names(fData(e))
lvls <- c("treatment: none", "treatment: dexamethasone for 4 hr")
es <- e[,e$condition %in% lvls]
es$condition <- factor(es$condition, levels=lvls)

library(limma)
design <- model.matrix(~ es$condition)
fit <- lmFit(es, design=design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, genelist=fData(es)$GENE_SYMBOL)
tt

# ROAST method
# Immune response
idx <- grep("GO:0006955", fData(es)$GO_ID)
length(idx)
r1 <- roast(es, idx, design)
# ?roast
r1

# Testing multiple gene sets
library(org.Hs.eg.db)
org.Hs.egGO2EG
go2eg <- as.list(org.Hs.egGO2EG)
head(go2eg)
govector <- unlist(go2eg)
head(govector)
golengths <- sapply(go2eg, length)
head(golengths)
head(fData(es)$GENE)
idxvector <- match(govector, fData(es)$GENE)
head(idxvector)
table(is.na(idxvector))
idx <- split(idxvector, rep(names(go2eg), golengths))
go2eg[[1]]
fData(es)$GENE[idx[[1]]]

idxclean <- lapply(idx, function(x) x[!is.na(x)])
idxlengths <- sapply(idxclean, length)
head(idxlengths)
idxsub <- idxclean[idxlengths > 10]
length(idxsub)
# multiple roast
r2 <- mroast(es, idxsub, design)
head(r2)
r2 <- r2[order(r2$PValue.Mixed),]

library(GO.db)
columns(GO.db)
keytypes(GO.db)
GOTERM[[rownames(r2)[1]]]
r2tab <- select(GO.db, keys=rownames(r2)[1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]
r2 <- r2[order(r2$PValue),]
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Up"][1:10],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]
r2tab <- select(GO.db, keys=rownames(r2)[r2$Direction == "Down"][1:5],
                columns=c("GOID","TERM","DEFINITION"), 
                keytype="GOID")
r2tab[,1:2]



library(leukemiasEset)
library(preprocessCore)
library(genefilter)
library(qvalue)
library(org.Hs.eg.db)
data(leukemiasEset)
table(leukemiasEset$LeukemiaType)
boxplot(exprs(leukemiasEset)) # note different distributions across samples
lenorm = normalize.quantiles(exprs(leukemiasEset)) # quantile normalization
rownames(lenorm) = rownames(exprs(leukemiasEset)) # add row names
colnames(lenorm) = colnames(exprs(leukemiasEset)) # add column names
exprs(leukemiasEset) = lenorm # overwrite initial exprs
boxplot(exprs(leukemiasEset)) # all samples are quantile normalized

cml_normal = leukemiasEset[, leukemiasEset$LeukemiaType %in% c("CML", "NoL")]
cml_normal$LeukemiaType = factor(cml_normal$LeukemiaType)
cml_normal$LeukemiaType = relevel(cml_normal$LeukemiaType, "NoL")
factor(cml_normal$LeukemiaType)    # levels are now NoL, CML

library(org.Hs.eg.db)
symbols = mapIds(org.Hs.eg.db, rownames(cml_normal), keytype = "ENSEMBL", column="SYMBOL")
ids = mapIds(org.Hs.eg.db, rownames(cml_normal), keytype = "ENSEMBL", column="ENTREZID")
fData(cml_normal)$symbol = symbols
fData(cml_normal)$gene_id = ids
fData(cml_normal)
