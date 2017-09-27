source("https://bioconductor.org/biocLite.R")
biocLite("MethylSeekR")
library(MethylSeekR)
source("http://www.bioconductor.org/biocLite.R")
biocLite("BSgenome")
library(BSgenome)
library(data.table)
available.genomes()
biocLite("BSgenome.Mmusculus.UCSC.mm10")

####################
#Analysis Neil DKO##
####################
setwd("/media/tandrean/Elements/PhD/Whole.Genome.Bisulphite/Script/")

#Granges Object Wild Type
mESC.WT <- fread("control1.myCpG.chr1.methylseekeR.txt")
colnames(mESC.WT)[1] <- "chr"
colnames(mESC.WT)[2] <- "start"
colnames(mESC.WT)[3] <- "end"
colnames(mESC.WT)[4] <- "strand"
colnames(mESC.WT)[5] <- "T"
colnames(mESC.WT)[6] <- "M"
mESC.WT$T <- as.numeric(mESC.WT$T)
mESC.WT$M <- as.numeric(mESC.WT$M)
Gr.obj.mESC <- GRanges(seqnames=mESC.WT$chr,ranges=IRanges(mESC.WT$start,mESC.WT$end),strand=TRUE,T=mESC.WT$T,M=mESC.WT$M)


#Granges Object Neil DKO
Neil.DKO <- fread("test1.myCpG.chr1.methylseekeR.txt")
colnames(Neil.DKO)[1] <- "chr"
colnames(Neil.DKO)[2] <- "start"
colnames(Neil.DKO)[3] <- "end"
colnames(Neil.DKO)[4] <- "strand"
colnames(Neil.DKO)[5] <- "T"
colnames(Neil.DKO)[6] <- "M"
Neil.DKO$T <- as.numeric(Neil.DKO$T)
Neil.DKO$M <- as.numeric(Neil.DKO$M)
Gr.obj.Neil.DKO <- GRanges(seqnames=Neil.DKO$chr,ranges=IRanges(Neil.DKO$start,Neil.DKO$end),strand=TRUE,T=Neil.DKO$T,M=Neil.DKO$M)

#Granges object SNPs
SNPs <- fread("SNP.chr1.mm10.bed")
colnames(SNPs)[1] <- "chr"
colnames(SNPs)[2] <- "start"
colnames(SNPs)[3] <- "end"
colnames(SNPs)[4] <- "strand"
Gr.obj.SNPs <- GRanges(seqnames=SNPs$chr,ranges=IRanges(SNPs$start,SNPs$end),strand=TRUE)

library("BSgenome.Mmusculus.UCSC.mm10")
sLengths=seqlengths(Mmusculus)
set.seed(123)

#Remove SNPs
Gr.obj.mESC <- removeSNPs(Gr.obj.mESC,Gr.obj.SNPs)
Gr.obj.Neil.DKO <- removeSNPs(Gr.obj.Neil.DKO,Gr.obj.SNPs)

plotAlphaDistributionOneChr(m=Gr.obj.mESC, chr.sel="chr1", num.cores=1)
plotAlphaDistributionOneChr(m=Gr.obj.Neil.DKO, chr.sel="chr1", num.cores=1)

library(parallel)
detectCores()

#Detect PMD segments in mESC
PMDsegments.gr.WT <- segmentPMDs(m=Gr.obj.mESC, chr.sel="chr1", seqLengths=sLengths, num.cores=1)
plotAlphaDistributionOneChr(m=subsetByOverlaps(Gr.obj.mESC, PMDsegments.gr.WT[values(PMDsegments.gr.WT)$type=="notPMD"]), chr.sel="chr1", num.cores=1)
plotPMDSegmentation(m=Gr.obj.mESC, numRegions = 5, pdfFilename = "5.Regions.wt.pdf" , segs=PMDsegments.gr.WT)
write.table(PMDsegments.gr.WT,"segment.PMD.wt.txt",quote=F,col.names = T,row.names = F, sep="\t")

#Detect PMD segments in Neil DKO
PMDsegments.gr.Neil.DKO <- segmentPMDs(m=Gr.obj.Neil.DKO, chr.sel="chr1", seqLengths=sLengths, num.cores=1)
plotAlphaDistributionOneChr(m=subsetByOverlaps(Gr.obj.Neil.DKO, PMDsegments.gr.Neil.DKO[values(PMDsegments.gr.Neil.DKO)$type=="notPMD"]), chr.sel="chr1", num.cores=1)
plotPMDSegmentation(m=Gr.obj.Neil.DKO,  numRegions = 20, pdfFilename = "20.Regions.Neil.DKO.pdf", segs=PMDsegments.gr.Neil.DKO)
write.table(PMDsegments.gr.Neil.DKO,"segment.PMD.Neil.DKO.txt",quote=F,col.names = T,row.names = F, sep="\t")

library(rtracklayer)

#Retrieve CpG islands
session <- browserSession()
genome(session) <- "mm10"
query <- ucscTableQuery(session, "cpgIslandExt")
CpGislands.gr <- track(query)
genome(CpGislands.gr) <- NA
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))

#Detect UMR and LMR in mESC WT
stats.wt <- calculateFDRs(m=Gr.obj.mESC, CGIs=CpGislands.gr,PMDs=PMDsegments.gr.WT, num.cores=1)
stats.wt
FDR.cutoff <- 5
m.sel <- 0.5
n.sel=as.integer(names(stats.wt$FDRs[as.character(m.sel), ] [stats.wt$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
n.sel
UMRLMRsegments.gr.WT <- segmentUMRsLMRs(m=Gr.obj.mESC, meth.cutoff=m.sel, nCpG.cutoff=n.sel, PMDs=PMDsegments.gr.WT, num.cores=1, myGenomeSeq=Mmusculus, seqLengths=sLengths)
plotFinalSegmentation(m=Gr.obj.mESC, segs=UMRLMRsegments.gr.WT, PMDs=PMDsegments.gr.WT,meth.cutoff=m.sel,numRegions = 20, pdfFilename = "20.Regions.UMR.LMR.pdf")
write.table(UMRLMRsegments.gr.WT,"segment.UMR.wt.txt",quote=F,col.names = T,row.names = F, sep="\t")

#Detect UMR and LMR in Neil DKO
stats.neil.dko <- calculateFDRs(m=Gr.obj.Neil.DKO, CGIs=CpGislands.gr,PMDs=PMDsegments.gr.Neil.DKO, num.cores=1)
stats.neil.dko
FDR.cutoff <- 5
m.sel <- 0.5
n.sel=as.integer(names(stats.neil.dko$FDRs[as.character(m.sel), ] [stats.neil.dko$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
n.sel
UMRLMRsegments.gr.neil.dko <- segmentUMRsLMRs(m=Gr.obj.Neil.DKO, meth.cutoff=m.sel, nCpG.cutoff=n.sel, PMDs=PMDsegments.gr.Neil.DKO, num.cores=1, myGenomeSeq=Mmusculus, seqLengths=sLengths)
write.table(UMRLMRsegments.gr.neil.dko,"segment.UMR.Neil.DKO.txt",quote=F,col.names = T,row.names = F, sep="\t")
plotFinalSegmentation(m=Gr.obj.Neil.DKO, segs=UMRLMRsegments.gr.neil.dko, PMDs=PMDsegments.gr.Neil.DKO,meth.cutoff=m.sel,numRegions = 20, pdfFilename = "20.Regions.UMR.LMR.Neil.DKO.pdf")


#Segment Length Histogram
UMR.Length.mESC <- read.table("UMR.mESC.wt.txt")
UMR.Length.Neil.DKO <- read.table("UMR.Neil.DKO.txt")


UMR.Length.mESC$length <- UMR.Length.mESC$V4
UMR.Length.Neil.DKO$lenght <- UMR.Length.Neil.DKO$V4

UMR.Length.mESC <- as.data.frame(UMR.Length.mESC[,13])
UMR.Length.Neil.DKO <- as.data.frame(UMR.Length.Neil.DKO[,13])
names(UMR.Length.mESC)[1]<-"length"
names(UMR.Length.Neil.DKO)[1]<-"length"

UMR.Length.mESC$Class <- 'mESC'
UMR.Length.Neil.DKO$Class <- 'Neil.DKO'

UMR.Lengths <- rbind(UMR.Length.mESC,UMR.Length.Neil.DKO)
dim(UMR.Lengths)
library(ggplot2)
p <- ggplot(UMR.Lengths, aes(length, fill = Class)) + geom_density(alpha = 0.2) 
p +labs(title="Comparison of the Length UMR in Chr1",x ="Length UMR")
mean(UMR.Length.mESC$length)
mean(UMR.Length.Neil.DKO$length)



LMR.Length.mESC <- read.table("LMR.mESC.wt.txt")
LMR.Length.Neil.DKO <- read.table("LMR.Neil.DKO.txt")
LMR.Length.mESC$length <- LMR.Length.mESC$V4
LMR.Length.Neil.DKO$lenght <- LMR.Length.Neil.DKO$V4

LMR.Length.mESC <- as.data.frame(LMR.Length.mESC[,13])
LMR.Length.Neil.DKO <- as.data.frame(LMR.Length.Neil.DKO[,13])
names(LMR.Length.mESC)[1]<-"length"
names(LMR.Length.Neil.DKO)[1]<-"length"

LMR.Length.mESC$Class <- 'mESC'
LMR.Length.Neil.DKO$Class <- 'Neil.DKO'

LMR.Lengths <- rbind(LMR.Length.mESC,LMR.Length.Neil.DKO)
dim(LMR.Lengths)
library(ggplot2)
p <- ggplot(LMR.Lengths, aes(length, fill = Class)) + geom_density(alpha = 0.2) 
p +labs(title="Comparison of the Length LMR in Chr1",x ="Length LMR")
mean(LMR.Length.mESC$length)
mean(LMR.Length.Neil.DKO$length)


PMD.Length.mESC <- read.table("PMD.mESC.wt.txt")
PMD.Length.Neil.DKO <- read.table("PMD.Neil.DKO.txt")

PMD.Length.mESC$length <- PMD.Length.mESC$V4
PMD.Length.Neil.DKO$lenght <- PMD.Length.Neil.DKO$V4

head(PMD.Length.mESC)

PMD.Length.mESC <- as.data.frame(PMD.Length.mESC[,8])
PMD.Length.Neil.DKO <- as.data.frame(PMD.Length.Neil.DKO[,8])
names(PMD.Length.mESC)[1]<-"length"
names(PMD.Length.Neil.DKO)[1]<-"length"

PMD.Length.mESC$Class <- 'mESC'
PMD.Length.Neil.DKO$Class <- 'Neil.DKO'

PMD.Lengths <- rbind(PMD.Length.mESC,PMD.Length.Neil.DKO)
dim(PMD.Lengths)
library(ggplot2)
p <- ggplot(PMD.Lengths, aes(length, fill = Class)) + geom_density(alpha = 0.2) 
p +labs(title="Comparison of the Length PMD in Chr1",x ="Length PMD")
mean(PMD.Length.mESC$length)
mean(PMD.Length.Neil.DKO$length)

source("https://bioconductor.org/biocLite.R")
biocLite("genomation")
library(genomation)
library(methylKit)
setwd("/media/tandrean/Elements/PhD/Whole.Genome.Bisulphite/DMRs/Sample5.Vs.Sample6.7.8.Pooled")
hyper <- read.table("Hyper.Methytlated.DMRs.1000.bp.tiles.coverage.min.10.Clone5.Vs.Pooled.TKO.p.val.0.1.delta.0.15.txt",header=T)
myDiff25p <- as(hyper,"GRanges")
gene.obj=readTranscriptFeatures("mm10.refseq.genes.bed")
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
anno
cpg.obj=readFeatureFlank("cpgIslandExt.txt.bed",feature.flank.name=c("CpGi","shores"))
cpg.obj
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),cpg.obj$CpGi,cpg.obj$shores,feature.name="CpGi",flank.name="shores")
diffCpGann
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
getAssociationWithTSS(diffAnn)
diffAnn


setwd("/media/tandrean/Elements/PhD/Whole.Genome.Bisulphite/Script")
UMR.mESC <- read.table("UMR.mESC.wt.txt")
UMR.Neil.DKO <- read.table("UMR.Neil.DKO.txt")
head(UMR.mESC)
hist(UMR.mESC$V11)
hist(UMR.Neil.DKO$V11)


LMR.mESC <- read.table("LMR.mESC.wt.txt")
LMR.Neil.DKO <- read.table("LMR.Neil.DKO.txt")

LMR.mESC$meth <- LMR.mESC$V11
LMR.Neil.DKO$meth <-LMR.Neil.DKO$V11

LMR.mESC.meth <- as.data.frame(LMR.mESC[,13])
LMR.Neil.DKO.meth <- as.data.frame(LMR.Neil.DKO[,13])

names(LMR.mESC.meth)[1]<-"meth"
names(LMR.Neil.DKO.meth)[1]<-"meth"

LMR.mESC.meth$Class <- 'mESC'
LMR.Neil.DKO.meth$Class <- 'Neil.DKO'

LMR.meth <- rbind(LMR.mESC.meth,LMR.Neil.DKO.meth)

head(LMR.meth)
library(ggplot2)
library(methylKit)
p <- ggplot(data = LMR.meth, aes(x=Class, y=meth)) + geom_boxplot(aes(fill=Class))+ ggtitle("LMR values in WT and Neil DKO (chr1)")
p

setwd("/media/tandrean/Elements/PhD/Whole.Genome.Bisulphite/Script")
#LMR.Neil.DKO.in.mESC <- read.table("LMR.Neil.DKO.in.mESC.txt")
#boxplot(LMR.Neil.DKO.in.mESC$V23,LMR.Neil.DKO.in.mESC$V11)

UMR.Length.mESC <- read.table("segment.UMR.wt.txt",header=T)
UMR.Length.Neil.DKO <- read.table("segment.UMR.Neil.DKO.txt",header=T)

UMR.Length.mESC.gr <- makeGRangesFromDataFrame(UMR.Length.mESC)
UMR.Length.Neil.DKO.gr <- makeGRangesFromDataFrame(UMR.Length.Neil.DKO)

fmr.Neil.DKO <- gaps(UMR.Length.Neil.DKO.gr)
fmr.mESC <- gaps(UMR.Length.mESC.gr)

head(UMR.Length.Neil.DKO)
head(fmr.Neil.DKO)


dim(UMR.Length.Neil.DKO)
dim(UMR.Length.mESC)

length(intersect(fmr.Neil.DKO,fmr.mESC))


fmr.Neil.DKO <- fmr.Neil.DKO[chrom(fmr.Neil.DKO)=="chr1"]
fmr.mESC <- fmr.mESC[chrom(fmr.mESC)=="chr1"]

library(data.table)
mESC.WT <- fread("control1.myCpG.chr1.methylseekeR.txt")

colnames(mESC.WT)[1] <- "chr"
colnames(mESC.WT)[2] <- "start"
colnames(mESC.WT)[3] <- "end"
colnames(mESC.WT)[4] <- "strand"
colnames(mESC.WT)[5] <- "T"
colnames(mESC.WT)[6] <- "M"
mESC.WT$T <- as.numeric(mESC.WT$T)
mESC.WT$M <- as.numeric(mESC.WT$M)
mESC.WT$Meth <- as.numeric(mESC.WT$M/(mESC.WT$T+mESC.WT$M))

Gr.obj.mESC <- GRanges(seqnames=mESC.WT$chr,ranges=IRanges(mESC.WT$start,mESC.WT$end),strand=TRUE,Meth=mESC.WT$Meth)



Neil.DKO <- fread("test1.myCpG.chr1.methylseekeR.txt")
colnames(Neil.DKO)[1] <- "chr"
colnames(Neil.DKO)[2] <- "start"
colnames(Neil.DKO)[3] <- "end"
colnames(Neil.DKO)[4] <- "strand"
colnames(Neil.DKO)[5] <- "T"
colnames(Neil.DKO)[6] <- "M"
Neil.DKO$T <- as.numeric(Neil.DKO$T)
Neil.DKO$M <- as.numeric(Neil.DKO$M)
Neil.DKO$Meth <- as.numeric(Neil.DKO$M/(Neil.DKO$T+Neil.DKO$M))

Gr.obj.Neil.DKO <- GRanges(seqnames=Neil.DKO$chr,ranges=IRanges(Neil.DKO$start,Neil.DKO$end),strand=TRUE,Meth=Neil.DKO$Meth)

head(Gr.obj.Neil.DKO)
head(Gr.obj.mESC)




Sub.Neil.DKO <- subsetByOverlaps(Gr.obj.Neil.DKO,fmr.Neil.DKO)
Sub.mESC <- subsetByOverlaps(Gr.obj.mESC,fmr.mESC)


head(Gr.obj.mESC)

LMR.Neil.DKO <- subset(UMR.Length.Neil.DKO,type=LMR)
LMR.mESC <- subset(UMR.Length.mESC,type=LMR)

LMR.mESC <- GRanges(LMR.mESC)
LMR.Neil.DKO <- GRanges(LMR.Neil.DKO)

Sub.Neil.DKO.LMR.meth <- subsetByOverlaps(Gr.obj.Neil.DKO,LMR.Neil.DKO)
Sub.mESC.LMR.meth <- subsetByOverlaps(Gr.obj.mESC,LMR.mESC)
summary(Sub.mESC.LMR.meth$Meth)
summary(Sub.Neil.DKO.LMR.meth$Meth,Sub.mESC.LMR.meth$Meth)


#length FMR

fmr.Neil.DKO.df <- as.data.frame(fmr.Neil.DKO)
fmr.mESC.df <- as.data.frame(fmr.mESC)
fmr.Neil.DKO.df$length <- fmr.Neil.DKO.df$end- fmr.Neil.DKO.df$start
fmr.mESC.df$length <- fmr.mESC.df$end- fmr.mESC.df$start
hist(fmr.mESC.df$length)
hist(fmr.Neil.DKO.df$length)

head(fmr.mESC.df)
fmr.mESC.df <- as.data.frame(fmr.mESC.df[,6])
fmr.Neil.DKO.df <- as.data.frame(fmr.Neil.DKO.df[,6])


names(fmr.mESC.df)[1]<-"length"
names(fmr.Neil.DKO.df)[1]<-"length"

head(fmr.mESC.df)

fmr.mESC.df$Class <- 'mESC'
fmr.Neil.DKO.df$Class <- 'Neil.DKO'

fmr.Lengths <- rbind(fmr.mESC.df,fmr.Neil.DKO.df)
tail(fmr.Lengths)

library(ggplot2)
p <- ggplot(data = fmr.Lengths, aes(x=Class, y=length)) + geom_boxplot(aes(fill=Class))+ ggtitle("CpG Islands Meth Values at Promoters of mESC Vs Pooled Neil DKO (meth values > 0.1)")
p

p +labs(title="Comparison of the Length FMR in Chr1",x ="Length FMR")
mean(UMR.Length.mESC$length)
mean(UMR.Length.Neil.DKO$length)

