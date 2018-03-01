#Gene table was downloaded from this link in table browser of UCSC:
#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=658679805_JkE5CAPKU7ctMv8j3q3l1EaPx3er&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=0&hgta_regionType=genome&position=chr12%3A109540092-109547123&hgta_outputType=primaryTable&hgta_outFileName=gene.test
#From the genes where removed random genes, mitocondrial genes and genes with negative start position after being extended. Genes were extended of 10 kb up stream and down stream the TSS

#Extend the genes of 10 kb from the TSS
#Extract
downstream=10000
upstream=10000

cat RefSeq.annotation.noRandom.noChrM.bed | awk '{ if ($4 == "+") { print $1"\t"$2-'$upstream'"\t" $2+'$downstream'"\t"$4"\t"$5} else if ($4 == "-") { print $1"\t"$3+'$upstream'"\t"$3-'$downstream'"\t"$4"\t"$5}}' > RefSeq.10kb.extended.bed


#####
#Random DMRs were obtained from bedtools sampling the same amount of DMRs as the number of the significant ones (6904) with the same lenght
 bedtools shuffle -incl background.file.with.all.DNA.regions.tested.bed -i HyperMe.DMRs.G45.TKO.bed -g mm10.chr.sizes


#Intersection of HyperMe DMRs with imprinted genes (normal and random)

bedtools intersect -a Hyper.Definitive.merged.bed -b Imprinted.Genes.bed -wb | cut -f 8 | sort -u | wc -l
bedtools intersect -a Hyper.Definitive.merged.Random.from.background.bed -b Imprinted.Genes.bed -wb | cut -f 8 | sort -u | wc -l


##Test
#Enrichment of HyperMe DMRs at Imprinted Genes:
# H0 we have high probability to find HyperMe DMRs at imprinted genes by chance
# H1 HyperMe DMRs are at imprinted genes not by chance
M <- as.table(rbind(c(37, 87), c(7, 117)))
dimnames(M) <- list(DMRs = c("Real", "Random"), Imprinted.Genes = c("yes","no"))
M
chisq.test(M)

###################################################################################################

#Intersection of HyperMe DMRs with all genes (normal)
bedtools intersect -a Hyper.Definitive.merged.bed -b RefSeq.10kb.extended.Ordered.bed -wb | cut -f 8 | sort -u | wc -l 
bedtools intersect -a Hyper.Definitive.merged.bed -b Imprinted.Genes.bed -wb | cut -f 8 | sort -u | wc -l


#Enrichment of HyperMe DMRs at Imprinted Genes:
# H0 HperMe DMRs are located at several random classes of genes
# H1 HyperMe DMRs are at not distribuited random in the genome but are significantly present at imprinted genes
M <- as.table(rbind(c(37, 3614), c(87, 32124)))
dimnames(M) <- list(DMRs =c("Imprinted.yes", "Imprinted.no"), Genes.at.DMRs = c("Imprinted.yes","Imprinted.no"))
M
chisq.test(M)

####################################################################################################

#Intersection of Random HyperMe DMRs with all genes (random)
bedtools intersect -a Hyper.Definitive.merged.Random.from.background.bed -b RefSeq.10kb.extended.Ordered.bed -wb | cut -f 8 | sort -u | wc -l  
bedtools intersect -a Hyper.Definitive.merged.Random.from.background.bed -b Imprinted.Genes.bed -wb | cut -f 8 | sort -u | wc -l 


#Enrichment of HyperMe DMRs at Imprinted Genes:
# H0 Random HperMe DMRs are located at imprinted genes
# H1 Random HyperMe DMRs are not imprinted genes
M <- as.table(rbind(c(7, 2624), c(117, 33114)))
dimnames(M) <- list(Random.DMRs =c("Imprinted.yes", "Imprinted.no"), Genes.at.DMRs = c("Imprinted.yes","Imprinted.no"))
M
chisq.test(M)
