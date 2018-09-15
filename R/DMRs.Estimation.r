library(methylKit)

#This part to read the cytosine report files and select CpG with a given coverage threshold is taken from the methylkit developer Altuna: 
#https://gist.github.com/al2na/4839e615e2401d73fe51

readBismarkCytosineReport<-function(location,sample.id,assembly="unknown",treatment,
                                    context="CpG",min.cov=10){
  if(length(location)>1){
    stopifnot(length(location)==length(sample.id),
              length(location)==length(treatment))
  }
  
  result=list()
  for(i in 1:length(location)){
    df=fread.gzipped(location[[i]],data.table=FALSE)
    
    # remove low coverage stuff
    df=df[ (df[,4]+df[,5]) >= min.cov ,]
    
    
    
    
    # make the object (arrange columns of df), put it in a list
    result[[i]]= new("methylRaw",
                     data.frame(chr=df[,1],start=df[,2],end=df[,2],
                                strand=df[,3],coverage=(df[,4]+df[,5]),
                                numCs=df[,4],numTs=df[,5]),
                     sample.id=sample.id[[i]],
                     assembly=assembly,context=context,resolution="base"
    )
  }
  
  if(length(result) == 1){
    return(result[[1]])
  }else{
    
    new("methylRawList",result,treatment=treatment)
  }
}

# reads gzipped files,
fread.gzipped<-function(filepath,...){
  require(R.utils)
  require(data.table)
  
  
  
  # decompress first, fread can't read gzipped files
  if (R.utils::isGzipped(filepath)){
    
    if(.Platform$OS.type == "unix") {
      filepath=paste("zcat",filepath)
    } else {
      filepath <- R.utils::gunzip(filepath,temporary = FALSE, overwrite = TRUE,
                                  remove = FALSE)
    }
    
    
  }
  
  ## Read in the file
  fread(filepath,...)
  
}


########################################
###Read your Sample and Call DMRs#######
########################################
setwd("/gpfs/fs6/jgu-cbdm/andradeLab/scratch/tandrean/Data/WGBS/Gadd45.TKO/Gadd45.tko.Quality.Filter/Extraction/MethylKit")

#Files are the output of Bismarck i.e. the Cytosine Report Files (example Gadd45.tko2.CpG_report.txt) with a ToT CG of 43841737 in each file
file.list <- list("control1.myCpG.gz","control2.myCpG.gz","test2.myCpG.gz","test3.myCpG.gz")
myobj=readBismarkCytosineReport(file.list,sample.id=list("ctrl1","ctrl2","test2","test3"),assembly="mm10",treatment=c(0,0,1,1))
tiles <- tileMethylCounts(myobj,cov.bases = 2, win.size = 100,step.size = 100)
meth=unite(tiles,destrand=FALSE)
pdf('Samples.ctrl1.ctrl2.test2.test3.Correlation.Tiles.100.pdf')
getCorrelation(meth, plot = T)
dev.off()

pdf('Samples.ctrl1.ctrl2.test2.test3.PCA.Tiles.100.pdf')
PCASamples(meth)
dev.off()

covariates=data.frame(experiment=c(0,1,0,0))
myDiff <- calculateDiffMeth(meth,covariates=covariates,overdispersion=c("MN"),effect=c("predicted"),test=c("Chisq"),num.cores=4)
myDiff25p.hyper <- getMethylDiff(myDiff, difference = 30,qvalue = 0.05, type = "hyper")
myDiff25p.hypo <- getMethylDiff(myDiff, difference = 30,qvalue = 0.05, type = "hypo")
Hyper <- getData(myDiff25p.hyper)
Hypo <- getData(myDiff25p.hypo)



write.table(myDiff,"Background.Gadd45.TKO.cov.10.Tiles.100.delta.30.txt",quote=FALSE,col.names=T,row.names=F,sep="\t") ##For Backgroung
write.table(myDiff25p.hyper,"Hyper.DMRs.Gadd45.TKO.cov.10.Tiles.100.delta.30.FDR.5%.txt",quote=FALSE,col.names=T,row.names=F,sep="\t")
write.table(myDiff25p.hypo,"Hypo.DMRs.Gadd45.TKO.cov.10.Tiles.100.delta.30.FDR.5%.txt",quote=FALSE,col.names=T,row.names=F,sep="\t")

