

#pre-processing the curatedOvarianData package (this was completed in 2015.)
#Note that the processing functions had to be a bit different than the breast
#cancer datasets because the curatedOvarianData authors did not provide quite
#as many R-level data access functions.
library("Biobase")
library("Coincide")
library("curatedOvarianData")
datasetNames <- data(package="curatedOvarianData")
datasetNames <- datasetNames[3]$results[,"Item"]

#NOTE: I found it easier to just untar the package myself and access the 
#data files myself for pre-processing.
#Change these paths to your own:
#first one is where downloaded curatedOvarianData package - need /data file:
local_download_path <- "~/R/x86_64-pc-linux-gnu-library/curatedOvarianData/data/" 
saveDirGlobal <- "./ovarian_analysis/"
saveDir_20 <- "./ovarian_analysis/withTop20Genes"
saveDir_no20 <- "./ovarian_analysis/notTop20Genes"
outputFile <- "./ovarian_analysis/ovarian_proc_messages.txt"

#odd way to get a dataset list...I just downloaded the source file
#and pull out their code.
#hmmm...is this the case for all genes?? WTF?
expandProbesets <- function (eset, sep = "///"){
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  eset <- eset[order(sapply(x, length)), ]
  x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
  idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
  xx <- !duplicated(unlist(x))
  idx <- idx[xx]
  x <- unlist(x)[xx]
  eset <- eset[idx, ]
  #Katie's notes: featureNames is from Biobase, for eset objects.
  featureNames(eset) <- x
  eset
}


#locally:
#data(list=data(package=package.name)[[3]][,3])
#strEsets <- ls(pattern="^.*_eset$")
#ON SERVER: just download the tarbell to get package, code
#cd ./R/x86_64-redhat-linux-gnu-library
#sudo wget http://www.bioconductor.org/packages/release/data/experiment/src/contrib/curatedOvarianData_1.3.5.tar.gz
#tar xvzf curatedOvarianData_1.3.5.tar.gz
 strEsets <- list.files(local_download_path)
# #datalist is a text file...shouldn't be any other text files in this directory, only .rda files.

 #only want rda datasets
  strEsets <- strEsets[grep("rda",strEsets)]
 
# library("limma")
# names(strEsets) <- strsplit2(strEsets,split=".rda")[,1]
# #now load all of this data. equivalent of data() call above if can install package.
# #is having trouble reading from a connection?? what's up...
# #WTF..may have to manually step through this one..
# #weird...now suddently dataset 16 isn't working? why only THIS one?
# #makes me wonder if there is a space issue on the server.
# #I had to re-download the datasets...
for (strEset in strEsets){
  #loading each dataset:
  load(paste0(local_download_path,strEset))
  
}
#now call strEset again like above:
strEsets <- ls(pattern="^.*_eset$")

##now back to common code for server or local computer:

#UBD///GABBR1 in final output...
esets <- list()
for (strEset in strEsets){
  
  eset <- get(strEset)
  ##Split out rows without unique gene name
  eset <- expandProbesets(eset, sep = "///")
  if(length(grep("///",rownames(exprs(eset))))>0){
    
    stop("Not collapsing probes correcting.")
  }
  esets[[strEset]] <- eset
  
}



#pData(esets) already contains survival data.

#remove the extra TCGA datasets
#TCGA.RNASeqV2_eset                           Integrated genomic analyses of ovarian carcinoma.
#TCGA.mirna.8x15kv2_eset                      Integrated genomic analyses of ovarian carcinoma.
#TCGA_eset                                    Integrated genomic analyses of ovarian carcinoma.
indicesRemove <- na.omit(match(c("TCGA.RNASeqV2_eset","TCGA.mirna.8x15kv2_eset"),names(esets)))
esets <- esets[-indicesRemove]

#save for easier access next time:
#saveRDS(esets, file="./curatedOvarianData_esetList.rds",compress=TRUE)
#remove any healthy samples. do this for ALL studies in case some have healthy tissue.
#for example: TCGA has some normal samples
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"sample_type"])
indicesRemove <- c()

# work around bug that causes error when validating because assaydata is an environment
esets2 <- lapply(esets, function(x) x[x@phenoData$sample_type == "tumor",])
esets3 <- mapply(function(x, y) x[, y], esets, esets2)
# for(e in 1:length(esets)){
  # esets[[e]] <- esets[[e]][esets[[e]]@phenoData$sample_type =="tumor",]
  # #choose only primary tumor samples.
  # pDat <- pData(esets[[e]])[which(pData(esets[[e]])["sample_type"]=="tumor"),]
  # expr <- exprs(esets[[e]])[ ,which(pData(esets[[e]])["sample_type"]=="tumor")]
  # 
  # #if no samples left: remove
  # if(ncol(expr)==0){
  #   
  #   indicesRemove <- append(indicesRemove,e)
  #   
  # }else{
  #   
  #   exprs(esets[[e]]) <- expr
  #   pData(esets[[e]]) <- pDat
  #   protocolData(esets[[e]])@data <- data.frame(row.names=colnames(exprs(esets[[e]])))
    #check: make sure is still a valid object
    # validObject(esets[[e]])
# }
  
message(paste("dataset indices to remove : ",indicesRemove))

esets <- esets3
#gut check: in TCGA, should be only tumor samples now!
if(any(unique(pData(esets3[["TCGA_eset"]])[,"sample_type"]) != "tumor")){
  
  stop("\nNormal filtering appeared to not work.")
  
}

#check? any duplicated samples?
#it looks like unique_patient_ID ids duplicated samples, NOT the actual expression colnames?
#a little annoying - re-label the colnames to be this unique ID, if duplicated samples are found
#but esets can't have duplicated column names.
#in the end, looks like no samples were duplicated? 
for(e in 1:length(esets)){
  
  if(!all(is.na(pData(esets[[e]])[ ,"unique_patient_ID"]))){
    #for certain datasets: this is actually all NA values. so don't use then!
    if(any(duplicated(pData(esets[[e]])[,"unique_patient_ID"]))){
      
      #if(any(duplicated(colnames(exprs(esets[[e]]))))){
      
      
      message(paste0("Duplicated samples in study ",names(esets)[e]))
      
    }
    
  }
  
}

#something to think about later: filter by primary site? 
head(pData(esets[["TCGA_eset"]])[,"primarysite"])


#}
#want very lowly varying genes to be removed for Combat.
#featureDataFieldName: must match the gene symbol column title for these esets.
esets <- procExprSetList(exprSetList=esets,outputFileDirectory=saveDirGlobal,
                                  minVar=.001,featureDataFieldName="gene",uniquePDataID="unique_patient_ID")

#look at batches - looks like the only compress=TRUE batches (besides machine runs)
#are in the TCGA dataset.
for(e in 1:length(esets)){ 
  
  if(!all(is.na(pData(esets[[e]])[,"batch"]))){
    
    cat("\n",names(esets)[e]," contains labeled batches")
    cat("\n",unique(pData(esets[[e]])[,"batch"],"\n"))
    
  }
  
}
# GSE14764_eset  contains labeled batches
# 2005-01-25 2005-01-21 2005-01-28 2005-01-26 2005-03-02 2004-09-30 2004-10-01 2004-09-29 2006-08-11 2006-07-26 2006-08-21 2006-07-27 2006-08-18 2006-07-28 2006-08-19

# GSE18520_eset  contains labeled batches
# 2004-03-12 2004-04-08 2004-04-09 2004-07-20 2004-08-12 2004-08-13

# GSE19829.GPL570_eset  contains labeled batches
# 2009-08-14

# GSE19829.GPL8300_eset  contains labeled batches
# 2002-08-20 2003-09-09 2001-09-14 2003-09-18 2001-12-14

# GSE20565_eset  contains labeled batches
# 2006-06-01 2006-06-27 2006-06-30 2006-06-28 2006-06-29 2006-07-20 2008-03-06

# GSE2109_eset  contains labeled batches
# 2006-01-26 2006-03-06 2006-02-28 2006-02-07 2006-03-14 2006-01-24 2006-06-08 2006-04-18 2006-05-16 2006-04-20 2006-09-14 2006-09-12 2006-07-26 2006-07-28 2006-10-24 2006-11-09 2006-11-30 2006-10-31 2006-12-07 2006-10-10 2006-11-21 2007-03-07 2007-03-15 2007-01-12 2007-02-09 2007-03-09 2007-05-03 2007-05-01 2007-05-15 2007-05-18 2007-05-30 2007-06-12 2007-07-27 2007-09-07 2007-09-11 2007-09-05 2007-09-12 2008-02-15 2008-03-04 2008-02-21 2008-02-27 2008-05-16 2008-05-23 2008-05-13 2004-12-03 2004-12-04 2004-12-07 2005-03-22 2005-03-19 2005-03-11 2005-03-16 2005-03-17 2005-03-15 2005-03-10 2005-03-03 2005-02-11 2005-06-01 2005-04-26 2005-06-03 2005-04-29 2005-05-10 2005-06-08 2005-04-13 2005-06-17 2005-08-11 2005-09-07 2005-09-09 2005-08-05 2005-08-09 2005-09-13 2005-11-15 2005-12-02 2005-11-02 2005-11-04 2005-11-18

# GSE26193_eset  contains labeled batches
# 2006-06-01 2006-06-27 2009-03-18 2009-03-19 2006-06-30 2006-06-28 2006-06-29 2006-07-20 2008-03-06

# GSE26712_eset  contains labeled batches
# 2003-11-04 2003-11-06 2003-11-20 2003-11-07 2004-04-21 2003-11-21 2003-12-16 2004-04-20 2003-12-23 2003-12-24 2004-04-27 2003-11-05 2004-09-28 2005-07-27

# GSE30161_eset  contains labeled batches
# 2009-10-07 2009-10-09 2009-10-20 2009-10-08

# GSE44104_eset  contains labeled batches
# 2010-12-14 2010-10-14 2010-09-07 2010-12-10 2010-09-08

# GSE6008_eset  contains labeled batches
# 2002-08-23 2002-08-27 2002-08-13 2002-04-04 2002-04-09 2002-08-15 2002-08-22 2002-08-28 2002-08-29 2002-08-30 2002-04-10 2002-04-03 2002-04-12 2002-09-11

# GSE6822_eset  contains labeled batches
# 2000-12-21 2001-05-03 2001-05-29 2001-06-12 2001-09-25 2001-09-26 2001-09-27 2002-02-14 2002-04-17 2002-04-18 2002-07-18 2002-07-24 2002-10-20 2002-10-30 2002-11-13

# GSE9891_eset  contains labeled batches
# 2004-12-03 2005-01-12 2005-01-17 2005-02-21 2004-12-23 2005-05-05 2005-07-15 2005-07-06 2005-06-16 2005-06-08 2005-06-06 2005-01-24 2005-06-17 2006-04-13 2005-12-21 2005-11-09 2006-04-06 2005-08-24 2005-10-05 2005-08-18 2005-11-04 2005-08-05 2005-09-09 2005-01-31 2005-06-02 2005-12-15 2006-05-03 2006-01-31 2006-01-20 2005-06-24 2005-08-26 2005-07-29 2005-10-28 2005-07-20 2006-02-08 2006-06-22 2006-07-19 2006-04-28 2006-04-12 2005-11-11 2005-09-21 2006-06-06 2005-10-26 2006-02-28 2006-04-05 2005-08-03 2005-09-14 2005-05-27 2005-05-30 2005-09-16 2005-05-09 2005-03-17 2006-06-07 2005-05-25 2005-11-23 2006-07-07

# PMID15897565_eset  contains labeled batches
# 2002-09-20 2002-10-23 2002-11-12 2002-12-16 2002-12-21 2003-01-03 2003-05-30 2003-07-02

# PMID17290060_eset  contains labeled batches
# 2002-09-20 2002-10-23 2002-11-12 2002-12-16 2002-12-21 2003-01-03 2003-05-30 2004-06-22 2004-06-23 2004-05-27 2004-05-18 2004-05-21 2004-03-09 2004-03-16 2004-04-20

# PMID19318476_eset  contains labeled batches
# 2004-05-21 2004-06-22 2004-05-27 2004-05-18 2004-04-20 2004-06-23 2004-03-09 2004-03-16

# TCGA_eset  contains labeled batches
# 12 NA 15 9 13 24 22 21 17 40 18 19 11 14 27



#only split up TCGA into batches
#na.omit: if any NA batch variables, then unique will return one "NA" 
TCGAbatches <- na.omit(unique(pData(esets[["TCGA_eset"]])[,"batch"]))
#any samples with NA batch? these will be thrown out.
which(is.na(pData(esets[["TCGA_eset"]])[,"batch"]))
# [1] 29

#let's see if batches seem to have uneven distributions of known subtypes
#and/or grade...that will make me hesitant to use batch correction.
#hmm...so most are serious?
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"histological_type"])
#    ser
# 9   37
# 11  37
# 12  45
# 13  47
# 14  46
# 15  22
# 17  47
# 18  47
# 19  47
# 21  46
# 22  47
# 24  46
# 27   6
# 40  47
unique(pData(esets[["TCGA_eset"]])[,"histological_type"])
# [1] "ser" NA   
#yep...only 2 are NA, rest are serous
length(which(is.na(pData(esets[["TCGA_eset"]])[,"histological_type"])))
# [1] 2
#this matches the Nature paper saying the samples are only serous: http://www.nature.com/nature/journal/v474/n7353/full/nature10166.html

#and predominantly late-stage for all batches
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"summarystage"])
#    early late
# 9      1   36
# 11     0   37
# 12     2   42
# 13     3   44
# 14     0   46
# 15     0   22
# 17     4   43
# 18     4   43
# 19     7   40
# 21     2   42
# 22     9   38
# 24     7   38
# 27     1    5
# 40     3   43

#same trend for substage
table(pData(esets[["TCGA_eset"]])[,"batch"],pData(esets[["TCGA_eset"]])[,"substage"])
# b  c
# 9   0 30
# 11  0 32
# 12  3 30
# 13  4 33
# 14  0 35
# 15  5 17
# 17  4 38
# 18  1 40
# 19  6 37
# 21  3 33
# 22  1 42
# 24  2 34
# 27  0  6
# 40  2 40

outcomesAndCovariates <- pData(esets[["TCGA_eset"]])[,"batch",drop=FALSE]
rownames(outcomesAndCovariates) <- colnames(exprs(esets[["TCGA_eset"]]))
exprMatrix <- exprs(esets[["TCGA_eset"]])
rownames(exprMatrix) <- featureNames(esets[["TCGA_eset"]])

#p-value after batch correction:
batchCorrect <- Coincide::batchNormalization(countsMatrixNoNANoDup=exprMatrix,
                                   outcomesAndCovariates=outcomesAndCovariates,
                                   MinInBatch=4,
                                   combatModelFactorName=NULL,
                                   pvalueThresh=.05,
                                   batchColName="batch",
                                   outputFile="combatoutput.txt")
#original p-value:
# batchCorrect$beforePvalue
#[1] 6.129178e-05
#p-value after batch correction:
# batchCorrect$afterPvalue
#[1] 0.5019748
esets[["TCGA_eset"]] <- esets[["TCGA_eset"]][rownames(batchCorrect$GEN_Data_Corrected), colnames(batchCorrect$GEN_Data_Corrected)]
featureNames(esets[["TCGA_eset"]] ) <- rownames(batchCorrect$GEN_Data_Corrected)
validObject(esets[["TCGA_eset"]])


saveRDS(esets,file=paste0(saveDirGlobal,"esets_proc_TCGAcombat.rds"),compress=TRUE)

####
###select features

#now format just as a list of data matrices.
dataMatrixList <- exprSetListToMatrixList(esets,featureDataFieldName="gene")

names(dataMatrixList) <- names(esets)
###save this dataMatrixList now.
saveRDS(dataMatrixList,file=paste0(saveDirGlobal,"curatedOvarianData_procDataMatrixList.rds"))

###meta-features using top 20 by variance for each dataset:
#ran meta-feature analysis for 1000,500,1000,2000 features.
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_200.rds"),compress=TRUE)
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_500.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_1000.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=20,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_20,"/metaFeatures_2000.rds"),compress=TRUE)


###now do without top 20
#ran meta-feature analysis for 1000,500,1000,2000 features.
load(paste0(saveDirGlobal,"/curatedOvarianData_procDataMatrixList.rds"))
metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=200,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

#(saved each one after ran)
saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_200.rds"),compress=TRUE)

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=500,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_500.rds"),compress=TRUE)


metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=1000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_1000.rds"),compress=TRUE)


metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=2000,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_2000.rds"),compress=TRUE)

##now do 250, 300 - these were for my dissertation to answer the question
#if slightly tweaking the gene set size changed the clusters a lot.

metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=250,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=paste0(saveDir_no20,"/metaFeatures_250.rds"),compress=TRUE)


metaFeatures <- selectFeaturesMetaVariance_wrapper(dataMatrixList,rankMethod=c("mad"), 
                                                   numNAstudiesAllowedPerFeat=ceiling(length(dataMatrixList)/10),
                                                   numFeatSelectByGlobalRank=300,
                                                   numTopFeatFromEachDataset=0,fractNATopFeatAllowedPerDataset=.1,selectMethod=c("median"),
                                                   outputFile=outputFile)

message(paste0("Dropping the following datasets: ",paste(names(dataMatrixList)[metaFeatures$datasetListIndicesToRemove],collapse=" ")))

saveRDS(metaFeatures,file=pase0(saveDir_no20,"/metaFeatures_300.rds"),compress=TRUE)




