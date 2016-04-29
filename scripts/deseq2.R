#!/usr/local/bin/Rscript

suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(GGally))
source("/usr/local/bin/deseq2_functions.R")

indir = commandArgs(TRUE)[1]
outdir = commandArgs(TRUE)[2]
sample_info = commandArgs(TRUE)[3]
comparisons = commandArgs(TRUE)[4]
design = commandArgs(TRUE)[5]
gtf.file= commandArgs(TRUE)[6]

dir.create(outdir, showWarnings = FALSE)

# Reading counted reads files #
#-----------------------------#

files=list.files(indir, pattern = "_count.txt", full.names=TRUE)

data=read.table(files[1], row.names=1)
colnames(data)=files[1]

for (i in 2:length(files)){
  name=files[i]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1])
  data=cbind(data, sample[index,2])
  colnames(data)[i]=name
}
colnames(data)=gsub('_count.txt','',colnames(data))
colnames(data)=gsub('.*/','',colnames(data))
total_counts=colSums(data)

# Excluding reads not falling in features, ambiguous, too low qual, multiple locations
data.notcounted=data[which(rownames(data)=='__no_feature'):dim(data)[1],]

# Subsetting including reads
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),]
total.data.counted=colSums(data)


# Getting information about the groups #
#--------------------------------------#
sample_info = read.csv(sample_info)
sample_info=sample_info[order(sample_info$Group),]
group <- c(as.character(sample_info$Group))

# Matching the data matrix to the sample_info table #
#---------------------------------------------------#
data=data[,match(as.character(sample_info$SampleID),colnames(data))]

# Creating a DDS object to analyse data with DESEQ2 #
#--------------------------------------------------#

#colData<-data.frame(Group=sample_info$Group)
#rownames(colData)=sample_info$SampleID

#dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design= ~Group)

# Prefiltering
#dds <- dds[rowSums(counts(dds)) > 1, ]
#dds <- DESeq(dds)

# Exploring data  #
#-----------------#
# Transforming values
#rld <- rlog(dds, blind=FALSE)

# Heatmap genes
#select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
#heatmap.2(assay(rld)[select,], margins=c(8,12), scale=c('none'), density.info='density', trace='none', Rowv=FALSE, Colv=FALSE, dendogram='none')

# Heatmap sample dist
#y=assay(rld)
#pdf(paste0(outdir,'/Heatmap_allSamples', ".pdf"), width=9,height=9)
#heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none', margins=c(8,8), cex.lab=0.8)
#dev.off()

# PCA
#plotPCA(rld)


# Doing Differential Gene Expression  #
#-------------------------------------#
comparisons=read.csv(comparisons, header=TRUE)
if (design=='pairedSamples'){
  pairedDesign=TRUE
} else if(design=='non-pairedSamples'){
  pairedDesign=FALSE
}

multipleComparison(data,comparisons,pairedDesign, gtf.file)
