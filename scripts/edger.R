#!/usr/local/bin/Rscript

suppressMessages(library(edgeR))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(GGally))
source("/usr/local/bin/edger_functions.R")

indir = commandArgs(TRUE)[1]
outdir = commandArgs(TRUE)[2]
sample_info = commandArgs(TRUE)[3]
comparisons = commandArgs(TRUE)[4]
min.count = as.numeric(commandArgs(TRUE)[5]) # filtering: minimun number of reads a gene should have to be considered as expressed
min.nsamples = as.numeric(commandArgs(TRUE)[6]) # filtering: minimum number of samples having a gene expressed
design = commandArgs(TRUE)[7]
gtf.file= commandArgs(TRUE)[8]

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


# Creating a DGE object to analyse data with edgeR #
#--------------------------------------------------#
dge <- DGEList(counts=data, group=group)


# Looking at mds before filtering
dge.temp <- calcNormFactors(dge)
mycoldf=cbind(unique(as.character(dge.temp$samples$group)),1:length(unique(as.character(dge.temp$samples$group))))
mycol=c()
for (i in 1:length(as.character(dge.temp$samples$group))){
  mycol=c(mycol,mycoldf[which(mycoldf[,1]==as.character(dge.temp$samples$group)[i]),2])
}
pdf(paste0(outdir,'/mds_normalised_noFiltering', ".pdf"))
plotMDS(dge.temp, col=mycol, method="bcv")
dev.off()


# Calculating the filtering threshold #
#-------------------------------------#
smallest_lib=min(dge$samples$lib.size)
smallest_lib_pm=smallest_lib/1000000
min.cpm=min.count/smallest_lib_pm


# Filtering dge #
#---------------#
keep <- rowSums(cpm(dge)>min.cpm) >= min.nsamples
dge.all.samples <- dge[keep,]
dge.all.samples$samples$lib.size <- colSums(dge.all.samples$counts)


# Normalizing dge #
#-----------------#
dge.all.samples <- calcNormFactors(dge.all.samples)


# Exploring data  #
#-----------------#
mycoldf=cbind(unique(as.character(dge.all.samples$samples$group)),1:length(unique(as.character(dge.all.samples$samples$group))))
mycol=c()
for (i in 1:length(as.character(dge.all.samples$samples$group))){
  mycol=c(mycol,mycoldf[which(mycoldf[,1]==as.character(dge.all.samples$samples$group)[i]),2])
}
pdf(paste0(outdir,'/MDSPlot_allSamples_bcv', ".pdf"))
plotMDS(dge.all.samples, col=mycol, method="bcv")# , xlim=c(-1.4,1.4),ylim=c(-.7,.7))# 
dev.off()
pdf(paste0(outdir,'/MDSPlot_allSamples_logFc', ".pdf"))
plotMDS(dge.all.samples, col=mycol, method="logFC")# , xlim=c(-1.4,1.4),ylim=c(-.7,.7))# 
dev.off()

y = cpm(dge.all.samples,prior.count = 1, log=TRUE)
pdf(paste0(outdir,'/Heatmap_allSamples', ".pdf"), width=9,height=9)
heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none', cex.lab=0.8)
dev.off()

pdf(paste0(outdir,'/Libsize', ".pdf"), width=9)
par(mar=c(14.8,4.1,4.1,2.1))
barplot(dge.all.samples$samples$lib.size*1e-6, names=colnames(data), las=2, ylab="Library size (millions)") # library sizes
dev.off()


# Doing Differential Gene Expression  #
#-------------------------------------#
comparisons=read.csv(comparisons, header=TRUE)
if (design=='pairedSamples'){
  pairedDesign=TRUE
} else if(design=='non-pairedSamples'){
  pairedDesign=FALSE
}

multipleComparison(dge,comparisons,pairedDesign, min.count, min.nsamples, gtf.file)
