#!/usr/local/bin/Rscript

suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))

indir = commandArgs(TRUE)[1]
outdir = commandArgs(TRUE)[2]
gtf.file= commandArgs(TRUE)[3]

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

# Read GTF File
GTF <- import.gff(gtf.file, format="gtf", feature.type="gene")
df=as.data.frame(GTF)
df=subset(df,select=c(seqnames,start,end,width,strand,gene_id,gene_name,gene_biotype))
rownames(df)=df$gene_id
df=df[rownames(data),]

# Library proportion to different biotypes #
#------------------------------------------
type=unique(df$gene_biotype)
for (i in 1:length(type)){ 
	name=type[i]
	assign(name,colSums(data[which(df$gene_biotype==type[i]),])/colSums(data[,]))
}

matrix=matrix(ncol=ncol(data),nrow=length(type))
rownames(matrix)=type
colnames(matrix)=colnames(data)
for (i in 1:length(type)){
	matrix[type[i],]=eval(as.name(type[i]))
}
matrix=matrix[-which(rowSums(matrix)==0),]
matrix.melt=melt(matrix)
colnames(matrix.melt)[2]='Sample'
pdf(paste0(outdir,'/biotype_proportion.pdf'),width=12, height=7)
ggplot(matrix.melt, aes(x=X1, y=value, colour=Sample, shape=Sample)) + geom_point(size=1) + 
theme(axis.text.x=element_text(angle=90), axis.text=element_text(size=10, colour='black'), axis.title.x=element_blank(), legend.text=element_text(size=8)) + 
scale_shape_manual(values=0:ncol(data)) + ylab('Library Proportion')
dev.off()

# Library proportion #
#--------------------#
prop=t(t(data)/total.data.counted)
#boxplot(prop)

prop.adjusted=prop/df$width
#boxplot(prop.adjusted)

genes_most_expressed=sapply(1:ncol(prop.adjusted),function(x){which(prop.adjusted[,x]>.0000002)})
print(head(genes_most_expressed))
#names(genes_most_expressed)=colnames(prop.adjusted)
#genes_most_expressed=lapply(genes_most_expressed, function(x){df[names(x),]})
#genes_most_expressed=lapply(names(genes_most_expressed),function(x){cbind(genes_most_expressed[[x]],prop=prop[rownames(genes_most_expressed[[x]]),x])})
#names(genes_most_expressed)=colnames(prop.adjusted)

genes=unique(names(unlist(genes_most_expressed)))
prop.top=prop[genes,]
prop.top=cbind(prop.top,subset(df[genes,],select=c(seqnames,strand,gene_name,gene_biotype)))

#prop.top$variable=factor(prop.top$variable)
prop.top=melt(prop.top)
prop.top$variable=factor(prop.top$variable)
pdf(paste0(outdir,'/top_expressed_genes.pdf'),width=12, height=7)
ggplot(prop.top,aes(x=gene_name,y=value,group=variable,color=variable, shape=variable)) +
geom_point(size=1.5) + facet_wrap(~gene_biotype, scales='free_x')+ xlab('') + ylab('Library proportion') + theme(axis.text=element_text(size=6),axis.text.x=element_text(angle=90), axis.title=element_text(size=10), legend.text=element_text(size=8), legend.title=element_text(size=8)) + scale_shape_manual(values=0:ncol(data))
dev.off()
#legend.key.height=unit(.3, 'cm')

prop.top=prop.top[prop.top$gene_biotype=='protein_coding',]
pdf(paste0(outdir,'/top_expressed_genes_protein_coding.pdf'),width=12, height=7)
ggplot(prop.top,aes(x=gene_name,y=value,group=variable,color=variable, shape=variable)) +
geom_point(size=3) + facet_wrap(~gene_biotype, scales='free_x')+ xlab('') + ylab('Library proportion') + theme(axis.text=element_text(size=8, color='black'),axis.text.x=element_text(angle=90, vjust=0.5), axis.title=element_text(size=10), legend.text=element_text(size=8),  legend.title=element_text(size=8)) + scale_shape_manual(values=0:ncol(data))
dev.off()
#legend.key.height=unit(.3, 'cm'),
