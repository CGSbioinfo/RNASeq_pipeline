suppressMessages(library(edgeR))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(GGally))


multipleComparison=function(data,comparison,design, gtf.file){
  # Initialize vector to record number of DE genes
  significant01=c()
  significant05=c()
  
  # Get annotation
  GTF <- import.gff(gtf.file, format="gtf", feature.type="gene")
  df=GTF$gene_name
  names(df)=GTF$gene_id
  
  for (i in 1:nrow(comparisons)){
    print(i)
    # Subset data and setup environment
    temp_comparison=c(as.character(comparisons[i,1]), as.character(comparisons[i,2]))
    temp_sample_info=sample_info[which(sample_info$Group %in% temp_comparison),]
    temp_sample_info=droplevels(temp_sample_info)
    temp_data=data[,which(colnames(data)%in%temp_sample_info$SampleID)]
    group=group <- c(as.character(temp_sample_info$Group))
    newd=paste0(temp_comparison, collapse = "__VS__")
    dir.create(paste0(outdir,'/',newd), showWarnings=FALSE)
    
    temp_colData<-data.frame(Group=temp_sample_info$Group)
    rownames(temp_colData)=temp_sample_info$SampleID

    dds <- DESeqDataSetFromMatrix(countData= temp_data, colData=temp_colData, design= ~Group)
    
    #  Prefiltering #
    #---------------#
    dds <- dds[rowSums(counts(dds)) > 1, ]

    # Normalizing dds #
    #-----------------#
    dds <- estimateSizeFactors(dds)
    pdf(paste0(outdir,'/',newd,'/size_factors_', newd, ".pdf"))
    plot(sizeFactors(dds), colSums(counts(dds)))
    abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
    dev.off()
    
    # Exploring data  #
    #-----------------#
    rld <- rlog(dds, blind=FALSE)

    # Heatmap genes
    #select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
    #heatmap.2(assay(rld)[select,], margins=c(12,12), scale=c('none'), density.info='density', trace='none', Rowv=FALSE, Colv=FALSE, dendogram='none')

    # Heatmap sample dist
    #y=assay(rld)
    #pdf(paste0(outdir,'/',newd,'/Heatmap_', newd, ".pdf"))
    #heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none', margins=c(12,12), cex.lab=0.8)
    #dev.off()

    # PCA
    #pdf(paste0(outdir,'/',newd,'/PCA_', newd, ".pdf"))
    #plotPCA(rld, intgroup='Group')
    #dev.off()
    
    # scatter plot
    png(paste0(outdir,'/',newd,'/ScatterPlot_', newd, ".png"), width = 1360, height = 1360)
    pairs.panels(assay(rld), smooth=FALSE)
    dev.off()
    
    # Design
    if (pairedDesign==TRUE){
      sibship=factor(temp_sample_info$Sibship)
      treatment=factor(temp_sample_info$Group)
      design(dds) <- ~sibship + treatment # this is not tested !!
    } else {
      design(dds) <- ~Group
    }
    
    # DE
    dds <- DESeq(dds)
    res <- results(dds)
    res <- res[order(res$padj),]
    #pdf(paste0(outdir, '/',newd,'/MAplot_', newd, ".pdf"))
    #plotMA(res)
    #dev.off()
    head(res)
    res <- as.data.frame(res)

    detags <- rownames(res)
    res = cbind(res,counts(dds,normalized=TRUE)[detags,]) # chech individual cpm values for 
    res=cbind(GeneName=df[rownames(res)], res)
    write.csv(res,paste0(outdir, '/',newd,'/Results_', newd, ".csv"))
    
    significant01=c(significant01, sum(res$padj <0.01, na.rm=TRUE))
    names(significant01)[i]=newd
    significant05=c(significant05, sum(res$padj <0.05, na.rm=TRUE))
    names(significant05)[i]=newd

    # Dispersion
    pdf(paste0(outdir, '/',newd,'/Dispersion_', newd, ".pdf"))
    plotDispEsts(dds)
    dev.off()
  }
  summ=cbind(significant01,significant05)
  colnames(summ)=c('p<0.01',"p<0.05")
  write.csv(summ,paste0(outdir,"/significant_summ.csv"),quote=F)  
}

########### scatterplots
panel.cor.scale <- function(x, y, digits=2, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex * abs(r))
}
panel.cor <- function(x, y, digits=4, prefix="")
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = (cor(x, y,use="pairwise"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 2)
}
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE)
{if (smooth ){
  if (scale) {
    pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth)
  }
  else {pairs(x,diag.panel=panel.hist,lower.panel=panel.smooth)
  } #else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth)
}
else #smooth is not true
{ if (scale) {pairs(x,diag.panel=panel.hist)
} else {pairs(x,diag.panel=panel.hist, cex.labels = 1.3) }
} #end of else (smooth)
} #end of function
