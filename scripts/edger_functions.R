suppressMessages(library(edgeR))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(GGally))


multipleComparison=function(dge,comparison,design, min.count, min.nsamples, gtf.file){
  # Initialize vector to record number of DE genes
  significant01=c()
  significant05=c()
  
  # Get annotation
  GTF <- import.gff(gtf.file, format="gtf", asRangedData=F, feature.type="gene")
  df=GTF$gene_name
  names(df)=GTF$gene_id
  
  for (i in 1:nrow(comparisons)){
    
    # Subset data and setup environment
    temp_comparison=c(as.character(comparisons[i,1]), as.character(comparisons[i,2]))
    temp_sample_info=sample_info[which(sample_info$Group %in% temp_comparison),]
    temp_sample_info=droplevels(temp_sample_info)
    temp_data=data[,which(colnames(data)%in%temp_sample_info$SampleID)]
    group=group <- c(as.character(temp_sample_info$Group))
    newd=paste0(temp_comparison, collapse = "__VS__")
    dir.create(paste0(outdir,'/',newd), showWarnings=FALSE)
    
    dge <- DGEList(counts=temp_data, group=group)
    
    # Calculating the filtering threshold #
    #-------------------------------------#
    smallest_lib=min(dge$samples$lib.size)
    smallest_lib_pm=smalles_lib/1000000
    min.cpm=min.count/smallest_lib_pm
    
    # Filtering dge #
    #---------------#
    keep <- rowSums(cpm(dge)>min.cpm) >= min.nsamples
    dge <- dge[keep,]
    dge$samples$lib.size <- colSums(dge$counts)
    
    # Normalizing dge #
    #-----------------#
    dge <- calcNormFactors(dge)
    
    # Exploring data  #
    #-----------------#
    mycoldf=cbind(unique(as.character(dge$samples$group)),1:length(unique(as.character(dge$samples$group))))
    mycol=c()
    for (j in 1:length(as.character(dge$samples$group))){
      mycol=c(mycol,mycoldf[which(mycoldf[,1]==as.character(dge$samples$group)[j]),2])
    }
    pdf(paste0(outdir,'/',newd,'/MDS_', newd, ".pdf"))
    plotMDS(dge, col=mycol, method="bcv")
    dev.off()
    
    # Heatmap
    y = cpm(dge,prior.count = 1, log=TRUE)
    pdf(paste0(outdir,'/',newd,'/Heatmap_', newd, ".pdf"))
    heatmap.2(cor(y),scale=c('none'), density.info='density', trace='none')
    dev.off()
    
    # scatter plot
    y = cpm(dge,prior.count = 1, log=TRUE)
    png(paste0('differentialExpression/',newd,'/ScatterPlot_', newd, ".png"), width = 1360, height = 1360)
    pairs.panels(y, smooth=FALSE)
    dev.off()
    
    # Dispersion
    if (pairedDesign==TRUE){
      sibship=factor(temp_sample_info$Sibship)
      treatment=factor(temp_sample_info$Group)
      design<-model.matrix(~sibship+treatment)
      dge<-estimateGLMCommonDisp(dge, design)
      dge<-estimateGLMTrendedDisp(dge, design)
      dge<-estimateGLMTagwiseDisp(dge, design)
    } else {
      dge <- estimateCommonDisp(dge) # Common dispersio estimates the overall BCV of the dataset, averaged over all genes # BCV is the sqroot of the common dispersion
      dge <- estimateTagwiseDisp(dge) # Estimaes gene specific dispersions
    }
    
    pdf(paste0(outdir,'/',newd,'/Dispersion_', newd, ".pdf"))
    plotBCV(dge)
    dev.off()
    
    # Differential Expression
    if (pairedDesign==TRUE){
      fit<-glmFit(dge,design)
      lrt<-glmLRT(fit, coef=7)
      results=topTags(lrt, n=dim(lrt)[1])[[1]]
    } else {
      et <- exactTest(dge, pair=c(temp_comparison[2], temp_comparison[1]))
      results=topTags(et, n=dim(et)[1])[[1]]
    }
    detags <- rownames(results) # chech individual cpm values for top genes
    results = cbind(results,cpm(dge,log=TRUE)[detags,]) # chech individual cpm values for top genes
    
    results=cbind(GeneName=df[rownames(results)], results)
    write.csv(results,paste0(outdir, '/',newd,'/Results_', newd, ".csv"))
    
    significant01=c(significant01, table(results$FDR<0.01)["TRUE"])
    names(significant01)[i]=newd
    significant05=c(significant05, table(results$FDR<0.05)["TRUE"])
    names(significant05)[i]=newd
    
    if (pairedDesign==TRUE){
      de <- decideTestsDGE(lrt, p=0.05, adjust="BH") # total number of DE genes at 5% FDR
      detags <- rownames(dge)[as.logical(de)]
      pdf(paste0(outdir, '/',newd,'/MAplot_', newd, ".pdf"))
      plotSmear(lrt, de.tags=detags)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
    } else {
      de <- decideTestsDGE(et, p=0.05, adjust="BH") # total number of DE genes at 5% FDR
      detags <- rownames(dge)[as.logical(de)]
      pdf(paste0(outdir, '/',newd,'/MAplot_', newd, ".pdf"))
      plotSmear(et, de.tags=detags)
      abline(h = c(-1, 1), col = "blue")
      dev.off()
    }
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
