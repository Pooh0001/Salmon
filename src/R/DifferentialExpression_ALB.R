#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=Time,y=value,col=Virus,shape=Cells,group=Cells)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter,lfcThreshold=lfc,alpha=padj)
    } else {
        res <- results(dds,contrast=contrast,filter = filter,lfcThreshold=lfc,alpha=padj)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
        message(sprintf(paste(
            ifelse(sum(sel)==1,
                   "There is %s gene that is DE",
                   "There are %s genes that are DE"),
                "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
                        sum(sel),padj,
                        lfc,expression_cutoff))
    }
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
        if (sum(val) >0){
            warning(sprintf(paste(
                ifelse(sum(val)==1,
                       "There is %s DE gene that has",
                       "There are %s DE genes that have"),
                       "no vst expression in the selected samples"),sum(val)))
            sel[sel][val] <- FALSE
        }
        
        if(export){
            if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
            write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot & sum(sel)>1){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis_macrogen/bioQA/dds.rda"))

#' ## Normalisation for visualisation
#' 
#' This time the normalisation takes advantage of the model, unlike during the biological QA
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis_macrogen/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis_macrogen/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis_macrogen/DE/vst-aware.tsv"))

#' ## Gene of interests
#goi <- read_lines(here("doc/goi.txt"))
# stopifnot(all(goi %in% rownames(vst)))
# dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
#' 
#' Because the design is unduly complex to answer the question asked,
#' we simplify it
#' 
#' ### U4.4
ddsU4 <- dds[,dds$Cells=="U4.4"]
design(ddsU4) <- ~Virus
ddsU4 <- DESeq(ddsU4)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(ddsU4)

#' Check the different contrasts
resultsNames(ddsU4)

#' #### Results
U4.4 <- extract_results(dds=ddsU4,vst=vst,
                        contrast="Virus_yes_vs_no",
                        default_dir=here("data/analysis_macrogen/DE"),
                        default_prefix="U4.4_Virus_vs_Ctl_",
                        sample_sel=dds$Cells=="U4.4",
                        labels=ddsU4$Virus)

#' ### C6.36
ddsC6 <- dds[,dds$Cells=="C6.36"]
design(ddsC6) <- ~Virus * Time
ddsC6$Time <- relevel(ddsC6$Time,"6h")
ddsC6 <- DESeq(ddsC6)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(ddsC6)

#' Check the different contrasts
resultsNames(ddsC6)

#' #### Results
#' ##### 6h
C6.6h <- extract_results(dds=ddsC6,vst=vst,
                        contrast="Virus_yes_vs_no",
                        default_dir=here("data/analysis_macrogen/DE"),
                        default_prefix="C6.6h_Virus_vs_Ctl_",
                        sample_sel=dds$Cells=="C6.36" & dds$Time=="6h",
                        labels=ddsC6$Virus)

dds$Time <- relevel(dds$Time,"6h")
line_plot(dds=dds,vst=vst,gene_id=C6.6h$all)

#' ##### 24h
C6.24h <- extract_results(dds=ddsC6,vst=vst,
                          contrast=c(0,1,0,1),
                          default_dir=here("data/analysis_macrogen/DE"),
                          default_prefix="C6.24h_Virus_vs_Ctl_",
                          sample_sel=dds$Cells=="C6.36",
                          labels=paste(ddsC6$Virus,dds$Time))

#' ##### 24h (interaction)
C6.24h_vs_6h <- extract_results(dds=ddsC6,vst=vst,
                         contrast="Virusyes.Time24h",
                         default_dir=here("data/analysis_macrogen/DE"),
                         default_prefix="C6.24h_vs_C6.6h_",
                         sample_sel=dds$Cells=="C6.36" & dds$Time=="24h",
                         labels=ddsC6$Virus)

line_plot(dds=dds,vst=vst,gene_id=C6.6h$all)

#' ##### Time effect
#' Just a sanity check
T24h_vs_6h <- extract_results(dds=ddsC6,vst=vst,
                                contrast="Time_24h_vs_6h",
                                export=FALSE,
                                sample_sel=dds$Cells=="C6.36",
                                labels=paste(ddsC6$Virus,dds$Time))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


