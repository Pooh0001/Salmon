#' ---
#' title: "Biological QA of LamV, WNV and dual-infected U4.4 cells"
#' author: "Nicolas Delhomme and Pontus Öhlund"
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
  library(ggtree)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(treeio)
  library(tximport)
  library(vsn)
})



#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Sample information
#' ```{r CHANGEME1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
samples <- read_csv(here("doc/CHANGEME.csv"),
                    col_types=cols(.default=col_factor()))


#' # Raw data
filelist <- list.files(here("CHANGEME/"), 
                       recursive = TRUE, 
                       pattern = "quant.sf",
                       full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info
#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are in the same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
stopifnot(all(match(basename(dirname(filelist)),
                    samples$SampleID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$SampleID

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, replace with the following
#' txi <- suppressMessages(tximport(files = filelist, type = "salmon",txOut=TRUE))
#' ```
txi <- suppressMessages(tximport(files = filelist,
                                 type = "salmon",
                                 txOut=TRUE))
counts <- txi$counts

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by CHANGEME
#' ```{r CHANGEME5,eval=FALSE,echo=FALSE}
#' # In the following most often you need to replace CHANGEME by your
#' # variable of interest, i.e. the metadata represented as column in
#' # your samples object, e.g. SamplingTime
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Treatment)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

ggplot(dat,aes(x,y,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())


#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by CHANGEME. 
#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Treatment=samples$Treatment[match(ind,samples$SampleID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=Treatment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("CHANGEME/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("CHANGEME/analysis/salmon/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ Treatment)

save(dds,file=here("CHANGEME/analysis/salmon/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y",ylim=c(0.75,1.25))


#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=2

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables aas needed.
#' ```
nlevel=nlevels(dds$Treatment) * nlevels(dds$Time)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)

#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Treatment,shape=Time,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Sequencing depth
#' Number of genes expressed per condition at different cutoffs
conds <- factor(paste(dds$Treatment,dds$Time))
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' ### Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 1

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,
                col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",labels=dds$Treatment)

#' ### Hierarchical clustering
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 1000, parallel = TRUE)

#' plot the clustering with bp and au
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' fancy
tree2 <- full_join(as.treedata(hm.pvclust),
                   as.data.frame(colData(dds)) %>% mutate(label=conds), by='label')

ggtree(tree2) + 
  geom_tippoint(aes(shape=Time)) +
  geom_tiplab(size=3,aes(label=Treatment),hjust=-0.5) +
  geom_nodelab(aes(label=au,color=au),hjust=-0.5,size=3) + 
  scale_color_continuous(low="red", high="blue")


#' bootstrapping results as a table
print(hm.pvclust, digits=3)


#' ## Conclusion
#' The data quality looks good. The sequencing depth is homogeneous. The PCA shows the 
#' expected clustering of the data. There is enough variance in the data to expect
#' differential expression.
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
