#' ---
#' title: "Macrogen data Biological QA"
#' author: "Nicolas Delhomme & Juliette Hayer"
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
  library(parallel)
  library(pander)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Metadata
#' Sample information
samples <- read_csv(here("doc/macrogen-samples.csv"),
                      col_types=cols(.default=col_factor()))

                              
#' JH: I use the Expression Profile.gene.txt from StringTie and I select only 
#' the Feature_GID and the read_count columns of each sample.
#' 
#' I transform the column Feature_GID into the row names => this will be my 
#' 'counts' object
#' 
#' ND I added a filter to remove small RNA ("misc_RNA", "rRNA", "snRNA", "snoRNA","tRNA"  )
counts <- suppressWarnings(suppressMessages(read_tsv(here("data/macrogen/ALB_RNAseq_for_Nico/R_analysis_JH/salmon/Expression_Profile.GCF_006496715.gene.txt"), 
                   col_names = TRUE, col_types = cols(.default=col_guess())) %>% 
  filter(Type %in% c("protein_coding","pseudogene","lncRNA","transcribed_pseudogene")) %>% 
  select(c("Feature_GID",ends_with("Read_Count")))  %>% 
  column_to_rownames("Feature_GID")))

#' Sanity check to ensure that the data is sorted according to the sample info
#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#' # This step is to validate that the salmon files are in the same order as 
#' # described in the samples object. If not, then they need to be sorted
#' ````
stopifnot(all(colnames(counts) == samples$Sample_ID))

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, colouring by CHANGEME
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Virus)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

ggplot(dat,aes(x,y,fill=Time)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

ggplot(dat,aes(x,y,fill=Cells)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage shows a very high shoulder on the right side.
#' 
#' A lot of genes are lowly expressed. This is probably indicating too high a sequencing depth.
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by the different meta-data. 
#' 
#' There is no obvious deviation between samples
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Virus=samples$Virus[match(ind,samples$Sample_ID)]) %>% 
  mutate(Time=samples$Time[match(ind,samples$Sample_ID)]) %>%
  mutate(Cells=samples$Cells[match(ind,samples$Sample_ID)])

ggplot(dat,aes(x=values,group=ind,col=Virus)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

ggplot(dat,aes(x=values,group=ind,col=Time)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

ggplot(dat,aes(x=values,group=ind,col=Cells)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis_macrogen/bioQA"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis_macrogen/bioQA/raw-unormalised-gene-expression_data.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
## JH: I do not have txi so I use DESeqDataSetFromMatrix instead
## JH: I block the effect of the Cells in the design, as we do not want yet to model it in the DEG 1 to 4

dds <- DESeqDataSetFromMatrix(
  counts,
  colData = samples,
  design = ~ Cells + Virus * Time
)

save(dds,file=here("data/analysis_macrogen/bioQA/dds.rda"))

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
#' The sequencing depth is relatively uniform -15% / + 30% of the average
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Assess whether there might be a difference in library size linked to a
#' given metadata
#' 
#' We should keep an eye on the Virus=yes samples as they seem to have 
#' a slightly deeper sequencing depth than the no counterpart.
boxplot(split(sizes,dds$Virus),las=2,
        main="Sequencing libraries size factor by treatment with Virus")

boxplot(split(sizes,dds$Time),las=2,
        main="Sequencing libraries size factor by Time")

boxplot(split(sizes,dds$Cells),las=2,
        main="Sequencing libraries size factor by Cell liness")

#' It does not seem all so bad, a couple samples are outliers for the "yes" class.
plot(sizes,log10(colSums(counts(dds))),ylab="log10 raw depth",xlab="scaling factor",
     col=rainbow(n=nlevels(dds$Virus))[as.integer(dds$Virus)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$Virus)),
       legend=levels(dds$Virus),cex=0.6)

#' ## Variance Stabilising Transformation
#' 
#' For visualisation, in addition to the correction for the sequencing depth above,
#' we want to correct for the mean-variance relationship that exists in RNA-Seq data,
#' hence we apply a transformation to render the data homoscedastic
#' 
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately, impressively so even.
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=3

#' An the number of possible combinations
nlevel=nlevels(dds$Virus) * nlevels(dds$Time) * nlevels(dds$Cells)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
#' 
#' A lot of the variance is explained by the first components, while all samples
#' are needed to explain all the variance in the data
#' 
#' The inflection point arise early (4) after which the increase is linearly constant
#' 
#' All these are promising.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  
#' ### 2D
#' The first components separates the cell types. The second is the time within C6.36
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Cells,text=Sample_ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Sequencing depth
#' Number of genes expressed per condition at different cutoffs
#' 
#' Here we try to assess whether there would be any link between sequencing depth and the 
#' variable of interest. 
#' 
#' It would seem not to be the case or limited to highly expressed genes in the case of
#' a Virus yes _vs._ no comparison
conds <- factor(paste(dds$Cells,dds$Time))
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=dds$Virus,
                                nrep=6)
#' ### Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 1

#' * Heatmap of "all" genes
#' The heatmap is showing the same as the PCA, Time and Cells are preponderant
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = paste(conds,dds$Virus),
          col=hpal)

#' There are a few interesting things here. U4.4 further separates based on virus presence.
#' 
#' C6.36 does not at 24h and neither does it at 6h unless there would be a sample swap...
plot(as.hclust(hm$colDendrogram),xlab="",sub="",labels=paste(conds,dds$Virus))

#' ### Hierarchical clustering
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 1000, parallel = TRUE)

#' There is little doubt about the clustering, the au score is 
#' always close to the maximum (100 for a 0-100 range)
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' ## Conclusion
#' The data is of good quality. It might have been sequenced too deeply.
#' 
#' The data might not satisfy the expected study design, as the Virus 
#' effect is minimal and only appears in the U4.4 cell type.
#' 
#' A putative sample swap for the C6.36 cells at 6h could be investigated, but it 
#' might just be a random effect due to the little variance effect created by the
#' Virus variable overall.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
