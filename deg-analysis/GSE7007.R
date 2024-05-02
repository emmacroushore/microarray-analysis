# sample script used to fetch raw publicly-available microarray data for reanalysis
# BioProject:	 	PRJNA99237
# GEO Accession: GSE7007 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7007)
# Citation: Tirode F, Laud-Duval K, Prieur A, Delorme B et al. Mesenchymal stem cell features of Ewing tumors. 
#           Cancer Cell 2007 May;11(5):421-9. PMID: 17482132

# version info: R 4.1.2, Biobase 2.54.0, GEOquery 2.62.2, limma 3.50.3
################################################################
#  data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)
library(tidyverse)

```{r read_from_GEO, message = FALSE, warnings = FALSE, error = FALSE}
# load series and platform data from GEO
gset <- getGEO("GSE7007", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "X44552233001XXXXXXXXXXXXXXXXXXXXXXXXXXX"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

## getting just sample, expression, probe data from GSE
ex <- exprs(gset)

# log2 transformation
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

## extract sample annotations
pdata <- pData(gset)

## subset relevant conditions (cell line, treatment) from pdata
cond <- pdata %>%
  dplyr::select(title, geo_accession) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  dplyr::mutate(title = stringr::str_replace_all(string = title,
                      pattern = "_rep1", replacement = "")) %>%
  dplyr::mutate(title = stringr::str_replace_all(string = title,
                      pattern = "_rep2", replacement = "")) %>%
  as.data.frame()

## make design matrix
cond$title <- factor(cond$title, levels=c("EW24_shCTL", "EW24_shEF1",
                                          "SKNMC_shCTL","SKNMC_shEF1", 
                                          "A673_siCT", "A673_siEF1"))
# set rownames
rownames(cond) <- cond$geo_accession

# create model
model <- model.matrix(~0 + title, cond)

# fit linear model to gene expression data
fit <- lmFit(gset, model)

# store annotation data on genes
annot <- fit$genes

# subset probe id, gene title, and gene symbol
annot_tidy <- annot %>%
  dplyr::select(1:3)

# make the contrasts between conditions of interest
cont.matrix <- makeContrasts(EW24 = titleEW24_shEF1 - titleEW24_shCTL, 
                             SKNMC = titleSKNMC_shEF1 - titleSKNMC_shCTL, 
                             A673 = titleA673_siEF1 - titleA673_siCT, 
                             levels=model)


# apply contrasts to fitted model
fit2 <- contrasts.fit(fit, cont.matrix)

# eBayes statistical testing
fit2 <- eBayes(fit2)

# results summary with topTable
topTable(fit2)

# subsetting individual cell lines and their differentially regulated genes between conditions (p < 0.05)
# note: per limma documentation, the "logFC" column produced by topTable represents log2FC
#       see: https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf (p. 61)

# A673 cell line
all.A673 <- topTable(fit2, coef = "A673", number = Inf, adjust.method = "fdr")

all.A673.tidy <- all.A673 %>%
  dplyr::select(1:3,22,26) %>%
  group_by(Gene.symbol) %>%
  dplyr::filter(adj.P.Val == min(adj.P.Val))

colnames(all.A673.tidy) <- c("probe","description","symbol","log2FC","adj.p.val")

deg.A673 <- dplyr::filter(all.A673.tidy, adj.p.val <= 0.05)
upreg.A673 <- dplyr::filter(deg.A673, log2FC > 1)

# SK-N-MC cell line
all.SKNMC <- topTable(fit2, coef = "SKNMC",number = Inf, adjust.method = "fdr")

all.SKNMC.tidy <- all.SKNMC %>%
  dplyr::select(1:3,22,26)%>%
  group_by(Gene.symbol) %>%
  dplyr::filter(adj.P.Val == min(adj.P.Val))

colnames(all.SKNMC.tidy) <- c("probe","description","symbol","log2FC","adj.p.val")

deg.SKNMC <- dplyr::filter(all.SKNMC.tidy, adj.p.val < 0.05)
upreg.SKNMC <- dplyr::filter(deg.SKNMC, log2FC > 1)

# EW-24 cell line
all.EW24 <- topTable(fit2, coef = "EW24",number = Inf, adjust.method = "fdr")

all.EW24.tidy <- all.EW24 %>%
  dplyr::select(1:3,22,26) %>%
  group_by(Gene.symbol) %>%
  dplyr::filter(adj.P.Val == min(adj.P.Val))

colnames(all.EW24.tidy) <- c("probe","description","symbol","log2FC","adj.p.val")

# subset all differentially regulated genes (adj.p.val < 0.05)
deg.EW24 <- dplyr::filter(all.EW24.tidy, adj.p.val < 0.05)

# subset differentially upregulated genes (adj.p.val < 0.05, log2FC > 1)
upreg.EW24 <- dplyr::filter(deg.EW24, log2FC > 1)



