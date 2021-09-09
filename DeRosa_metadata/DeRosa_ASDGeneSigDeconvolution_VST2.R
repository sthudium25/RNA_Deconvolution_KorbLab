library(jaffelab)
library(GEOquery)
library(biomaRt)
library(RColorBrewer)
library(tidyverse)
library(stringr)

# Load GEO data from desired experiment
analysisName = "DeRosa_Deconvolution_log2FPKM_ASDSig"
samplesUsed = read_tsv("DeRosa_metadata/DeRosa_BioSamples", col_names = F)
GEO_data = getGEO("GSE124308")[[1]]

# If downloading count data from the internet, place URL here and save file as a .rda
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102741&format=file&file=GSE102741%5Flog2RPKMcounts%2Erda%2Egz",
#              destfile = "rdas/Wright2017_hASD_GSE102741.rda")
# This results in the object named 'countData'
#load("rdas/Wright2017_hASD_GSE102741.rda")

# Read in count data
countData = read.csv("GeneCountFiles/DeRosa2017_Deconvolution_FPKM.csv", row.names = 1) 

# If using Salmon TPMs, remove the transcript ID
rownames(countData) = gsub("\\..*", "", rownames(countData)) 

## Rename and transform expression data
GEO_data_pd= pData(GEO_data)
str(GEO_data_pd)
colnames(GEO_data_pd) = ss(colnames(GEO_data_pd), ":")
colnames(GEO_data_pd)= gsub(" ", "_", colnames(GEO_data_pd))

## Rename 'group' column to Dx
GEO_data_pd = GEO_data_pd %>% 
  rename(Dx = disease_state) 
GEO_data_pd$Dx = as_factor(GEO_data_pd$Dx) 
GEO_data_pd$Dx = fct_rev(GEO_data_pd$Dx)
levels(GEO_data_pd$Dx) = c("Control", "ASD")
GEO_data_pd$iPSc_source = as_factor(GEO_data_pd$iPSc_source)
GEO_data_pd$reprogramming_method = as_factor(GEO_data_pd$reprogramming_method)

GEO_data_pd_filtered = GEO_data_pd %>% 
  filter(geo_accession %in% samplesUsed$X1)

## Create signature gene set from mouse ensembl gene list
korb_signature = read_tsv("genelists/korbAllis_ASDGeneSig_Down_ensembl.txt", col_names = F)
referenceSignature = korb_signature %>% deframe()

ensembl=useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
MMtoHG = getBM(attributes = c('ensembl_gene_id',
                              'hsapiens_homolog_ensembl_gene', 
                              'hsapiens_homolog_associated_gene_name'),
               mart = ensembl)

signatureIndex = rownames(countData) %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% referenceSignature]

## Generate eigengene (original)
pca_analysis = prcomp(t(countData[signatureIndex,]))
eigenGene = prcomp(t(countData[signatureIndex,]))$x[,1]
pcaVars = getPcaVars(prcomp(t(log2(countData[signatureIndex,]+1))))[1]
pcaVars_alt = getPcaVars(prcomp(t(countData[signatureIndex,])))[1]
loads = pca_analysis$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(names(loads),
                                                                  MMtoHG$hsapiens_homolog_ensembl_gene)] 

## Log2 + 1 scaled FPKM
pca_analysis = prcomp(t(log2CountData[signatureIndex,]))
eigenGene = prcomp(t(log2CountData[signatureIndex,]))$x[,1]
pcaVars = getPcaVars(prcomp(t(log2(log2CountData[signatureIndex,]+1))))[1]
loads = pca_analysis$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(names(loads),
                                                                  MMtoHG$hsapiens_homolog_ensembl_gene)]

principalComps = prcomp(t(log2CountData))
f = summary(lm(eigenGene ~ Dx + principalComps$x[,1:5], data=GEO_data_pd_filtered))
f

pdf(sprintf("figures/%s_korbASDSigDown_deconv_uncleaned1.pdf", analysisName), h=6, w=5)
par(mar=c(5, 6, 2, 2), cex.axis=2, cex.lab=2)
palette(brewer.pal(5, "Set1"))
boxplot(eigenGene ~ GEO_data_pd_filtered$Dx, outline=FALSE,
        ylim = range(eigenGene), xlab="",
        ylab=paste0("ASDSig Eigengene"))
points(eigenGene ~ jitter(as.numeric(GEO_data_pd_filtered$Dx), amount=0.15), 
        pch = 21, bg=factor(GEO_data_pd_filtered$Dx))
legend("topright", paste0("p=", signif(f$coef[2,4],3)),cex=1.5)
title(main = sprintf("%s\nacct Disease + 5 PCs / uncleaned", analysisName))
dev.off()

# clean up the eigenGene data


fPC = lm(eigenGene ~ Dx + principalComps$x[,1:3], 
         data=GEO_data_pd_filtered)
summary(fPC)
mod = model.matrix(~ Dx + principalComps$x[,1:3], 
                   data=GEO_data_pd_filtered)
ff = summary(lm(eigenGene ~ mod - 1))
ff
cleanEigen = cleaningY(matrix(eigenGene, nr=1), mod, P=2)

# Regenerate plot from cleaned data
pdf(sprintf("figures/%s_korbASDSigDown_deconv1.pdf", analysisName),h=6,w=5)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen[1,] ~ GEO_data_pd_filtered$Dx, outline = FALSE,
        ylim = range(cleanEigen[1,]),xlab="",
        ylab = paste0("ASD Gene Sig Eigengene (Adj)"))
points(cleanEigen[1,] ~ jitter(as.numeric(factor(GEO_data_pd_filtered$Dx)), amount=0.15),
       pch = 21, bg=factor(GEO_data_pd_filtered$Dx))
legend("topright", paste0("p=", signif(ff$coef[2,4],3)),cex=1.5)
title(main = sprintf("%s\nacct for Dx + 3 PCs \ cleaned", analysisName))
dev.off()

# Generate correlation matrix

GEO_data_pd_filtered$ados_total_score_at_biopsy = as.numeric(GEO_data_pd_filtered$ados_total_score_at_biopsy)
GEO_data_pd_filtered$age_at_biopsy = as.numeric(GEO_data_pd_filtered$age_at_biopsy)
GEO_subset = GEO_data_pd_filtered %>% select(age_at_biopsy, Dx, ados_total_score_at_biopsy)
GEO_subset = GEO_subset %>% mutate(ASD = if_else(Dx == "ASD", 1, 0)) %>% select(-Dx)
temp = as.data.frame(mod[,3:12])
GEO_subset = GEO_subset %>% bind_cols(temp)
cor_mat = cor(GEO_subset, use="complete.obs")

GEO_subset_noDx = GEO_subset %>% select(-ASD)
cor_mat_noDx = cor(GEO_subset_noDx,use="complete.obs")

pheatmap(cor_mat_noDx)



