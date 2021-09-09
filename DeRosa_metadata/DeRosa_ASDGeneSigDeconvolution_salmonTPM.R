library(jaffelab)
library(GEOquery)
library(biomaRt)
library(RColorBrewer)
library(tidyverse)
library(stringr)

# Load GEO data from desired experiment
analysisName = "DeRosaD35_Deconvolution_SalmonTPMs"
samplesUsed = read_tsv("DeRosa_metadata/DeRosa_BioSamples", col_names = F)
GEO_data = getGEO("GSE124308")[[1]]

# If downloading count data from the internet, place URL here and save file as a .rda
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102741&format=file&file=GSE102741%5Flog2RPKMcounts%2Erda%2Egz",
#              destfile = "rdas/Wright2017_hASD_GSE102741.rda")
# This results in the object named 'countData'
#load("rdas/Wright2017_hASD_GSE102741.rda")

# Read in count data
countData = read.csv("GeneCountFiles/DeRosa2017_Deconvolution_SalmonTPM.csv", row.names = 1) 
rownames(countData) = gsub("\\..*", "", rownames(countData)) 

## Rename and transform expression data
GEO_data_pd= pData(GEO_data)
str(GEO_data_pd)
colnames(GEO_data_pd) = ss(colnames(GEO_data_pd), ":")
colnames(GEO_data_pd)= gsub(" ", "_", colnames(GEO_data_pd))

## Add a new factor column for disease status
Dx = GEO_data_pd %>%
  select(contains("disease")) %>% 
  mutate(Dx = factor(ifelse(. == 'healthy control', 'CONT', 'ASD'), 
                     levels = c("CONT", "ASD"))) %>%
  select(Dx)
GEO_data_pd = GEO_data_pd %>% mutate(Dx)

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

## Generate eigengene
eigenGene = prcomp(t(countData[signatureIndex,]))$x[,1]
pcaVars = getPcaVars(prcomp(t(log2(countData[signatureIndex,]+1))))[1]
pcaVars_alt = getPcaVars(prcomp(t(countData[signatureIndex,])))[1]
loads = prcomp(t(countData[signatureIndex,]))$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(names(loads),
                                                                  MMtoHG$hsapiens_homolog_ensembl_gene)] 
boxplot(eigenGene ~ GEO_data_pd_filtered$Dx)

# clean up the eigenGene data
f = lm(eigenGene ~ Dx, data=GEO_data_pd_filtered)

summary(f)

principalComps = prcomp(t(countData))
fPC = lm(eigenGene ~ Dx + principalComps$x[,1:7], data=GEO_data_pd_filtered)

mod = model.matrix(~Dx + principalComps$x[,1:7], data=GEO_data_pd_filtered)
ff = summary(lm(eigenGene ~ mod -1))
ff
cleanEigen = cleaningY(matrix(eigenGene, nr=1), mod, P=2)

# Regenerate plot from cleaned data
pdf(sprintf("figures/%s_korbASDSigDown_deconvolution.pdf", analysisName),h=6,w=5)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen[1,] ~ GEO_data_pd_filtered$Dx, outline = FALSE,
        ylim = range(cleanEigen[1,]),xlab="",
        ylab = paste0("CAG Eigengene (Adj)"))
points(cleanEigen[1,] ~ jitter(as.numeric(factor(GEO_data_pd_filtered$Dx)), amount=0.15),
       pch = 21, bg=factor(GEO_data_pd_filtered$Dx))
legend("topright", paste0("p=", signif(ff$coef[2,4],3)),cex=1.5)
dev.off()

