library(jaffelab)
library(GEOquery)
library(biomaRt)
library(RColorBrewer)
library(tidyverse)

# Load GEO data from desired experiment
wrightASD = getGEO("GSE102741")[[1]]
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102741&format=file&file=GSE102741%5Flog2RPKMcounts%2Erda%2Egz",
              destfile = "rdas/Wright2017_hASD_GSE102741.rda")
# This results in the object named 'geneRpkm2'
load("rdas/Wright2017_hASD_GSE102741.rda")

## Rename and transform expression data
wrightASD_pd = pData(wrightASD)

colnames(wrightASD_pd)[54:64] = ss(colnames(wrightASD_pd)[54:64], ":")
colnames(wrightASD_pd)[54:64] = gsub(" ", "_", colnames(wrightASD_pd)[54:64] )

## Add a new factor column for disease status
wrightASD_pd$Dx = factor(ifelse(wrightASD_pd$disease_status == "Healthy control", "CONT", "ASD"),
                         levels = c("CONT", "ASD"))


## Create signature gene set from mouse ensembl gene list
korb_signature = read_tsv("genelists/korbAllis_ASDGeneSig_Down_ensembl.txt", col_names = F)
DownSigGenes = korb_signature %>% deframe()

ensembl=useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
MMtoHG = getBM(attributes = c('ensembl_gene_id',
                              'hsapiens_homolog_ensembl_gene', 
                              'hsapiens_homolog_associated_gene_name'),
               mart = ensembl)

signatureIndex = rownames(geneRpkm2) %in% MMtoHG$hsapiens_homolog_ensembl_gene[MMtoHG$ensembl_gene_id %in% DownSigGenes]

## Generate eigengene
eigenGene = prcomp(t(geneRpkm2[signatureIndex,]))$x[,1]
pcaVars = getPcaVars(prcomp(t(log2(geneRpkm2[signatureIndex,]+1))))[1]

loads = prcomp(t(geneRpkm2[signatureIndex,]))$rot[,1]
names(loads) = MMtoHG$hsapiens_homolog_associated_gene_name[match(names(loads),
                                                                  MMtoHG$hsapiens_homolog_ensembl_gene)] 
boxplot(eigenGene ~ wrightASD_pd$Dx)

# clean up the eigenGene data
f = lm(eigenGene ~ Dx + as.numeric(rin) + 
         as.numeric(percentexonicmapping) + Sex, data=wrightASD_pd)

summary(f)

principalComps = prcomp(t(geneRpkm2))
fPC = lm(eigenGene ~ Dx + as.numeric(rin) + 
           as.numeric(percentexonicmapping) + Sex + principalComps$x[,1:9], data=wrightASD_pd)

mod = model.matrix(~Dx + as.numeric(rin) + 
                     as.numeric(percentexonicmapping) + Sex + principalComps$x[,1:9], data=wrightASD_pd)
ff = summary(lm(eigenGene ~ mod - 1))

cleanEigen = cleaningY(matrix(eigenGene, nr=1), mod, P=2)

# Regenerate plot from cleaned data
pdf("figures/2.wright2017_korbASDSigDown_deconvolution.pdf",h=6,w=5)
par(mar=c(5,6,2,2), cex.axis=2,cex.lab=2)
palette(brewer.pal(5,"Set1"))
boxplot(cleanEigen[1,] ~ wrightASD_pd$Dx, outline = FALSE,
        ylim = range(cleanEigen[1,]),xlab="",
        ylab = paste0("CAG Eigengene (Adj)"))
points(cleanEigen[1,] ~ jitter(as.numeric(factor(wrightASD_pd$Dx)), amount=0.15),
       pch = 21, bg=factor(wrightASD_pd$Dx))
legend("topright", paste0("p=", signif(ff$coef[2,4],3)),cex=1.5)
dev.off()

