library(maftools)

snv_dat <- read.delim("sample.maf", header=T) 

#read the mutation file in to build a maf file that maftools oncoplot can recgnize as input 
#rearrange the columns so that maftools can recognize the input

colnames(snv_dat) <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1","Hugo_Symbol","Variant_Classification", "Amino_Acid_Change","DNA_VAF","SIFT", "PolyPhen") # this has to be adapted to the dataset
snv_dat$NCBI_Build <- "GRCh37"
snv_dat$Tumor_Seq_Allele2 <- snv_dat$Tumor_Seq_Allele1
snv_dat$Variant_Type <- "SNP"
head(snv_dat)
cols <- c("Hugo_Symbol", "NCBI_Build", "Variant_Classification", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Variant_Type")
snv_dat <- snv_dat[,cols]

write.table(snv_dat, "sample_mod.maf", sep = "\t", quote=F, row.names=F)

clinical <- read.delim("clinical_classification.tsv", sep = "\t")
clinical <- clinical[c('sampleID', 'X_cohort')]
colnames(clinical) <- c('Tumor_Sample_Barcode', "Disease")
head(clinical)

maf_dat <- read.maf(maf = "sample_mod.maf", vc_nonSyn = c("RNA"), clinicalData = clinical)
#it is important that even if there are no non-synonymous mutations (for ex working with non-coding mutations), the mutations must be labelled as non-synonymous for maftools to plot them
#in this case I'm labelling variants classified as 'RNA' as non-synonymous

vc_cols = RColorBrewer::brewer.pal(n=5, name='Dark2')
names(vc_cols) = c('RNA')

oncoplot(maf=maf_dat, top = 10, colors = vc_cols, draw_titv = TRUE, gene_mar = 6, fontSize = 0.1, sepwd_genes=-1, sepwd_samples=-0.1, removeNonMutated = TRUE,  clinicalFeatures=c('Disease'))

