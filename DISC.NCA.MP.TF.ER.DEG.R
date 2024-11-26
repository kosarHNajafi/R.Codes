The results are not satisfying applying DEG tha taking the sig ones and procedd only with them is like applying DEG twice.
Dropped
#---Load libraries-----------------
library(biomaRt)
library(RTN)

set.seed(123)
#---1.DEG for MP Genes-------------------

#METABRIC Sample Annotation Retrieved Manually
metabric_disc <- read.delim("~/NetworkComponentAnalysis/data/METABRIC Clinical.txt",header = TRUE,row.names = 1)
rownames(metabric_disc) <- gsub("\\-", "_", rownames(metabric_disc))
metabric_disc <- metabric_disc[order(rownames(metabric_disc)),]
View(metabric_disc)

DISC.MP <- read.delim("~/NCA/data/DISC.MP.Genes.txt",row.names = 1,header = TRUE)
colnames(DISC.MP) <- gsub("\\.", "_", colnames(DISC.MP))
View(DISC.MP)

DISC.MP.Only <- read.delim("~/NetworkComponentAnalysis/data/Disc.Just.MP.Genes.txt",header = TRUE,row.names = 1)
colnames(DISC.MP.Only) <- gsub("\\.", "_", colnames(DISC.MP.Only))
View(DISC.MP.Only)


disc_sample_mp_ids <- colnames(DISC.MP)
View(disc_sample_mp_ids)
length(disc_sample_mp_ids)

sort(disc_sample_mp_ids)
View(disc_sample_mp_ids)
length(disc_sample_mp_ids)

all.equal(disc_sample_mp_ids,colnames(DISC.MP))

# Find the common sample IDs between the two datasets
common_samples_mp_disc <- intersect(rownames(metabric_disc), disc_sample_mp_ids)
common_samples_mp_disc <- sort(common_samples_mp_disc)
View(common_samples_mp_disc)
length(common_samples_mp_disc)

# Subset and reorder both datasets to only include common samples
counts.mp.disc <- DISC.MP[, common_samples_mp_disc]
all.equal(rownames(DISC.MP), rownames(counts.mp.disc))
View(counts.mp.disc)
dim(counts.mp.disc)

sampleAnnotation.disc.mp <- metabric_disc[common_samples_mp_disc, ]
dim(sampleAnnotation.disc.mp)

# Verify that they are aligned
all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc.mp))  # Should return TRUE

all(colnames(counts.mp.disc) == rownames(sampleAnnotation.disc.mp))

# Create the design matrix for the linear model
# Assuming your label column is named "ER_status" with values "ERpos" and "ERneg"
design.mp <- model.matrix(~ 0 + factor(sampleAnnotation.disc$ER.Status))
colnames(design.mp) <- levels(factor(sampleAnnotation.disc$ER.Status))
rownames(design.mp) <- rownames(sampleAnnotation.disc)
View(design.mp)

# Fit the linear model using limma
fit.mp <- lmFit(counts.mp.disc, design.mp)

# Create contrast matrix to compare ERpos vs ERneg
contrast_matrix.mp <- makeContrasts(Positive_vs_Negative = Positive - Negative, levels = design.mp)

# Apply the contrast matrix
fit2.mp <- contrasts.fit(fit.mp, contrast_matrix.mp)

# Empirical Bayes adjustment
fit2.mp <- eBayes(fit2.mp)

# Extract results (log2 fold changes, p-values, etc.)
phenotype.mp <- topTable(fit2.mp, coef = "Positive_vs_Negative", adjust.method = "BH", number = Inf)

# Save results to a file if needed
write.table(phenotype.mp, file = file.path("~/NetworkComponentAnalysis/Output/","DEG_results.mp.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# Filter genes with significant adjusted p-values and a logFC threshold
hits.mp <- subset(phenotype.mp, adj.P.Val < 0.05 & abs(logFC) > 1)

# Save results to a file if needed
write.table(hits.mp, file = file.path("~/NetworkComponentAnalysis/Output/","DEG_sig_results.mp.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# View the differentially expressed genes
View(hits.mp)
dim(hits.mp)

#---2.DEG for TF Genes---------------

#Metabric_disc is the same
DISC.TF.Lmbert <- read.delim("~/NetworkComponentAnalysis/data/Disc.TF.Lambert.Curated.txt",header = TRUE, row.names = 2)
disc.tf <- DISC.TF.Lmbert[,-1]
colnames(disc.tf) <- gsub("\\.","_",colnames(disc.tf))
disc.tf <- disc.tf[order(rownames(disc.tf)),order(colnames(disc.tf))]
disc_sample_tf_ids <- colnames(disc.tf)
View(disc_sample_tf_ids)
length(disc_sample_mp_ids)

sort(disc_sample_tf_ids)
View(disc_sample_tf_ids)
length(disc_sample_tf_ids)

all.equal(disc_sample_tf_ids,colnames(disc.tf))

# Find the common sample IDs between the two datasets
common_samples_tf_disc <- intersect(rownames(metabric_disc), disc_sample_tf_ids)
common_samples_tf_disc <- sort(common_samples_tf_disc)
View(common_samples_tf_disc)
length(common_samples_tf_disc)

# Subset and reorder both datasets to only include common samples
counts.tf.disc <- disc.tf[, common_samples_tf_disc]
all.equal(rownames(disc.tf), rownames(counts.tf.disc))
View(counts.tf.disc)
dim(counts.tf.disc)

sampleAnnotation.disc.tf <- metabric_disc[common_samples_tf_disc, ]
dim(sampleAnnotation.disc.tf)

# Verify that they are aligned
all.equal(colnames(counts.tf.disc), rownames(sampleAnnotation.disc.tf))  # Should return TRUE

all(colnames(counts.tf.disc) == rownames(sampleAnnotation.disc.tf))

# Create the design matrix for the linear model
# Assuming your label column is named "ER_status" with values "ERpos" and "ERneg"
design.tf <- model.matrix(~ 0 + factor(sampleAnnotation.disc.tf$ER.Status))
colnames(design.tf) <- levels(factor(sampleAnnotation.disc.tf$ER.Status))
rownames(design.tf) <- rownames(sampleAnnotation.disc.tf)
View(design.tf)

# Fit the linear model using limma
fit.tf <- lmFit(counts.tf.disc, design.tf)

# Create contrast matrix to compare ERpos vs ERneg
contrast_matrix.tf <- makeContrasts(Positive_vs_Negative = Positive - Negative, levels = design.tf)

# Apply the contrast matrix
fit2.tf <- contrasts.fit(fit.tf, contrast_matrix.tf)

# Empirical Bayes adjustment
fit2.tf <- eBayes(fit2.tf)

# Extract results (log2 fold changes, p-values, etc.)
phenotype.tf <- topTable(fit2.tf, coef = "Positive_vs_Negative", adjust.method = "BH", number = Inf)

# Save results to a file if needed
write.table(phenotype.tf, file = file.path("~/NetworkComponentAnalysis/Output/","DEG_results.tf.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# Filter genes with significant adjusted p-values and a logFC threshold
hits.tf <- subset(phenotype.tf, adj.P.Val < 0.05 & abs(logFC) > 1)

# Save results to a file if needed
write.table(hits.tf, file = file.path("~/NetworkComponentAnalysis/Output/","DEG_sig_results.tf.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# View the differentially expressed genes
View(hits.tf)
dim(hits.tf)

DISC.MP.TF <- read.delim("~/NetworkComponentAnalysis/data/DISC.MP.TF.Curated.txt",header = TRUE, row.names = 1)
colnames(DISC.MP.TF) <- gsub("\\.", "_", colnames(DISC.MP.TF))
View(DISC.MP.TF)

sig.disc.mp.tf.Genes <- c(rownames(hits.mp),rownames(hits.tf))
disc.sig.mp.tf.common <- intersect(sig.disc.mp.tf.Genes,rownames(DISC.MP.TF))
disc.mp.tf.sig <- DISC.MP.TF[disc.sig.mp.tf.common,]
View(disc.mp.tf.sig)

#Sort data by columnnames(Sample IDs) and rownames(Gene IDs)
disc.mp.tf.sig <- disc.mp.tf.sig[order(rownames(disc.mp.tf.sig)),order(colnames(disc.mp.tf.sig))]
View(disc.mp.tf.sig)
head(disc.mp.tf.sig)
str(disc.mp.tf.sig)
dim(disc.mp.tf.sig)

print("Step 1 completed: Loading Metabolic and TF Gene expression in Discovery set")


#---2.Available METABRIC Clinical File as sampleAnnotation.disc =colData----------------------

Metabric_Manual_disc <- read.delim("~/NCA/data/METABRIC Clinical.txt", row.names = 1)
View(Metabric_Manual_disc)

Metabric_Manual_disc <- Metabric_Manual_disc[order(rownames(Metabric_Manual_disc)),]
View(Metabric_Manual_disc)

# Replace dots with underscores in column names if needed
rownames(Metabric_Manual_disc) <- gsub("\\-", "_", rownames(Metabric_Manual_disc))
View(Metabric_Manual_disc)

# Load required library
library(dplyr)

# Assuming your dataset is named Metabric_Manual_disc
# Create a new binary dataframe based on transformations
Metabric_Manual_disc <- Metabric_Manual_disc %>%
  transmute(
    IDs = rownames(Metabric_Manual_disc),
    Cohort = Cohort,
    
    OS.time = Overall.Survival..Months.,
    OS.event = ifelse(Overall.Survival.Status == "1:DECEASED", 1, 0),
    DSS.event = ifelse(Patient.s.Vital.Status == "Died of Disease", 1, 0),
    RFS.event = ifelse(Relapse.Free.Status == "1:Recurred", 1, 0),
    RFS.time = Relapse.Free.Status..Months.,
    
    Grade = Neoplasm.Histologic.Grade,
    Size = Tumor.Size,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    
    
    # Subtypes for LumA, LumB, Basal, Her2, Normal based on Pam50 subtype
    LumA = ifelse(Pam50...Claudin.low.subtype == "LumA", 1, 0),
    LumB = ifelse(Pam50...Claudin.low.subtype == "LumB", 1, 0),
    Basal = ifelse(Pam50...Claudin.low.subtype == "Basal", 1, 0),
    Her2 = ifelse(Pam50...Claudin.low.subtype == "Her2", 1, 0),
    Normal = ifelse(Pam50...Claudin.low.subtype == "Normal", 1, 0),
    
    # ER and PR status (positive and negative)
    `ER+` = ifelse(ER.Status == "Positive", 1, 0),
    `ER-` = ifelse(ER.Status == "Negative", 1, 0),
    
    
    # Histologic Grade categories G1, G2, G3
    G1 = ifelse(Neoplasm.Histologic.Grade == 1, 1, 0),
    G2 = ifelse(Neoplasm.Histologic.Grade == 2, 1, 0),
    G3 = ifelse(Neoplasm.Histologic.Grade == 3, 1, 0),
    
    # Hormone Therapy (HT)
    HT = ifelse(Hormone.Therapy == "YES", 1, 0),
    
  )
View(Metabric_Manual_disc)

# Print head of the colAnnotation for verification
head(Metabric_Manual_disc)

print("Step 2 completed: Manual METABRIC Clinical data")

#---3.disc_sample_mp.tf_ids Sample/Column Annotation-----------------------------
disc_sample_mp.tf_ids <- colnames(disc.mp.tf)
View(disc_sample_mp.tf_ids)
length(disc_sample_mp.tf_ids)

sort(disc_sample_mp.tf_ids)
View(disc_sample_mp.tf_ids)
length(disc_sample_mp.tf_ids)

all.equal(disc_sample_mp.tf_ids,colnames(disc.mp.tf))

# Find the common sample IDs between the two datasets
common_samples_disc <- intersect(rownames(Metabric_Manual_disc), disc_sample_mp.tf_ids)
common_samples_disc <- sort(common_samples_disc)
View(common_samples_disc)
length(common_samples_disc)

# Subset and reorder both datasets to only include common samples
counts.mp.tf.disc <- disc.mp.tf[, common_samples_disc]
all.equal(rownames(disc.mp.tf), rownames(counts.mp.tf.disc))
View(counts.mp.tf.disc)
dim(counts.mp.tf.disc)

sampleAnnotation.disc <- Metabric_Manual_disc[common_samples_disc, ]
dim(sampleAnnotation.disc)

# Verify that they are aligned
all.equal(colnames(counts.mp.tf.disc), rownames(sampleAnnotation.disc))  # Should return TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.mp.tf.disc), rownames(sample_metadata))
#sampleAnnotation.disc <- sample_metadata[tempindex, ]

#Check again

all.equal(colnames(counts.mp.tf.disc), rownames(sampleAnnotation.disc))

print("Step 3 completed: Discovery SampleAnnotation Prepared")

#---4.Load Ensemble Gene Annotation----------------------------------------------

#SHOULD I KEEP OR REMOVE?
#1.Accessing the data available in Ensembl by biomaRT
browseVignettes("biomaRt")

#2.Selecting an Ensembl BioMart database and dataset

# Connect to Ensembl database, and query human genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


#2.1##Step1: Identifying the database you need

#find the names of the BioMart services Ensembl is currently providing
##and Can be used to connect to the desired BioMart database:
listEnsembl()
#biomart argument should be given a valid name from the output of listEnsembl()
ensembl <- useEnsembl(biomart = "genes")


#2.2##Step 2: Choosing a dataset

#Within the Ensembl dataset each species is a different dataset.
#look at which datasets are available in the selected BioMart
#by using the function listDatasets()
datasets <- listDatasets(ensembl)
View(datasets)

##Because the listDatasets are so long
##in ensembl mart find anything that contain the term ‘hsapiens’

searchDatasets(mart = ensembl, pattern = "hsapiens")

#If you've been through these all before;select a both the database and dataset in one step
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)


#2.3##Ensembl mirror sites
#listMarts() to find the biomart
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")
#UP TO HERE SHALL I KEEP OR CAN I REMOVE IT???
#maintain consistent annotation throughout the duration of a project.
#2.4##Using archived versions of Ensembl, with no arguments
listEnsemblArchives()
listEnsembl(version = 113)
searchDatasets(mart = ensembl_113, pattern = "hsapiens")
#to ensure that script you write now will return exactly the same results in the future; 
##copy the URL from listEnsemblArchive() or www.ensembl.org in host Argument
ensembl_113 <- useEnsembl(biomart = "genes", 
                          dataset = "hsapiens_gene_ensembl",
                          host = "https://oct2024.archive.ensembl.org",
                          version = 113)

print("Step 8 completed: Ensembl Gene Information loaded")

#---5.disc_mp.tf Gene Annotation/build a biomaRt query ---------------------------
# Extract gene annotations for the corresponding genes from `ensembl_113`
filters = listFilters(ensembl_113)
View(filters)

View(listAttributes(ensembl_113))

gene_annot_mp.tf_disc <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name","description"),
  filters = "external_gene_name",
  values = rownames(disc.mp.tf),
  mart = ensembl_113
)

#WITH ENTREZ 3185 AND WITHOUT THAT 3152, DON'T MENTION "entrezgene_id".
View(gene_annot_mp.tf_disc)

# Check for duplicated external_gene_names
table(duplicated(gene_annot_mp.tf_disc$external_gene_name))

write.table(gene_annot_mp.tf_disc,file = file.path("~/NCA.ER/","Gene_annot_mp.tf.disc.txt"),sep = "\t",row.names = FALSE)

#Remove Duplicates Manually
gene_annot_mp.tf_disc <- read.delim("~/NCA.ER/NCA.MP.TF.DISC/GeneAnnotationCurated_Ensembl113.Used.txt")

dim(gene_annot_mp.tf_disc)
length(gene_annot_mp.tf_disc)
View(gene_annot_mp.tf_disc)
dim(counts.mp.tf.disc)

#Only Validated IDs Taken

# Extract gene IDs from both datasets
ensemble_ids <- gene_annot_mp.tf_disc$external_gene_name
View(ensemble_ids)

ensemble_ids <- sort(ensemble_ids)
length(ensemble_ids)

disc_mp.tf_gene_ids <- rownames(counts.mp.tf.disc)
View(disc_mp.tf_gene_ids)

common_genes <- intersect(ensemble_ids, disc_mp.tf_gene_ids) # Should match or be close to 1420
common_genes <- sort(common_genes)
all.equal(ensemble_ids, common_genes)
length(common_genes)

# Subset and reorder both datasets to only include common samples
counts.mp.tf.disc <- counts.mp.tf.disc[common_genes, ]
View(counts.mp.tf.disc)

all.equal(ensemble_ids,rownames(counts.mp.tf.disc))

#Set rownames for gene annotation, %n% doesn't let your data go NA 
gene_annot_mp.tf_disc <- gene_annot_mp.tf_disc[gene_annot_mp.tf_disc$external_gene_name %in% common_genes, ]
rownames(gene_annot_mp.tf_disc) <- gene_annot_mp.tf_disc$external_gene_name
View(gene_annot_mp.tf_disc)
all.equal(gene_annot_mp.tf_disc$external_gene_name, rownames(counts.mp.tf.disc))
all.equal(gene_annot_mp.tf_disc$external_gene_name,rownames(gene_annot_mp.tf_disc))
all.equal(gene_annot_mp.tf_disc$external_gene_name,common_genes)
all.equal(common_genes,ensemble_ids)

#Sort counts.mp.tf.disc and gene_annot_mp.tf_disc based on rownames which is gene names
counts.mp.tf.disc <- counts.mp.tf.disc[order(rownames(counts.mp.tf.disc)), ]
gene_annot_mp.tf_disc <- gene_annot_mp.tf_disc[order(rownames(gene_annot_mp.tf_disc)), ]

all.equal(rownames(counts.mp.tf.disc), rownames(gene_annot_mp.tf_disc))  # Should return TRUE if aligned
dim(counts.mp.tf.disc)
dim(gene_annot_mp.tf_disc)
View(gene_annot_mp.tf_disc)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(gene_annot_mp.tf_disc)
colnames(gene_annot_mp.tf_disc)[colnames(gene_annot_mp.tf_disc) == "external_gene_name"] <- "SYMBOL"
colnames(gene_annot_mp.tf_disc)[colnames(gene_annot_mp.tf_disc) == "ensembl_gene_id"] <- "ENSEMBL"

# Verify the change
colnames(gene_annot_mp.tf_disc)


# One final check:
stopifnot(rownames(gene_annot_mp.tf_disc) == rownames(counts.mp.tf.disc), # features
          rownames(sampleAnnotation.disc) == colnames(counts.mp.tf.disc)) # samples

print("Step 5 completed: GeneAnnotation/RowAnnotation done for discovery set")

#---6.Create disc_mp.tf list/SummarizedExperiment-----------------------
counts.mp.tf.disc <- as.matrix(counts.mp.tf.disc)
View(counts.mp.tf.disc)
str(counts.mp.tf.disc)

disc_mp.tf <- list(
  expData = counts.mp.tf.disc,                       # expData as an assay
  rowAnnotation = gene_annot_mp.tf_disc,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.disc           # Sample annotations (column metadata)
)

View(disc_mp.tf)
###DEBUGGING ERRORS###
#> rtni_disc_mp.tf.NUM <- tni.constructor(expData = disc_mp.tf.NUM$expData, 
#                                     +                              regulatoryElements = regulatoryElements, 
#                                     +                              rowAnnotation = disc_mp.tf.NUM$rowAnnotation, 
#                                     +                              colAnnotation = disc_mp.tf.NUM$colAnnotation)
#-Preprocessing for input data...
#Error: 'expData' should be a numeric matrix with genes on rows and 
#samples on cols!

#Values are not numeric

# Find elements that will turn into NA when coerced to numeric
non_numeric_elements <- counts.mp.tf.disc[is.na(as.numeric(counts.mp.tf.disc))]

# Display non-numeric elements to understand what they are
print(non_numeric_elements)

#---TO RETURN NUMERIC-BUT IT GIVES THE ERROR:
#Warning message:
#  In matrix(as.numeric(counts.mp.tf.disc), nrow = nrow(counts.mp.tf.disc),  :
#             NAs introduced by coercion

#SummarizedExperiment of rtni_disc_mp.tf.NUM_mp.tf_se
#rtni_disc_mp.tf.NUM_mp.tf_se <- SummarizedExperiment(
#  assays = counts.mp.tf.disc,
#  rowData = gene_annot_mp.tf_disc,
#  colData = sampleAnnotation.disc
#)

#View(rtni_disc_mp.tf.NUM_mp.tf_se)

# Check dimensions and alignment
#dim(counts.mp.tf.disc.NUM)

#dim(disc_mp.tf.NUM_se$expData)  # Should match dimensions of counts.mp.tf.disc
all.equal(colnames(disc_mp.tf.NUM$expData), rownames(sampleAnnotation.disc))  # Should return TRUE
all.equal(disc_mp.tf.NUM$expData, counts.mp.tf.disc)
all.equal(disc_mp.tf.NUM$colAnnotation, sampleAnnotation.disc)

all.equal(disc_mp.tf$rowAnnotation,disc_mp.tf.NUM$rowAnnotation)
#[1] TRUE
all.equal(disc_mp.tf$colAnnotation,disc_mp.tf.NUM$colAnnotation)
#[1] TRUE
all.equal(disc_mp.tf$expData,disc_mp.tf.NUM$expData)
#[1] "Modes: character, numeric"                        
#[2] "'is.NA' value mismatch: 14 in current 0 in target"
###DEBUGGING ERRORS###

print("Step 6 completed: list/SummarizedExperiment file prepared for Discovery")

#---7.RegulatoryElements-----------------

# Load TF annotation
data("tfsData")

# Check TF annotation:
# Intersect TFs from Lambert et al. (2018) with gene annotation 
# from the gene expression of 90 metabolic pathway cohort
regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, disc_mp.tf$rowAnnotation$SYMBOL)
View(regulatoryElements)
regulatoryElements <- sort(regulatoryElements)
View(regulatoryElements)

print("Step 7 completed: RegulatoryElements Defined")

#---8.Run the TNI constructor with the extracted matrix for disc_mp.tf.NUM--------------
#This dataset consists of a list with 3 objects:
##a named gene expression matrix (tniData$expData),
##a data frame with gene annotations (tniData$rowAnnotation), 
##and a data frame with sample annotations (tniData$colAnnotation).
##alternatively, 'expData' can be a 'SummarizedExperiment' object
rtni_disc_mp.tf <- tni.constructor(expData = disc_mp.tf$expData, 
                             regulatoryElements = regulatoryElements, 
                             rowAnnotation = disc_mp.tf$rowAnnotation, 
                             colAnnotation = disc_mp.tf$colAnnotation)

save(rtni_disc_mp.tf, file = file.path("~/NCA.ER/NCA.MP.TF.DISC.NCA/", "rtni_mp.tf.disc.RData"))

###snow clustering returns zero results
### Compute the reference regulatory network by permutation and bootstrap analyses.
### Please set 'spec' according to your available hardware
#options(cluster=snow::makeCluster(spec=4, "SOCK"))  #???

#1.with pValuCutoff =1e-7 it gives 0 values, even in the summary
#rtni_disc_mp.tf.NUM <- tni.permutation(rtni_disc_mp.tf.NUM, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#2.with pValuCutoff= 1e-7 and npermutation = 1000 it gives 0 values, even in the sumary
#rtni_disc_mp.tf.NUM <- tni.permutation(rtni_disc_mp.tf.NUM,nPermutations = 1000, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#3.  nPermutations >= 1000 with running snow goes the same zero values
#4.  nPermutations >= 1000 without snow package, and pValueCutoff works!
rtni_disc_mp.tf <- tni.permutation(rtni_disc_mp.tf, nPermutations = 1000) 

#Unstable interactions are subsequently removed by bootstrap analysis,
##creates a consensus bootstrap network, referred here as refnet (reference network).
rtni_disc_mp.tf <- tni.bootstrap(rtni_disc_mp.tf)

#stopCluster(getOption("cluster"))  # avoid leak memory

# Compute the DPI-filtered regulatory network
rtni_disc_mp.tf <- tni.dpi.filter(rtni_disc_mp.tf, eps = NA)
tni.regulon.summary(rtni_disc_mp.tf)

#detailed information about a specific regulon
tni.regulon.summary(rtni_disc_mp.tf, regulatoryElements = "SATB1")

# Save the TNI object for subsequent analyses
save(rtni_disc_mp.tf, file = file.path("~/NCA.ER/NCA.MP.TF.DISC.NCA/", "rtni_disc_mp.tf.RData"))

print("Step 8 completed: Loaded rtni data and saved.")

regulons <- tni.get(rtni_disc_mp.tf, what = "regulons.and.mode", idkey = "SYMBOL")
View(regulons)
head(regulons)

help("tni.get")
#---9. Regulon activity profiles----------

library(Fletcher2013b)
library(pheatmap)
library(grid)
library(gridExtra)

# Load 'rtni1st' data object, which includes regulons and expression profiles
#here I will use: rtni_disc_mp.tf.NUM

# A list of transcription factors of interest ( #9 )
# Compute regulon activity for individual samples
rtni_disc_mp.tf <- tni.gsea2(rtni_disc_mp.tf, regulatoryElements = rtni_disc_mp.tf@regulatoryElements)
metabric_regact_disc <- tni.get(rtni_disc_mp.tf, what = "regulonActivity")
View(metabric_regact_disc)

# Get sample attributes from the 'rtni_disc_mp.tf.NUM' dataset
metabric_annot_disc <- tni.get(rtni_disc_mp.tf, "colAnnotation")

# Get ER+/- and PAM50 attributes for pheatmap
attribs_disc <- c("ER+","ER-")
metabric_annot_disc <- metabric_annot_disc[,attribs_disc]

save(rtni_disc_mp.tf, file = file.path("~/NCA.ER/NCA.MP.TF.DISC.NCA/", "rtni_disc_mp.tf.gsea.RData"))

# Step 1: Identify samples with ER+ = 1
ER_positive_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 1]
ER_negative_samples <- rownames(metabric_annot_disc)[metabric_annot_disc$`ER+` == 0]

# Step 2: Order columns in metabric_regact_disc with ER+ samples on the left
ordered_sample_names <- c(ER_positive_samples, ER_negative_samples)
metabric_regact_disc_ordered <- metabric_regact_disc$differential[ordered_sample_names,]

# Define custom colors for each category in annotation_col
ER_annotation <- metabric_annot_disc[,c("ER+","ER-")]
ER_annotation_colors <- list(
  "ER+" = c("0" = "lightgrey", "1" = "blue"),
  "ER-" = c("0" = "lightgrey", "1" = "red")
  )

pdf("~/NCA.ER/NCA.MP.TF.DISC.NCA/Disc.NCA.MP.TF.Heatmap.pdf", width = 10, height = 10)

# Plot regulon activity profiles
disc.heatmap <- pheatmap(t(metabric_regact_disc_ordered), 
                         main="Discovery Set (n=988 samples)",
                         annotation_col = ER_annotation,
                         show_colnames = FALSE, annotation_legend = FALSE, 
                         clustering_method = "ward.D2", fontsize_row = 3,
                         clustering_distance_rows = "correlation",
                         cluster_cols = FALSE,
                         legend = TRUE,
                         annotation_colors = ER_annotation_colors,
                         fontsize_col = 3, fontsize = 6,
                         border_color = NA,
)

grid.text("Regulons", x= 0.97 , y=0.3, rot=270)

dev.off()

disc.session.info <- sessionInfo()

#---10.Transcriptional Network Analysis (TNA)-----------
Metabric_Manual_disc_tna <- read.delim("~/NCA/data/METABRIC Clinical.txt", row.names = 1)
View(Metabric_Manual_disc_tna)

Metabric_Manual_disc_tna <- Metabric_Manual_disc_tna[order(rownames(Metabric_Manual_disc_tna)),]
View(Metabric_Manual_disc_tna)

# Replace dots with underscores in column names if needed
rownames(Metabric_Manual_disc_tna) <- gsub("\\-", "_", rownames(Metabric_Manual_disc_tna))
View(Metabric_Manual_disc_tna)

disc_sample_mp.tf_ids_tna <- colnames(disc.mp.tf)
View(disc_sample_mp.tf_ids_tna)
length(disc_sample_mp.tf_ids_tna)

sort(disc_sample_mp.tf_ids_tna)
View(disc_sample_mp.tf_ids_tna)
length(disc_sample_mp.tf_ids_tna)

all.equal(disc_sample_mp.tf_ids_tna,colnames(disc.mp.tf))

# Find the common sample IDs between the two datasets
common_samples_disc_tna <- intersect(rownames(Metabric_Manual_disc_tna), disc_sample_mp.tf_ids_tna)
common_samples_disc_tna <- sort(common_samples_disc_tna)
View(common_samples_disc_tna)
length(common_samples_disc_tna)

sampleAnnotation.disc_tna <- Metabric_Manual_disc_tna[common_samples_disc_tna, ]
dim(sampleAnnotation.disc_tna)

# Verify that they are aligned
all.equal(colnames(counts.mp.tf.disc), rownames(sampleAnnotation.disc_tna))  # Should return TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.mp.tf.disc), rownames(sample_metadata))
#sampleAnnotation.disc <- sample_metadata[tempindex, ]

#Check again

all.equal(colnames(counts.mp.tf.disc), rownames(sampleAnnotation.disc_tna))
all(disc_mp.tf$colAnnotation$IDs == rownames(sampleAnnotation.disc_tna))

# Check that row names of sample_annotations match column names of expression_data
if (!all(rownames(sampleAnnotation.disc_tna) == colnames(disc_mp.tf$expData))) {
  stop("Mismatch between sample annotation rownames and expression data colnames!")
}

# Create the design matrix for the linear model
# Assuming your label column is named "ER_status" with values "ERpos" and "ERneg"
design <- model.matrix(~ 0 + factor(sampleAnnotation.disc_tna$ER.Status))
colnames(design) <- levels(factor(sampleAnnotation.disc_tna$ER.Status))
rownames(design) <- rownames(sampleAnnotation.disc_tna)
View(design)
all.equal(as.vector(design[,"Positive"]),as.vector(sampleAnnotation.disc$`ER+`))
all(design[,"Positive"] == sampleAnnotation.disc$`ER+`)
all(design[,"Negative"] == sampleAnnotation.disc$`ER-`)

# Fit the linear model using limma
fit <- lmFit(counts.mp.tf.disc, design)

# Create contrast matrix to compare ERpos vs ERneg
contrast_matrix <- makeContrasts(Positive_vs_Negative = Positive - Negative, levels = design)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes adjustment
fit2 <- eBayes(fit2)

# Extract results (log2 fold changes, p-values, etc.)
phenotype <- topTable(fit2, coef = "Positive_vs_Negative", adjust.method = "BH", number = Inf)

# Save results to a file if needed
write.table(phenotype, file = "differential_expression_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# Output the top results for inspection
View(phenotype)
dim(phenotype)

# Order phenotype alphabetically by row names
phenotype <- phenotype[order(rownames(phenotype)), ]

all(rownames(phenotype) == rownames(counts.mp.tf.disc))
all(rownames(phenotype) == gene_annot_mp.tf_disc$SYMBOL)

# Extract 'logFC' as a named numeric vector
logFC_disc_mp.tf <- setNames(phenotype$logFC, rownames(phenotype))
View(logFC_disc_mp.tf)

all.equal(names(logFC_disc_mp.tf),rownames(phenotype))
all.equal(as.vector(logFC_disc_mp.tf),phenotype$logFC)

# Filter genes with significant adjusted p-values and a logFC threshold
hits <- subset(phenotype, adj.P.Val < 0.05 & abs(logFC) > 1)

# View the differentially expressed genes
View(hits)
dim(hits)

all(rownames(phenotype) == gene_annot_mp.tf_disc$SYMBOL)

#Error: NOTE: all names in 'phenotype' should be available in col1 of 'phenoIDs'!
library(dplyr)
gene_annot_mp.tf_disc_tna <- gene_annot_mp.tf_disc %>%
  select(SYMBOL, everything())

View(gene_annot_mp.tf_disc_tna)
all(gene_annot_mp.tf_disc$SYMBOL == gene_annot_mp.tf_disc_tna$SYMBOL)
all(gene_annot_mp.tf_disc$ENSEMBL == gene_annot_mp.tf_disc_tna$ENSEMBL)

tna.disc_mp.tf <- list(
  phenotype = logFC_disc_mp.tf,
  phenoID = gene_annot_mp.tf_disc_tna,
  hits = rownames(hits)
)

View(tna.disc_mp.tf)
# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype

#set.seed(123)  Since I didnt run everything from scratch for this part I used set.seed() in saved Rproject

#CHECK "rtnaData"
rtna_disc.mp.tf <- tni2tna.preprocess(object = rtni_disc_mp.tf, 
                           phenotype = tna.disc_mp.tf$phenotype, 
                           hits = tna.disc_mp.tf$hits, 
                           phenoIDs = tna.disc_mp.tf$phenoID)
# Run the MRA method
rtna_disc.mp.tf <- tna.mra(rtna_disc.mp.tf)

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_disc.mp.tf <- tna.get(rtna_disc.mp.tf, what="mra", ntop = -1)
View(mra_disc.mp.tf)

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_disc.mp.tf <- tna.gsea1(rtna_disc.mp.tf, nPermutations=1000)

# Get GSEA results
gsea1_disc.mp.tf <- tna.get(rtna_disc.mp.tf, what="gsea1", ntop = -1)
head(gsea1_disc.mp.tf)
View(gsea1_disc.mp.tf)

# Filter for significant TFs
gsea1_disc.mp.tf.sig <- gsea1_disc.mp.tf[gsea1_disc.mp.tf$Adjusted.Pvalue <= 0.05, ]
#View(gsea1_disc.mp.tf.sig)

# Specify TFs to include in the plot
sig_tfs_disc_mp.tf <- rownames(gsea1_disc.mp.tf.sig)

# Plot GSEA results for significant TFs
tna.plot.gsea1(
  rtna_disc.mp.tf,
  labPheno = "abs(log2 fold changes)",
  tfs = sig_tfs_disc_mp.tf, # Include only significant TFs
  file = "gsea1_plot_disc_mp.tf_sig.pdf",
  filepath = "~/NCA.ER/NCA.MP.TF.DISC.NCA/"
)

# Plot GSEA results
tna.plot.gsea1(rtna_disc.mp.tf,
               labPheno="abs(log2 fold changes)", 
               ntop = 5,
               file = "gsea1_plot_top5", 
               ylimPanels = c(0.0,3,0.0,0.0),
               heightPanels = c(1,1,3))


#---11.Extract Regulons-------------------------

writeLines(capture.output(sessionInfo()), "sessionInfo_disc_mp.tf.txt")

write.table(regulons$BCL11A,file = file.path("~/NCA.ER/NCA.MP.TF.DISC.NCA/Regulons/","BCL11A.Regulon.txt"),sep = "\t")
# Load necessary library to get all the regulons
library(dplyr)

# Assuming `regulons` is a list where each element is named by the regulon (e.g., `regulons$DNMT1`)
for (regulon_name in names(regulons)) {
  regulon_data <- regulons[[regulon_name]]  # Access each regulon as a vector
  
  # Check if the regulon has more than one gene
  if (length(regulon_data) > 1) {
    # Prepare a data frame for output, including the gene names and values
    output_df <- data.frame(Gene = names(regulon_data), Value = regulon_data)
    
    # Prepare the output file path
    output_file <- paste0("~/NCA.ER/NCA.MP.TF.DISC.NCA/Regulons/", regulon_name, ".DISC.txt")
    
    # Write the regulon data to a file
    write.table(output_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}
