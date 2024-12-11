#---Install Required Packages---------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install biomaRt if not already installed(for gene annotation section)
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}

BiocManager::install("GenomicRanges")
#---Load libraries-----------------
library(GenomicRanges)
library(SummarizedExperiment)
library(MetaGxBreast)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(RTN)
library(snow)

set.seed(123)
###################################################################################
#---Metabolic pathway genes without it's overlap with TF as a seperate disc_mp[list] and TF gene expressions as disc_tf[list]-------
#---1.Load Gene Expression Data of 90 metabolic pathways = a matrix for counts/assays---------------------------------

#The mutal genes with TFs are being deleted
#load txt. Format
DISC.MP.Genes <- read.delim("~/NCA.ER/Disc.JustMP.Annotated2/Disc.Just.MP.Genes.txt")
View(DISC.MP.Genes)

#Preparing disc.mp.genes 
disc.mp.genes <- DISC.MP.Genes[,-1]
rownames(disc.mp.genes) <- DISC.MP.Genes[,1]

# Replace dots with underscores in column names if needed
colnames(disc.mp.genes) <- gsub("\\.", "_", colnames(disc.mp.genes))

#Sort data by columnnames(Sample IDs) and Rownames(Gene Symbols)
disc.mp.genes <- disc.mp.genes[,order(colnames(disc.mp.genes))]
disc.mp.genes <- disc.mp.genes[order(rownames(disc.mp.genes)),]
View(disc.mp.genes)
str(disc.mp.genes)
dim(disc.mp.genes)
print("Step 1 completed: Loading Metabolic Gene expression in Discovery set")

#---2.Available METABRIC Clinical File as sampleAnnotation.disc = colData----------------------
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

# Print head of the colAnnotation for verification
head(Metabric_Manual_disc)

print("Step 3 completed: Manual METABRIC Clinical data")

#---3.disc_mp Sample Annotation-----------------------------
disc_mp_sample_ids <- colnames(disc.mp.genes)
View(disc_mp_sample_ids)
length(disc_mp_sample_ids)

sort(disc_mp_sample_ids)
View(disc_mp_sample_ids)
length(disc_mp_sample_ids)

# Find the common sample IDs between the two datasets
common_samples_disc <- intersect(rownames(Metabric_Manual_disc), disc_mp_sample_ids)
common_samples_disc <- sort(common_samples_disc)
View(common_samples_disc)
length(common_samples_disc)

#
all.equal(colnames(disc.mp.genes),colnames(disc.tf.genes))

# Subset and reorder both datasets to only include common samples
#1:
counts.mp.disc <- disc.mp.genes[, common_samples_disc]
all.equal(rownames(disc.mp.genes), rownames(counts.mp.disc))
View(counts.mp.disc)
dim(counts.mp.disc)

#2:
counts.tf.disc <- disc.tf.genes[,common_samples_disc]
all.equal(rownames(disc.tf.genes),rownames(counts.tf.disc))
View(counts.tf.disc)
dim(counts.tf.disc)


all.equal(colnames(counts.mp.disc),colnames(counts.tf.disc))

#SampleAnnotation:
sampleAnnotation.disc <- Metabric_Manual_disc[common_samples_disc, ]
dim(sampleAnnotation.disc)

# Verify that they are aligned
all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc))  # Should return TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.mp.disc), rownames(sample_metadata))
#sampleAnnotation.disc <- sample_metadata[tempindex, ]

#Check again

all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc))

print("Step 4 completed: Discovery SampleAnnotation Prepared")

#---4.Load Ensembl Gene Annotation----------------------------------------------
#1.Accessing the data available in Ensembl by biomaRT
#2.Selecting an Ensembl BioMart database and dataset

# Connect to Ensembl database, and query human genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


#2.1##Step1: Identifying the database you need
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
ensembl

#2.2##Step 2: Choosing a dataset
#we look at which datasets are available in the selected BioMart
#by using the function listDatasets()
datasets <- listDatasets(ensembl)

##if we want to find the details of any datasets
##in our ensembl mart that contain the term ‘hsapiens’
##we could do the following
searchDatasets(mart = ensembl, pattern = "hsapiens")

ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)


#2.3##Ensembl mirror sites

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

#2.4##Using archived versions of Ensembl
listEnsemblArchives()
listEnsembl(version = 113)
ensembl_113 <- useEnsembl(biomart = "genes", 
                          dataset = "hsapiens_gene_ensembl",
                          host = "https://oct2024.archive.ensembl.org",
                          version = 113)

View(listAttributes(ensembl_113))

print("Step 4 completed: Ensembl Gene Information loaded")

#---5.disc_mp Gene Annotation ---------------------------
# Extract gene annotations for the corresponding genes from `ensembl_113`
gene_annot_mp_disc <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand","description"),
  filters = "external_gene_name",
  values = rownames(disc.mp.genes),
  mart = ensembl_113
)

View(gene_annot_mp_disc)

table(duplicated(gene_annot_mp_disc$external_gene_name))

#Since it has Duplicated values, we will select the correct ID manually and then read into RStudio again
write.table(gene_annot_mp_disc,file = file.path("~/NCA.ER/","Gene_Ensembl_113_Retrieved_NO.TF.txt"),sep = "\t",row.names = FALSE)
gene_annot_mp_disc <- read.delim("~/NCA.ER/Gene_Ensembl_113_Used_NO.TF.txt")
View(gene_annot_mp_disc)

# Save gene_annot_mp_disc as a txt file
write.table(gene_annot_mp_disc, file = file.path("~/NCA.ER/Disc.JustMP.Annotated2/","gene_annot_mp_disc_NO.TF_113.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Extract gene IDs from both datasets
ensemble_ids <- gene_annot_mp_disc$external_gene_id
View(ensemble_ids)

ensemble_ids <- sort(ensemble_ids)
length(ensemble_ids)

disc_mp_gene_ids <- rownames(counts.mp.disc)
View(disc_mp_gene_ids)

common_genes <- intersect(ensemble_ids, disc_mp_gene_ids) # Should match or be close to 1420
common_genes <- sort(common_genes)
all.equal(ensemble_ids, common_genes)
length(common_genes)

# Subset and reorder both datasets to only include common samples
counts.mp.disc <- counts.mp.disc[common_genes, ]
View(counts.mp.disc)

all.equal(ensemble_ids,rownames(counts.mp.disc))

gene_annot_mp_disc <- gene_annot_mp_disc[gene_annot_mp_disc$external_gene_name %in% common_genes, ]
rownames(gene_annot_mp_disc) <- gene_annot_mp_disc$external_gene_name
View(gene_annot_mp_disc)
all.equal(gene_annot_mp_disc$external_gene_name, rownames(counts.mp.disc))

#Sort counts.mp.disc and gene_annot_mp_disc based on rownames which is gene names
counts.mp.disc <- counts.mp.disc[order(rownames(counts.mp.disc)), ]
gene_annot_mp_disc <- gene_annot_mp_disc[order(rownames(gene_annot_mp_disc)), ]

all.equal(rownames(counts.mp.disc), rownames(gene_annot_mp_disc))  # Should return TRUE if aligned
dim(counts.mp.disc)
dim(gene_annot_mp_disc)
View(gene_annot_mp_disc)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(gene_annot_mp_disc)

# Suppose you have a data frame named 'df' with a column called 'external_gene_id'
colnames(gene_annot_mp_disc)[colnames(gene_annot_mp_disc) == "external_gene_name"] <- "SYMBOL"
colnames(gene_annot_mp_disc)[colnames(gene_annot_mp_disc) == "ensembl_gene_id"] <- "ENSEMBL"

# Verify the change
colnames(gene_annot_mp_disc)


# One final check:
stopifnot(rownames(gene_annot_mp_disc) == rownames(counts.mp.disc), # features
          rownames(sampleAnnotation.disc) == colnames(counts.mp.disc)) # samples

print("Step 5 completed: GeneAnnotation/RowAnnotation done for discovery set")

#---6.Create disc_mp list/SummarizedExperiment-----------------------
counts.mp.disc <- as.matrix(counts.mp.disc)
View(counts.mp.disc)
str(counts.mp.disc)

disc_mp <- list(
  expData = counts.mp.disc,                       # expData as an assay
  rowAnnotation = gene_annot_mp_disc,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.disc           # Sample annotations (column metadata)
)

View(disc_mp)

#SummarizedExperiment of rtni_disc_mp_se
#rtni_disc_mp_se <- SummarizedExperiment(
#  assays = counts.mp.disc,
#  rowData = gene_annot_mp_disc,
#  colData = sampleAnnotation.disc
#)

#View(rtni_disc_mp_se)

# Check dimensions and alignment
dim(counts.mp.disc)
#dim(disc_mp_se$expData)  # Should match dimensions of counts.mp.disc
all.equal(colnames(disc_mp$expData), rownames(sampleAnnotation.disc))  # Should return TRUE
all.equal(disc_mp$expData, counts.mp.disc)
all.equal(disc_mp$colAnnotation, sampleAnnotation.disc)

print("Step 6 completed: list/SummarizedExperiment file prepared for Discovery")

#---7.RegulatoryElements: Load Gene Expression Data of Lambert TFs-----------------

# Load TF annotation
data("tfsData")

#Manually aligned the expression of Lamber TF genes and loaded in R as a txt. Format
DISC.TF.Genes <- read.delim("~/NCA.ER/Disc.JustMP.Annotated2/DISC.TF.Genes.Lambert.txt")
View(DISC.TF.Genes)

#Preparing disc.tf.genes 
disc.tf.genes <- DISC.TF.Genes[,-c(1,2)]
rownames(disc.tf.genes) <- DISC.TF.Genes[,2]

# Replace dots with underscores in column names if needed
colnames(disc.tf.genes) <- gsub("\\.", "_", colnames(disc.tf.genes))

#Sort data by columnnames(Sample IDs) and Rownames(Gene Symbols)
disc.tf.genes <- disc.tf.genes[,order(colnames(disc.tf.genes))]
disc.tf.genes <- disc.tf.genes[order(rownames(disc.tf.genes)),]
View(disc.tf.genes)
str(disc.tf.genes)
dim(disc.tf.genes)

all.equal(colnames(counts.mp.disc),colnames(counts.tf.disc))

gene_annot_tf_disc <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand","description"),
  filters = "external_gene_name",
  values = rownames(disc.tf.genes),
  mart = ensembl_113
)

table(duplicated(gene_annot_tf_disc$external_gene_name))

#Manually Selecting the Unique IDs
write.table(gene_annot_tf_disc,file = file.path("~/NCA.ER/Disc.Justtf.Annotated2/","Gene_Ensembl_113_Retrieved_TF.txt"),sep = "\t",row.names = FALSE)
gene_annot_tf_disc <- read.delim("~/NCA.ER/Disc.Justtf.Annotated2/Gene_Ensembl_113_Used_TF.txt")

table(duplicated(gene_annot_tf_disc$external_gene_name))

# Extract gene IDs from both datasets
ensemble_ids_tf <- gene_annot_tf_disc$external_gene_name
View(ensemble_ids_tf)

ensemble_ids_tf <- sort(ensemble_ids_tf)
length(ensemble_ids_tf)

disc_tf_gene_ids <- rownames(counts.tf.disc)
View(disc_tf_gene_ids)

common_genes_tf <- intersect(ensemble_ids_tf, disc_tf_gene_ids) # Should match or be close to 1420
common_genes_tf <- sort(common_genes_tf)
all.equal(ensemble_ids_tf, common_genes_tf)
length(common_genes)

# Subset and reorder both datasets to only include common satfles
counts.tf.disc <- counts.tf.disc[common_genes_tf, ]
View(counts.tf.disc)

all.equal(ensemble_ids_tf,rownames(counts.tf.disc))

gene_annot_tf_disc <- gene_annot_tf_disc[gene_annot_tf_disc$external_gene_name %in% common_genes_tf, ]
rownames(gene_annot_tf_disc) <- gene_annot_tf_disc$external_gene_name
View(gene_annot_tf_disc)
all.equal(gene_annot_tf_disc$external_gene_name, rownames(counts.tf.disc))

#Sort counts.tf.disc and gene_annot_tf_disc based on rownames which is gene names
counts.tf.disc <- counts.tf.disc[order(rownames(counts.tf.disc)), ]
gene_annot_tf_disc <- gene_annot_tf_disc[order(rownames(gene_annot_tf_disc)), ]

all.equal(rownames(counts.tf.disc), rownames(gene_annot_tf_disc))  # Should return TRUE if aligned
dim(counts.tf.disc)
dim(gene_annot_tf_disc)
View(gene_annot_tf_disc)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(gene_annot_tf_disc)

# Suppose you have a data frame named 'df' with a column called 'external_gene_id'
colnames(gene_annot_tf_disc)[colnames(gene_annot_tf_disc) == "external_gene_name"] <- "SYMBOL"
colnames(gene_annot_tf_disc)[colnames(gene_annot_tf_disc) == "ensembl_gene_id"] <- "ENSEMBL"

# Verify the change
colnames(gene_annot_tf_disc)


counts.tf.disc <- as.matrix(counts.tf.disc)
View(counts.tf.disc)
str(counts.tf.disc)

disc_tf <- list(
  expData = counts.tf.disc,                       # expData as an assay
  rowAnnotation = gene_annot_tf_disc,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.disc           # Satfle annotations (column metadata)
)

View(disc_tf)

# Check dimensions and alignment
dim(counts.tf.disc)

#dim(disc_tf_se$expData)  # Should match dimensions of counts.tf.disc
all.equal(colnames(disc_tf$expData), rownames(sampleAnnotation.disc))  # Should return TRUE
all.equal(disc_tf$expData, counts.tf.disc)
all.equal(disc_tf$colAnnotation, sampleAnnotation.disc)

# One final check:
stopifnot(rownames(gene_annot_tf_disc) == rownames(counts.tf.disc), # features
          rownames(sampleAnnotation.disc) == colnames(counts.tf.disc)) # satfles

print("Step 7 completed: GeneAnnotation/RowAnnotation done for discovery set")

#---8.Run the TNI constructor with the extracted matrix for disc_mp--------------
#This dataset consists of a list with 3 objects:
##a named gene expression matrix (tniData$expData),
##a data frame with gene annotations (tniData$rowAnnotation), 
##and a data frame with sample annotations (tniData$colAnnotation).
##alternatively, 'expData' can be a 'SummarizedExperiment' object
rtni_disc <- tni.constructor(expData = disc_mp$expData, 
                             regulatoryElements = rownames(counts.tf.disc), 
                             rowAnnotation = disc_tf$rowAnnotation, 
                             colAnnotation = disc_mp$colAnnotation)

#-Preprocessing for input data...
#--Mapping 'expData' to 'rowAnnotation'...
#Error: all rownames in the expression data matrix should be available
#either in rownames or col1 of the row annotation!

#DUE TO ERROR THE PROCESS STOPPED HERE
#############################################################################################################################
#---1.Load Gene Expression Data of 90 metabolic pathways = a mtrix for counts/assays---------------------------------
#load txt. Format
DISC.MP.Genes <- read.delim("~/NCA/data/DISC.MP.Genes.txt")
View(DISC.MP.Genes)

#Preparing disc.mp.genes 
disc.mp.genes <- DISC.MP.Genes[,-1]
rownames(disc.mp.genes) <- DISC.MP.Genes[,1]

# Replace dots with underscores in column names if needed
colnames(disc.mp.genes) <- gsub("\\.", "_", colnames(disc.mp.genes))

#Sort data by columnnames(Sample IDs)
disc.mp.genes <- disc.mp.genes[,order(colnames(disc.mp.genes))]
disc.mp.genes <- disc.mp.genes[order(rownames(disc.mp.genes)),]
View(disc.mp.genes)
head(disc.mp.genes)
str(disc.mp.genes)
dim(disc.mp.genes)
print("Step 1 completed: Loading Metabolic Gene expression in Discovery set")

#---2.Available METABRIC Clinical File as sampleAnnotation.disc =colData----------------------
#Metabric Clinical file is uploaded in R.Codes Repository
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

# Print head of the colAnnotation for verification
head(Metabric_Manual_disc)

###Up to codes Above is the same for disc.tf.genes ###

print("Step 2 completed: Manual METABRIC Clinical data")

#---3.disc_mp Sample Annotation-----------------------------
disc_mp_sample_ids <- colnames(disc.mp.genes)
View(disc_mp_sample_ids)
length(disc_mp_sample_ids)

sort(disc_mp_sample_ids)
View(disc_mp_sample_ids)
length(disc_mp_sample_ids)

# Find the common sample IDs between the two datasets
common_samples_disc <- intersect(rownames(Metabric_Manual_disc), disc_mp_sample_ids)
common_samples_disc <- sort(common_samples_disc)
View(common_samples_disc)
length(common_samples_disc)

# Subset and reorder both datasets to only include common samples
counts.mp.disc <- disc.mp.genes[, common_samples_disc]
all.equal(rownames(disc.mp.genes), rownames(counts.mp.disc))
View(counts.mp.disc)
dim(counts.mp.disc)

sampleAnnotation.disc <- Metabric_Manual_disc[common_samples_disc, ]
dim(sampleAnnotation.disc)

# Verify that they are aligned
all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc))  # Should return TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.mp.disc), rownames(sample_metadata))
#sampleAnnotation.disc <- sample_metadata[tempindex, ]

#Check again

all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc))

print("Step 3 completed: Discovery SampleAnnotation Prepared")

#---4.Load Ensemble Gene Annotation----------------------------------------------
#1.Accessing the data available in Ensembl by biomaRT
#2.Selecting an Ensembl BioMart database and dataset

# Connect to Ensembl database, and query human genes
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


#2.1##Step1: Identifying the database you need
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
ensembl

#2.2##Step 2: Choosing a dataset
#we look at which datasets are available in the selected BioMart
#by using the function listDatasets()
datasets <- listDatasets(ensembl)

##if we want to find the details of any datasets
##in our ensembl mart that contain the term ‘hsapiens’
##we could do the following
searchDatasets(mart = ensembl, pattern = "hsapiens")

ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)


#2.3##Ensembl mirror sites

ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "asia")

#2.4##Using archived versions of Ensembl
listEnsemblArchives()
listEnsembl(version = 112)
ensembl_112 <- useEnsembl(biomart = "genes", 
                          dataset = "hsapiens_gene_ensembl",
                          host = "https://may2024.archive.ensembl.org",
                          version = 112)

View(listAttributes(ensembl_112))

print("Step 8 completed: Ensembl Gene Information loaded")

#---5.disc_mp Gene Annotation ---------------------------
# Extract gene annotations for the corresponding genes from `ensembl_112`
gene_annot_mp_disc <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name","description"),
  filters = "external_gene_name",
  values = rownames(disc.mp.genes),
  mart = ensembl_112
)

View(gene_annot_mp_disc)
# Sort gene_annot_mp_disc by external_gene_name
gene_annot_mp_disc <- gene_annot_mp_disc[order(gene_annot_mp_disc$external_gene_name), ]

# Check for duplicated external_gene_names
table(duplicated(gene_annot_mp_disc$external_gene_name))

# Remove duplicates based on 'external_gene_name'
gene_annot_mp_disc <- gene_annot_mp_disc[!duplicated(gene_annot_mp_disc$external_gene_name), ]
dim(gene_annot_mp_disc)
length(gene_annot_mp_disc)
View(gene_annot_mp_disc)
dim(counts.mp.disc)

# Save gene_annot_mp_disc as a TSV file
#write.table(gene_annot_mp_disc, "gene_annot_mp_disc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#write.table(gene_annot_mp_disc, "gene_annot_mp_disc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Extract gene IDs from both datasets
ensemble_ids <- gene_annot_mp_disc$external_gene_name
View(ensemble_ids)

ensemble_ids <- sort(ensemble_ids)
length(ensemble_ids)

disc_mp_gene_ids <- rownames(counts.mp.disc)
View(disc_mp_gene_ids)

common_genes <- intersect(ensemble_ids, disc_mp_gene_ids) # Should match or be close to 1420
common_genes <- sort(common_genes)
all.equal(ensemble_ids, common_genes)
length(common_genes)

# Subset and reorder both datasets to only include common samples
counts.mp.disc <- counts.mp.disc[common_genes, ]
View(counts.mp.disc)

all.equal(ensemble_ids,rownames(counts.mp.disc))

gene_annot_mp_disc <- gene_annot_mp_disc[gene_annot_mp_disc$external_gene_name %in% common_genes, ]
rownames(gene_annot_mp_disc) <- gene_annot_mp_disc$external_gene_name
View(gene_annot_mp_disc)
all.equal(gene_annot_mp_disc$external_gene_name, rownames(counts.mp.disc))

#Sort counts.mp.disc and gene_annot_mp_disc based on rownames which is gene names
counts.mp.disc <- counts.mp.disc[order(rownames(counts.mp.disc)), ]
gene_annot_mp_disc <- gene_annot_mp_disc[order(rownames(gene_annot_mp_disc)), ]

all.equal(rownames(counts.mp.disc), rownames(gene_annot_mp_disc))  # Should return TRUE if aligned
dim(counts.mp.disc)
dim(gene_annot_mp_disc)
View(gene_annot_mp_disc)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(gene_annot_mp_disc)
# Suppose you have a data frame named 'df' with a column called 'external_gene_id'
colnames(gene_annot_mp_disc)[colnames(gene_annot_mp_disc) == "external_gene_name"] <- "SYMBOL"
colnames(gene_annot_mp_disc)[colnames(gene_annot_mp_disc) == "ensembl_gene_id"] <- "ENSEMBL"
# Verify the change
colnames(gene_annot_mp_disc)


# One final check:
stopifnot(rownames(gene_annot_mp_disc) == rownames(counts.mp.disc), # features
          rownames(sampleAnnotation.disc) == colnames(counts.mp.disc)) # samples

print("Step 5 completed: GeneAnnotation/RowAnnotation done for discovery set")

#---6.Create disc_mp list/SummarizedExperiment-----------------------
counts.mp.disc <- as.matrix(counts.mp.disc)
View(counts.mp.disc)
str(counts.mp.disc)

disc_mp <- list(
  expData = counts.mp.disc,                       # expData as an assay
  rowAnnotation = gene_annot_mp_disc,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.disc           # Sample annotations (column metadata)
)

View(disc_mp)

#SummarizedExperiment of rtni_disc_mp_se
#rtni_disc_mp_se <- SummarizedExperiment(
#  assays = counts.mp.disc,
#  rowData = gene_annot_mp_disc,
#  colData = sampleAnnotation.disc
#)

#View(rtni_disc_mp_se)

# Check dimensions and alignment
dim(counts.mp.disc)
#dim(disc_mp_se$expData)  # Should match dimensions of counts.mp.disc
all.equal(colnames(disc_mp$expData), rownames(sampleAnnotation.disc))  # Should return TRUE
all.equal(disc_mp$expData, counts.mp.disc)
all.equal(disc_mp$colAnnotation, sampleAnnotation.disc)

print("Step 6 completed: list/SummarizedExperiment file prepared for Discovery")

#---7.RegulatoryElements-----------------
# Load TF annotation
data("tfsData")
# Check TF annotation:
# Intersect TFs from Lambert et al. (2018) with gene annotation 
# from the gene expression of 90 metabolic pathway cohort
regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, disc_mp$rowAnnotation$SYMBOL)
View(regulatoryElements)
regulatoryElements <- sort(regulatoryElements)
View(regulatoryElements)

print("Step 7 completed: RegulatoryElements Defined")

#---8.Run the TNI constructor with the extracted matrix for disc_mp--------------
#This dataset consists of a list with 3 objects:
##a named gene expression matrix (tniData$expData),
##a data frame with gene annotations (tniData$rowAnnotation), 
##and a data frame with sample annotations (tniData$colAnnotation).
##alternatively, 'expData' can be a 'SummarizedExperiment' object
rtni_disc <- tni.constructor(expData = disc_mp$expData, 
                             regulatoryElements = regulatoryElements, 
                             rowAnnotation = disc_mp$rowAnnotation, 
                             colAnnotation = disc_mp$colAnnotation)

###snow clustering returns zero results
### Compute the reference regulatory network by permutation and bootstrap analyses.
### Please set 'spec' according to your available hardware
#options(cluster=snow::makeCluster(spec=4, "SOCK"))  #???

#1.with pValuCutoff =1e-7 it gives 0 values, even in the sumary
#rtni_disc <- tni.permutation(rtni_disc, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#2.with pValuCutoff= 1e-7 and npermutation = 1000 it gives 0 values, even in the sumary
#rtni_disc <- tni.permutation(rtni_disc,nPermutations = 1000, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#3.  nPermutations >= 1000 with running snow goes the same zero values
#4.  nPermutations >= 1000 without snow package, and pValueCutoff works!
rtni_disc <- tni.permutation(rtni_disc, nPermutations = 1000) 

#Unstable interactions are subsequently removed by bootstrap analysis,
##creates a consensus bootstrap network, referred here as refnet (reference network).
rtni_disc <- tni.bootstrap(rtni_disc)

#stopCluster(getOption("cluster"))  # avoid leak memory

# Compute the DPI-filtered regulatory network
rtni_disc <- tni.dpi.filter(rtni_disc, eps = NA)
tni.regulon.summary(rtni_disc)

#detailed information about a specific regulon
tni.regulon.summary(rtni_disc, regulatoryElements = "DNMT1")

# Save the TNI object for subsequent analyses
save(rtni_disc, file = file.path("~/NCA/scripts/", "rtni_disc091124.Edges1761RData"))

print("Step 8 completed: Loaded rtni data and saved.")

regulons <- tni.get(rtni_disc, what = "regulons.and.mode", idkey = "SYMBOL")
View(regulons)
head(regulons)

###To extract regulons all in one txt.###
# Find the maximum number of genes across all regulons
max_genes <- max(sapply(regulon.NA, function(x) if (is.null(x)) 0 else length(x)))

# Create a list of data frames with equal row lengths
regulon_list <- lapply(names(regulon.NA), function(regulon_name) {
  regulon_data <- regulon.NA[[regulon_name]]
  
  # Handle empty or NULL regulons
  if (is.null(regulon_data) || length(regulon_data) == 0) {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  } else if (is.vector(regulon_data)) {
    # Ensure matching lengths for Gene and Value
    genes <- names(regulon_data)
    values <- regulon_data
    if (length(genes) == 0) genes <- rep(NA, length(values))
    df <- data.frame(
      Gene = genes,
      Value = values,
      stringsAsFactors = FALSE
    )
  } else if (is.matrix(regulon_data) || is.data.frame(regulon_data)) {
    df <- data.frame(
      Gene = rownames(regulon_data),
      Value = regulon_data[, 1], # Assuming values are in the first column
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  }
  
  # Extend to max_genes rows if needed
  if (nrow(df) < max_genes) {
    df <- rbind(df, data.frame(
      Gene = rep(NA, max_genes - nrow(df)),
      Value = rep(NA, max_genes - nrow(df))
    ))
  }
  
  # Rename columns with regulon names
  colnames(df) <- c(paste0(regulon_name, "_Gene"), paste0(regulon_name, "_Value"))
  return(df)
})

# Combine all into one data frame
regulon_df <- do.call(cbind, regulon_list)

# Write to file
write.table(regulon_df, file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/", "Disc.regulon.NA.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

#---9. Regulon activity profiles[skipped]----------

library(Fletcher2013b)
library(pheatmap)
library(grid)
library(gridExtra)

# Load 'rtni1st' data object, which includes regulons and expression profiles
#here I will use: rtni_disc

# A list of transcription factors of interest ( #9 )
# Compute regulon activity for individual samples
rtni_disc <- tni.gsea2(rtni_disc, regulatoryElements = rtni_disc@regulatoryElements)
metabric_regact_disc <- tni.get(rtni_disc, what = "regulonActivity")
View(metabric_regact_disc)

# Get sample attributes from the 'rtni_disc' dataset
metabric_annot_disc <- tni.get(rtni_disc, "colAnnotation")

# Get ER+/- and PAM50 attributes for pheatmap
attribs_disc <- c("LumA","LumB","Basal","Her2","Normal","ER+","ER-")
metabric_annot_disc <- metabric_annot_disc[,attribs_disc]

save(rtni_disc, file = file.path("~/NCA/scripts/", "rtni_disc091124.RData"))


# Define custom colors for each category in annotation_col
annotation_colors <- list(
  "ER+" = c("0" = "lightgray", "1" = "blue"),
  "ER-" = c("0" = "lightgray", "1" = "red"),
  "Normal" = c("0" = "lightgray", "1" = "green"),
  "Her2" = c("0"="lightgrey", "1"= "yellow"),
  "Basal"= c("0"="lightgrey","1"="darkorange"),
  "LumB" = c("0" = "lightgrey", "1"="pink"),
  "LumA"= c("0"="lightgrey", "1"="skyblue")
)

pdf("~/NCA/output/Disc.NCA.mp.091124", width = 10, height = 10)

# Plot regulon activity profiles
disc.heatmap <- pheatmap(t(metabric_regact_disc$differential), 
                         main="Discovery Set (n=988 samples)",
                         row_title = "Regulons",
                         row_title_side = "right",
                         annotation_col = metabric_annot_disc, 
                         show_colnames = FALSE, annotation_legend = FALSE, 
                         clustering_method = "ward.D2", fontsize_row = 6,
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         legend = TRUE,
                         annotation_colors = annotation_colors,
                         border_color = NA,
)

grid.text("Regulons", x= 0.97 , y=0.3, rot=270)

dev.off()
sessionInfo()
#---10.Transcriptional Network Analysis (TNA)-----------
 #---10.1--MP genes-----------

#Since the samples are the same sampleAnnotation.disc_tna will be used
Disc.MP.Genes <- read.delim("~/NCA.ER/data/DISC.MP.Genes.txt",header = TRUE, row.names = 1)
View(Disc.MP.Genes)

#Order by rownames and colnames
Disc.MP.Genes <- Disc.MP.Genes[order(rownames(Disc.MP.Genes)),order(colnames(Disc.MP.Genes))]

# Replace dots with underscores in column names if needed
colnames(Disc.MP.Genes) <- gsub("\\.", "_", colnames(Disc.MP.Genes))

View(Disc.MP.Genes)
dim(Disc.MP.Genes) #[1] 1420  993

counts.mp.disc <- as.matrix(disc.mp[,common_samples_disc])
identical(colnames(counts.mp.disc),colnames(counts.mp.tf.disc))
str(counts.mp.disc)
#num [1:1420, 1:988] 6.24 5.23 6.78 5.58 5.71 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:1420] "A4GALT" "A4GNT" "AACS" "AADAC" ...
#..$ : chr [1:988] "MB_0005" "MB_0006" "MB_0008" "MB_0014" ...

sampleAnnotation.disc_tna_mp <- METABRIC_Manual_Disc[common_samples_disc, ]
dim(sampleAnnotation.disc_tna_mp) #[1] 988  36
View(sampleAnnotation.disc_tna_mp)


# Verify that they are aligned
all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc_tna_mp))  # Should return TRUE

print("Step 10.3 completed: Discovery SampleAnnotation Prepared")

# Check that row names of sample_annotations match column names of expression_data
if (!all(rownames(sampleAnnotation.disc_tna_mp) == colnames(counts.mp.disc))) {
  stop("Mismatch between sample annotation rownames and expression data colnames!")
}

# Load required package
library(limma)

# Create the design matrix for the linear model
# Assuming your label column is named "ER_status" with values "ERpos" and "ERneg"
design_mp_pos <- model.matrix(~ 0 + factor(sampleAnnotation.disc_tna_mp$ER.Status))
colnames(design_mp_pos) <- levels(factor(sampleAnnotation.disc_tna_mp$ER.Status))
rownames(design_mp_pos) <- rownames(sampleAnnotation.disc_tna_mp)
View(design_mp_pos)

all.equal(as.vector(design_mp_pos[,"Positive"]),as.vector(sampleAnnotation.disc$`ER+`))
all(design_mp_pos[,"Positive"] == sampleAnnotation.disc$`ER+`)
all(design_mp_pos[,"Negative"] == sampleAnnotation.disc$`ER-`)

# Fit the linear model using limma
fit_mp_pos <- lmFit(counts.mp.disc, design_mp_pos)
View(fit_mp_pos)

# Create contrast matrix to compare ERpos vs ERneg
contrast_matrix_mp_pos <- makeContrasts(Positive_vs_Negative = Positive - Negative, levels = design_mp_pos)

# Apply the contrast matrix
fit2_mp_pos <- contrasts.fit(fit_mp_pos, contrast_matrix_mp_pos)

# Empirical Bayes adjustment
fit2_mp_pos <- eBayes(fit2_mp_pos)

# Extract results (log2 fold changes, p-values, etc.)
phenotype_mp_pos <- topTable(fit2_mp_pos, coef = "Positive_vs_Negative", adjust.method = "BH", number = Inf)

# Order phenotype alphabetically by row names
phenotype_mp_pos <- phenotype_mp_pos[order(rownames(phenotype_mp_pos)), ]

# Save results to a file if needed
write.table(phenotype_mp_pos,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/","DEG.Disc.MP.txt"),sep = "\t")

# Output the top results for inspection
View(phenotype_mp_pos)
dim(phenotype_mp_pos)  #[1] 1420    6

all(rownames(phenotype_mp_pos) == rownames(counts.mp.disc))

# Filter genes with significant adjusted p-value [& abs(logFC) > 1]
hits_disc.mp <- subset(phenotype_mp_pos, adj.P.Val < 0.05)

# View the differentially expressed genes
View(hits_disc.mp)
dim(hits_disc.mp) #[1] 976   6

#Error: NOTE: all names in 'phenotype' should be available in col1 of 'phenoIDs'!
library(dplyr)
gene_annot_mp_disc_tna <- gene_annot_mp.tf_disc %>%
  select(SYMBOL, everything())

gene_annot_mp_disc_tna_common <- intersect(gene_annot_mp_disc_tna$SYMBOL,rownames(hits_disc.mp))
gene_annot_mp_disc_tna <- gene_annot_mp_disc_tna[gene_annot_mp_disc_tna_common,]
View(gene_annot_mp_disc_tna)


all(gene_annot_mp_disc_tna$SYMBOL == gene_annot_mp_disc_tna_common)

hits_disc.mp <- hits_disc.mp[gene_annot_mp_disc_tna_common,]

# Extract 'logFC' as a named numeric vector
logFC_disc_mp <- setNames(hits_disc.mp$logFC, rownames(hits_disc.mp))
View(logFC_disc_mp)

all.equal(names(logFC_disc_mp),rownames(hits_disc.mp))
all.equal(as.vector(logFC_disc_mp),hits_disc.mp$logFC)

tna.disc_mp <- list(
  phenotype = logFC_disc_mp,
  phenoID = gene_annot_mp_disc_tna,
  hits = rownames(hits_disc.mp)
)

View(tna.disc_mp)
save(phenotype_mp_pos,logFC_disc_mp,hits_disc.mp,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/","Differential.Exp.Disc.MP.RData"))

# Input 1: 'object', a TNI object with regulons
# Input 2: 'phenotype', a named numeric vector, usually log2 differential expression levels
# Input 3: 'hits', a character vector, usually a set of differentially expressed genes
# Input 4: 'phenoIDs', an optional data frame with gene anottation mapped to the phenotype

#set.seed( Since I didnt run everything from scratch for this part I used set.seed() in saved Rproject

#rtni_disc_mp.tf.NA after dpi.filter
#CHECK "rtnaData"                       
rtna_disc <- tni2tna.preprocess(object = rtni_disc, 
                                phenotype = disc_tf$expData, 
                                hits = tnaData$hits, 
                                phenoIDs = tnaData$phenoIDs)
#######
#20k genes in discovery and validation datasets 
#-Preprocessing for input data...
#--Mapping 'phenotype' to 'phenoIDs'...
#--Mapping 'hits' to 'phenoIDs'...
#--Mapping 'transcriptionalNetwork' annotation to 'phenotype'...
#--Checking agreement between 'transcriptionalNetwork' and 'phenotype'... 5.3% ! 
# Error: NOTE: 94.7% of 'transcriptionalNetwork' targets not represented in the 'phenotype'!

#######
#First MP+TF gene expression to build the rtni, with an argument eps = NA in tni.dpi.filter,
#then using limma R Package the differential expression of metabolic genes were calculated 
#and genes with adj.pvalue<0.05 were selected to build the rtna; THE RESULTS WHILE BUILDING THE RTNA WAS AS BELOW:

#-Preprocessing for input data...
#--Mapping 'phenotype' to 'phenoIDs'...
#--Mapping 'hits' to 'phenoIDs'...
#-Mapping 'transcriptionalNetwork' annotation to 'phenotype'...
#--Checking agreement between 'transcriptionalNetwork' and 'phenotype'... 33.4% ! 
#  --Extracting regulons...
#-Preprocessing complete!
#  
#  Warning message:
#  NOTE: 66.6% of 'transcriptionalNetwork' targets not represented in the 'phenotype'! 

#First; all genes in discovery dataset Second 2611 metabolic genes(combined version) for builiding the TNA  
# load("~/NCA.ER/NCA.NEW.MP/DISC.NCA.NEW.MP/5.Disc.New.mp.5tf.RData") to ~/NCA.ER/NCA.Disc.All.Genes/NCA.Disc.All.Gene.Rproject/.RData]                       
#rtna_discovery <- tni2tna.preprocess(object = rtni_Disc_all.NA, 
#+                                    phenotype = tna.disc_mp$phenotype, 
#+                                    hits = tna.disc_mp$hits, 
#+                                    phenoIDs = tna.disc_mp$phenoID)
#-Preprocessing for input data...
#--Mapping 'phenotype' to 'phenoIDs'...
#--Mapping 'hits' to 'phenoIDs'...
#--Mapping 'transcriptionalNetwork' annotation to 'phenotype'...
#--Checking agreement between 'transcriptionalNetwork' and 'phenotype'... 9.3% ! 
#Error: NOTE: 90.7% of 'transcriptionalNetwork' targets not represented in the 'phenotype'!
# save.image("~/NCA.ER/NCA.Disc.All.Genes/All.Then.MP.RData")                    
                        
# Run the MRA method
rtna_disc.mp <- tna.mra(rtna_disc.mp)
#-Performing master regulatory analysis...
#--For 643 regulons...
#|==============================================================================| 100%
#Master regulatory analysis complete

# Get MRA results;
#..setting 'ntop = -1' will return all results, regardless of a threshold
mra_disc.mp <- tna.get(rtna_disc.mp, what="mra", ntop = -1)
View(mra_disc.mp)
write.table(mra_disc.mp,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/","DISC.MRA.MP.txt"),sep = "\t")

# Run the GSEA method
# Please set nPermutations >= 1000
rtna_disc.mp <- tna.gsea1(rtna_disc.mp, nPermutations=1000)
#-Performing gene set enrichment analysis...
#--For 438 regulons...
#|==============================================================================| 100%
#-Gene set enrichment analysis complete 

rtna_disc.mp.15 <- tna.gsea1(rtna_disc.mp, nPermutations=1000,minRegulonSize = 15)
#-Performing gene set enrichment analysis...
#--For 438 regulons...
#|==============================================================================| 100%
#-Gene set enrichment analysis complete 

rtna_disc.mp.15.filtermethod <- tna.gsea1(rtna_disc.mp, nPermutations=1000,sizeFilterMethod = "posANDneg",minRegulonSize = 15)
#-Performing gene set enrichment analysis...
#--For 12 regulons...
#|==============================================================================| 100%
#-Gene set enrichment analysis complete 

# Get GSEA results
gsea1_disc.mp <- tna.get(rtna_disc.mp, what="gsea1", ntop = -1)
head(gsea1_disc.mp)
#Regulon Regulon.Size Observed.Score     Pvalue Adjusted.Pvalue
#ENSG00000065978    YBX1           72           0.67 1.0634e-08      7.7625e-07
#ENSG00000091831    ESR1           56           0.72 1.4048e-08      8.7900e-07
#ENSG00000107485   GATA3           48           0.73 1.7182e-07      8.3618e-06
#ENSG00000119866  BCL11A           47           0.69 1.7526e-06      6.9787e-05
#ENSG00000173894    CBX2           58           0.64 1.0092e-05      3.6836e-04
#ENSG00000125850   OVOL2           32           0.74 2.1966e-05      7.4008e-04

View(gsea1_disc.mp)

gsea1_disc.mp.15.filtered <- tna.get(rtna_disc.mp.15.filtermethod, what="gsea1", ntop = -1)
head(gsea1_disc.mp.15.filtered)
#Regulon Regulon.Size Observed.Score     Pvalue Adjusted.Pvalue
#ENSG00000065978    YBX1           72           0.67 1.6002e-08      1.0013e-06
#ENSG00000091831    ESR1           56           0.72 1.9495e-08      1.0674e-06
#NSG00000107485   GATA3           48           0.73 1.9519e-07      8.5492e-06
#ENSG00000119866  BCL11A           47           0.69 8.5550e-07      3.4064e-05
#ENSG00000125850   OVOL2           32           0.74 6.5153e-06      2.3781e-04
#ENSG00000173894    CBX2           58           0.64 1.2187e-05      4.0821e-04

write.table(gsea1_disc.mp,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/","DISC.tna.gsea1.mp.txt"),sep = "\t",row.names = TRUE, col.names = TRUE)

# Filter for significant TFs
gsea1_disc.mp.sig <- gsea1_disc.mp[gsea1_disc.mp$Adjusted.Pvalue <= 0.05, ]
View(gsea1_disc.mp.sig)
write.table(gsea1_disc.mp.sig,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/","DISC.sig.tna.gsea1.mp.txt"),sep = "\t")

# Plot GSEA results for significant TFs
tna.plot.gsea1(
  rtna_disc.mp.15.filtermethod,
  ntop = 5,
  labPheno = "abs(log2 fold changes)",
  tfs = rownames(gsea1_disc.mp.15.filtered), # Include only significant TFs
  file = "gsea1_filtered_disc_mp_top5.pdf",
  filepath = "~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/"
)

# Plot GSEA results
tna.plot.gsea1(rtna_disc.mp,
               labPheno="abs(log2 fold changes),Disc", 
               ntop = 5,
               file = "tna.gsea1_no_top5.mp.disc", 
               filepath = "~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/",
               ylimPanels = c(0.0,3.5,0.0,2)
)


# Run the GSEA-2T method
# Please set nPermutations >= 1000
rtna_disc.mp <- tna.gsea2(rtna_disc.mp, nPermutations = 1000)
#-Performing two-tailed GSEA analysis...
#--For 438 regulons...
#|==============================================================================| 100%
#|==============================================================================| 100%
#GSEA2 analysis complete 
View(rtna_disc.mp)

# Get GSEA-2T results
gsea2_disc.mp <- tna.get(rtna_disc.mp, what = "gsea2", ntop = -1)
head(gsea2_disc.mp$differential)
#Regulon Regulon.Size Observed.Score     Pvalue Adjusted.Pvalue
#ENSG00000065978    YBX1           72          -0.84 0.00031238       0.0040515
#ENSG00000153207  AHCTF1           18          -1.34 0.00099900       0.0040515
#ENSG00000131668   BARX1           15          -1.91 0.00099900       0.0040515
#ENSG00000123685   BATF3           14          -0.79 0.00099900       0.0040515
#ENSG00000119866  BCL11A           47          -1.46 0.00099900       0.0040515
#ENSG00000134107 BHLHE40           14           0.78 0.00099900       0.0040515

write.table(gsea2_disc.mp,file = file.path("~/NCA.ER/NCA.MP.TF.Finalized/NCA.Disc.MP.tf/tni.dpi.filter.epsNA/","DISC.tna.gsea2.mp.txt"),sep = "\t",row.names = TRUE, col.names = TRUE)

# Plot GSEA-2T results
tna.plot.gsea2(rtna_disc.mp, labPheno="log2 fold changes,Disc", tfs="YBX1", file = "Disc.tna.gsea2_mp")

