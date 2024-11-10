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

#---9. Regulon activity profiles----------

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
rtna_disc <- tni2tna.preprocess(object = rtni_disc, 
                                phenotype = disc_tf$expData, 
                                hits = tnaData$hits, 
                                phenoIDs = tnaData$phenoIDs)