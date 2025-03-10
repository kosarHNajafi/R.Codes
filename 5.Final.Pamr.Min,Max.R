#############################
### Section 0: Setup & Libraries
#############################
rm(list = ls(all = TRUE))
# Optionally, load a previously saved RData if needed:
load("E:/5.Ferroptosis/Pam.r/Pathifier/Disc.Pathifier.pamr/2.Disc-pathifier-2025-03-04.RData")

# Define cohort, method, and file prefix (num)
cohort_name <- "Disc"
method <- "pathifier"
num <- "th=4"

# Load requNired libraries
library(pamr)

#############################
### Section 1: Load & Inspect Training Data
#############################

# Load training expression data 
data1 <- as.matrix(read.table("E:/5.Ferroptosis/Pam.r/data/Disc.Training.txt",
                              check.names = FALSE, row.names = 1, header = TRUE, sep = "\t"))
# Order data (by genes and samples) to ensure consistency
data1 <- data1[order(rownames(data1)), order(colnames(data1))]
str(data1)

# Compute gene-wise min and max from training data
gene_min <- apply(data1, 1, min)
gene_max <- apply(data1, 1, max)

# Load tumor_results_df and ensure ordering
tumor_results_df <- tumor_results_df[order(rownames(tumor_results_df)), order(colnames(tumor_results_df))]
# Check that sample IDs in tumor_results_df match the columns of data1
if (!all.equal(tumor_results_df$SampleID, colnames(data1))) {
  stop("Sample IDs do not match between tumor_results_df and data1!")
}
# Extract cluster labels
vecto <- tumor_results_df$Cluster
str(vecto)  # Confirm it is a character vector of proper length

# Start a single PDF device to store all plots
pdf(file = paste0(num,"-",cohort_name,"-",method, Sys.Date(), ".pdf"))

#############################
### Section 2: Train PAM Model
#############################

# Prepare data list for PAM: features (x) and labels (y)
mydata <- list(x = data1, y = vecto)
mydata$genenames <- rownames(mydata$x)
mydata$geneid <- rownames(mydata$x)

# Train PAM model using the raw (ranged) data
mytrain <- pamr.train(mydata)
# Cross-validation to evaluate model performance (both versions)
mycv  <- pamr.cv(mytrain, mydata)
# Plot CV curves and centroid plots for both model versions
pamr.plotcv(mycv)
# Plot centroid plots for selected thresholds
thresholds_to_plot <- c(0, 0.2, 0.5, 0.69, 1,1.5, 2,2.5, 3,3.5, 4)
for (th in thresholds_to_plot) {
  pamr.plotcen(mytrain, mydata, threshold = th)
  title(paste("Value of threshold =", th),cex.main = 1, font.main = 1)
}

#############################
### Section 3: Prediction Phase & Diagnostics
#############################

# Define the datasets for prediction
datasets <- c("Disc","Valid","Combined","TCGA", "Cell")
epsilon <- 1e-6  # To avoid division by zero if needed

# Loop through each prediction dataset
for (cohort in datasets) {
  
  cat("Predicting for cohort:", cohort, "\n")
  
  # Load the new dataset
  new_data_path <- paste0("E:/5.Ferroptosis/Pam.r/data/", cohort, ".Training.txt")
  newx <- as.matrix(read.table(new_data_path, sep = "\t", check.names = FALSE, row.names = 1, header = TRUE))
  newx <- newx[order(rownames(newx)), order(colnames(newx))]
  # Identify common genes between training data and new dataset
  common_gene <- intersect(rownames(data1), rownames(newx))
  if (length(common_gene) == 0) {
    stop("No common genes found between training data and cohort: ", cohort)
  }
  
  # Report any missing genes
  missing_genes <- setdiff(rownames(data1), rownames(newx))
  if (length(missing_genes) > 0) {
    cat("Missing genes in new data for cohort", cohort, ": ", paste(missing_genes, collapse = ", "), "\n")
  }
  
  # Subset the new data to the common genes
  newx_common <- newx[common_gene, , drop = FALSE]
  
  # --- Transform new data so its range matches the training data's range (gene-wise) ---
  newx_common_adj <- newx_common  # Copy the data
  for(i in 1:nrow(newx_common_adj)) {
    train_min <- gene_min[ common_gene[i] ]
    train_max <- gene_max[ common_gene[i] ]
    new_gene_vals <- newx_common_adj[i, ]
    new_min <- min(new_gene_vals)
    new_max <- max(new_gene_vals)
    
    if(new_max - new_min < epsilon) {
      newx_common_adj[i, ] <- train_min
    } else {
      newx_common_adj[i, ] <- (new_gene_vals - new_min) / (new_max - new_min) * (train_max - train_min) + train_min
    }
  }
  assign(paste0("newx_common_adj_", cohort), newx_common_adj, envir = .GlobalEnv)
  write.table(newx_common_adj,file = paste0(num,"-",cohort_name,"-",cohort,"-",method,"_newx_common_adj.txt"),sep = "\t")
  # Use a chosen threshold (here, using the last value of th from the loop) for prediction
  prediction_result <- pamr.predict(mytrain, newx_common_adj, threshold = 4)
  
  # Create a dataframe with prediction results
  prediction_df <- data.frame(SampleID = colnames(newx_common_adj),
                              PredictedLabel = as.character(prediction_result))
  assign(paste0("prediction_df_", cohort), prediction_df, envir = .GlobalEnv)
  
  # Diagnostic 5: Check predicted label proportions
  cat("Prediction summary for", cohort, ":\n")
  print(table(prediction_df$PredictedLabel))
  write.table(prediction_df,file = paste0(num,"-",cohort_name,"-",cohort,"-",method,"-result.txt"),sep = "\t",row.names = FALSE)
}  # End of prediction loop

#############################
### Section 4: Save Environment & Final Diagnostics
#############################

# Close the PDF device so that all plots are saved in one file
dev.off()

# Optionally, save the R session state for later inspection
save(list = ls(), file = paste0(num,"-",cohort_name,"-",cohort,"-",method, "-prediction_", Sys.Date(), ".RData"))

# End of script.
