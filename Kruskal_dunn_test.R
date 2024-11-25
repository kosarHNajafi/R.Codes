#DISCOVERY <- read.delim("E:/Ferroptosis/DISCOVERY.txt")
#View(DISCOVERY)
rownames(DISCOVERY) <- DISCOVERY[,2]
disc <- DISCOVERY[,-c(1,2)]
View(disc)
# Load required libraries
library(dunn.test)
library(dplyr)
library(purrr)

# Function to run Dunn's test with Sidak adjustment and ties correction
run_dunn_test <- function(disc, gene) {
  # Select data for the gene
  selected_data <- disc %>% 
    select(Label = 1, Gene = !!sym(gene)) %>% 
    filter(!is.na(Gene))
  
  # Check if there are enough samples in each group
  if (nrow(selected_data) < 3 || length(unique(selected_data$Label)) < 2) {
    message(paste("Insufficient data for gene:", gene))
    return(data.frame(
      Gene = gene, kruskal_chi_squared = NA, kruskal_p_value = NA,
      dunn_HF_vs_IF_z = NA, dunn_HF_vs_LF_z = NA, dunn_IF_vs_LF_z = NA,
      dunn_HF_vs_IF_p = NA, dunn_HF_vs_LF_p = NA, dunn_IF_vs_LF_p = NA,
      dunn_HF_vs_IF_adj_p = NA, dunn_HF_vs_LF_adj_p = NA, dunn_IF_vs_LF_adj_p = NA,
      stringsAsFactors = FALSE))
  }
  
  # Run Kruskal-Wallis test
  kw_result <- kruskal.test(Gene ~ Label, data = selected_data)
  kruskal_chi_squared <- kw_result$statistic
  kruskal_p_value <- kw_result$p.value
  
  # Run Dunn's test with Sidak correction
  dunn_result <- tryCatch({
    dunn.test(selected_data$Gene, selected_data$Label, method = "bonferroni", kw = TRUE)
  }, error = function(e) {
    message(paste("Error with gene:", gene, "-", e$message))
    return(NULL)
  })
  
  # Handle Dunn's test failure
  if (is.null(dunn_result) || length(dunn_result$comparisons) == 0) {
    message(paste("Dunn's test returned no results for gene:", gene))
    return(data.frame(
      Gene = gene, kruskal_chi_squared = NA, kruskal_p_value = NA,
      dunn_HF_vs_IF_z = NA, dunn_HF_vs_LF_z = NA, dunn_IF_vs_LF_z = NA,
      dunn_HF_vs_IF_p = NA, dunn_HF_vs_LF_p = NA, dunn_IF_vs_LF_p = NA,
      dunn_HF_vs_IF_adj_p = NA, dunn_HF_vs_LF_adj_p = NA, dunn_IF_vs_LF_adj_p = NA,
      stringsAsFactors = FALSE))
  }
  
  # Extract z-scores and p-values
  comparisons <- dunn_result$comparisons
  z_scores <- abs(dunn_result$Z)
  p_values <- dunn_result$P.adjusted
  unadjusted_p_values <- dunn_result$P
  
  # Initialize variables for the comparisons
  dunn_HF_vs_IF_z <- dunn_HF_vs_LF_z <- dunn_IF_vs_LF_z <- NA
  dunn_HF_vs_IF_p <- dunn_HF_vs_LF_p <- dunn_IF_vs_LF_p <- NA
  dunn_HF_vs_IF_adj_p <- dunn_HF_vs_LF_adj_p <- dunn_IF_vs_LF_adj_p <- NA
  
  # Assign values based on comparisons
  for (i in seq_along(comparisons)) {
    if (comparisons[i] == "HF - IF") {
      dunn_HF_vs_IF_z <- z_scores[i]
      dunn_HF_vs_IF_p <- unadjusted_p_values[i]  # Unadjusted p-value
      dunn_HF_vs_IF_adj_p <- p_values[i]   # Adjusted p-value
    } else if (comparisons[i] == "HF - LF") {
      dunn_HF_vs_LF_z <- z_scores[i]
      dunn_HF_vs_LF_p <- unadjusted_p_values[i]
      dunn_HF_vs_LF_adj_p <- p_values[i]
    } else if (comparisons[i] == "IF - LF") {
      dunn_IF_vs_LF_z <- z_scores[i]
      dunn_IF_vs_LF_p <- unadjusted_p_values[i]
      dunn_IF_vs_LF_adj_p <- p_values[i]
    }
  }
  
  # Return results in a data frame
  return(data.frame(
    Gene = gene, kruskal_chi_squared = kruskal_chi_squared, kruskal_p_value = kruskal_p_value,
    dunn_HF_vs_IF_z = dunn_HF_vs_IF_z, dunn_HF_vs_LF_z = dunn_HF_vs_LF_z, dunn_IF_vs_LF_z = dunn_IF_vs_LF_z,
    dunn_HF_vs_IF_p = dunn_HF_vs_IF_p, dunn_HF_vs_LF_p = dunn_HF_vs_LF_p, dunn_IF_vs_LF_p = dunn_IF_vs_LF_p,
    dunn_HF_vs_IF_adj_p = dunn_HF_vs_IF_adj_p, dunn_HF_vs_LF_adj_p = dunn_HF_vs_LF_adj_p, dunn_IF_vs_LF_adj_p = dunn_IF_vs_LF_adj_p,
    stringsAsFactors = FALSE))
}

# Iterate over all genes and run Dunn's test
results <- purrr::map_dfr(colnames(disc)[-1], function(gene) {
  message(paste("Processing gene:", gene))
  run_dunn_test(disc, gene)
})

# Remove duplicates: keep only the first occurrence
results <- results %>% distinct(Gene, .keep_all = TRUE)

# Save results to a file
write.table(results, "disc_bonferroni_dunn_test_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

