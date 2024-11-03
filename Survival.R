#DISC.ERP.Survival <- read.delim("D:/2.ER(KM_R)/DISC.ERP.Survival.txt")
#View(DISC.ERP.Survival)

#any(DISC.ERP.Survival$ER.Status == "negative")

#rownames(DISC.ERP.Survival) <- DISC.ERP.Survival[,1]
#disc.erp <- DISC.ERP.Survival[,-c(1,2)]
# Load necessary libraries
library(survival)
library(survminer)
library(cowplot) # For ggdraw

# set the correct file paths for your data files
# DISC.ERP.Survival <- read.delim("F:/2.ER(KM_R)/DISC.ERP.Survival.txt")

# Open a PDF file to save all Kaplan-Meier plots
pdf("ER.Data/Disc.ERpos.pdf", width = 10, height = 8) # Use a valid file path

# Create a flag to track if any plots are created
plot_created <- FALSE

# Loop over each PDS column (from column 7 to 16)
for (i in 7:16) {
  # Get the name of the current PDS and the data for that PDS column
  PDS_name <- colnames(disc.erp)[i]
  PDS <- disc.erp[[i]]
  
  # Ensure we skip columns with all NA values
  if (all(is.na(PDS))) {
    cat("All values are NA for PDS:", PDS_name, "\n")
    next
  }
  
  # Temporarily omit NA values for the purpose of calculating the average
  PDS_clean <- na.omit(PDS)
  
  # Calculate the average of the current PDS (ignoring NAs)
  PDS_avg <- mean(PDS_clean, na.rm = TRUE)
  cat("Processing PDS:", PDS_name, "with average:", PDS_avg, "\n")
  
  # Create two groups based on whether the PDS is above or below the average
  disc.erp$PDS_group <- ifelse(PDS > PDS_avg,
                               paste("High", PDS_name),
                               paste("Low", PDS_name))
  
  # ======================================
  # Kaplan-Meier Curve for Overall Survival (OS)
  OS_surv <- Surv(time = as.numeric(disc.erp$OS..Months.), event = as.numeric(disc.erp$OS.Status))
  OS_fit <- survfit(OS_surv ~ PDS_group, data = disc.erp, na.action = na.exclude)
  
  # Plot the Kaplan-Meier curve for OS
  plot_os <- ggsurvplot(OS_fit, data = disc.erp, pval = TRUE,
                        title = paste("Discovery.ERpos OS -", PDS_name),
                        risk.table = TRUE, conf.int = FALSE,
                        ylim = c(0, 1), # Set y-axis limits
                        ggtheme = theme_minimal() + theme(panel.grid.minor = element_blank()), # Clean theme
                        ylab = "Survival Probability",
                        pval.coord = c(0.8, 0.8)) # Adjust position for p-value
  
  # Ensure that plot limits are valid before calculating max
  if (!is.null(plot_os$plot$limits$y) && length(plot_os$plot$limits$y) > 0) {
    y_pos_os <- max(plot_os$plot$limits$y, na.rm = TRUE) * 0.75
  } else {
    y_pos_os <- 0.75 # Default position if limits are not valid
  }
  
  print(ggdraw(plot_os$plot) + theme(plot.margin = margin(0, 0, 0, 0)))
  plot_created <- TRUE
  
  # ======================================
  # Kaplan-Meier Curve for Disease-Specific Survival (DSS)
  DSS_surv <- Surv(time = as.numeric(disc.erp$OS..Months.), event = as.numeric(disc.erp$DSS))
  DSS_fit <- survfit(DSS_surv ~ PDS_group, data = disc.erp, na.action = na.exclude)
  
  plot_dss <- ggsurvplot(DSS_fit, data = disc.erp, pval = TRUE,
                         title = paste("Discovery.ERpos DSS -", PDS_name),
                         risk.table = TRUE, conf.int = FALSE,
                         ylim = c(0, 1), # Set y-axis limits
                         ggtheme = theme_minimal() + theme(panel.grid.minor = element_blank()), # Clean theme
                         ylab = "Survival Probability",
                         pval.coord = c(0.8, 0.8)) # Adjust position for p-value
  
  # Ensure that plot limits are valid before calculating max
  if (!is.null(plot_dss$plot$limits$y) && length(plot_dss$plot$limits$y) > 0) {
    y_pos_dss <- max(plot_dss$plot$limits$y, na.rm = TRUE) * 0.75
  } else {
    y_pos_dss <- 0.75 # Default position if limits are not valid
  }
  
  print(ggdraw(plot_dss$plot) + theme(plot.margin = margin(0, 0, 0, 0)))
  plot_created <- TRUE
  
  # ======================================
  # Kaplan-Meier Curve for Relapse-Free Survival (RFS)
  RFS_surv <- Surv(time = as.numeric(disc.erp$RFS..Months.), event = as.numeric(disc.erp$RFS))
  RFS_fit <- survfit(RFS_surv ~ PDS_group, data = disc.erp, na.action = na.exclude)
  
  plot_rfs <- ggsurvplot(RFS_fit, data = disc.erp, pval = TRUE,
                         title = paste("Discovery.ERpos RFS -", PDS_name),
                         risk.table = TRUE, conf.int = FALSE,
                         ylim = c(0, 1), # Set y-axis limits
                         ggtheme = theme_minimal() + theme(panel.grid.minor = element_blank()), # Clean theme
                         ylab = "Survival Probability",
                         pval.coord = c(0.8, 0.8)) # Adjust position for p-value
  
  # Ensure that plot limits are valid before calculating max
  if (!is.null(plot_rfs$plot$limits$y) && length(plot_rfs$plot$limits$y) > 0) {
    y_pos_rfs <- max(plot_rfs$plot$limits$y, na.rm = TRUE) * 0.75
  } else {
    y_pos_rfs <- 0.75 # Default position if limits are not valid
  }
  
  print(ggdraw(plot_rfs$plot) + theme(plot.margin = margin(0, 0, 0, 0)))
  plot_created <- TRUE
  
  # ======================================
  # Kaplan-Meier Curve for Metastasis-Free Survival (MFS)
  MFS_surv <- Surv(time = as.numeric(disc.erp$OS..Months.), event = as.numeric(disc.erp$MFS))
  MFS_fit <- survfit(MFS_surv ~ PDS_group, data = disc.erp, na.action = na.exclude)
  
  plot_mfs <- ggsurvplot(MFS_fit, data = disc.erp, pval = TRUE,
                         title = paste("Discovery.ERpos MFS -", PDS_name),
                         risk.table = TRUE, conf.int = FALSE,
                         ylim = c(0, 1), # Set y-axis limits
                         ggtheme = theme_minimal() + theme(panel.grid.minor = element_blank()), # Clean theme
                         ylab = "Survival Probability",
                         pval.coord = c(0.8, 0.8)) # Adjust position for p-value
  
  # Ensure that plot limits are valid before calculating max
  if (!is.null(plot_mfs$plot$limits$y) && length(plot_mfs$plot$limits$y) > 0) {
    y_pos_mfs <- max(plot_mfs$plot$limits$y, na.rm = TRUE) * 0.75
  } else {
    y_pos_mfs <- 0.75 # Default position if limits are not valid
  }
  
  print(ggdraw(plot_mfs$plot) + theme(plot.margin = margin(0, 0, 0, 0)));
  plot_created <- TRUE;
}

# Close the PDF device only if plots were created
if (plot_created) {
  dev.off();
} else {
  cat("No plots were created.\n");
}