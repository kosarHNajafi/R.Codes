#DISC.ERN.Survival <- read.delim("G:/My Drive/ER/Survival Analysis/Pathifier Survival Analysis_R/DISCOVERY/KM/ERN/DISC.ERN.Survival.txt")
#View(DISC.ERN.Survival)

#rownames(DISC.ERN.Survival) <- DISC.ERN.Survival[,1]
#disc.ERneg <- DISC.ERN.Survival[,-c(1,2)]

# Load necessary libraries  
library(survival)  
library(survminer)  
library(cowplot)  # For ggdraw


# Open a PDF file to save all Kaplan-Meier plots  
pdf("ER.Data/DISC.ERNeg.TOP10.(Pretty)Survival.pdf", width = 10, height = 8)  # Use a valid file path  

# Create a flag to track if any plots are created
plot_created <- FALSE

# Loop over each PDS column (from column 7 to 16)  
for (i in 7:16) {  
  
  # Get the name of the current PDS and the data for that PDS column  
  PDS_name <- colnames(disc.ERneg)[i]  
  PDS <- disc.ERneg[[i]]  
  
  # Ensure we skip columns with all NA values  
  if (!all(is.na(PDS))) {  
    
    # Temporarily omit NA values for the purpose of calculating the average  
    PDS_clean <- na.omit(PDS)  
    
    # Calculate the average of the current PDS (ignoring NAs)  
    PDS_avg <- mean(PDS_clean, na.rm = TRUE)  
    
    # Create two groups based on whether the PDS is above or below the average  
    disc.ERneg$PDS_group <- ifelse(PDS > PDS_avg,   
                                 paste("High", PDS_name),   
                                 paste("Low", PDS_name))  
    
    # ======================================  
    # Kaplan-Meier Curve for Overall Survival (OS)  
    OS_surv <- Surv(time = as.numeric(disc.ERneg$OS.Months.), event = as.numeric(disc.ERneg$OS.status))  
    OS_fit <- survfit(OS_surv ~ PDS_group, data = disc.ERneg, na.action = na.exclude)  
    
    # Plot the Kaplan-Meier curve for OS without confidence intervals  
    plot_os <- ggsurvplot(OS_fit, data = disc.ERneg, pval = TRUE,   
                          title = paste("Kaplan-Meier Curve in Disc.ERneg for OS -", PDS_name),  
                          risk.table = TRUE, conf.int = FALSE)  
    
    # Use ggdraw to center the plot
    print(ggdraw(plot_os$plot) + theme(plot.margin = margin(0, 0, 0, 0)))  
    plot_created <- TRUE  
    
    # ======================================  
    # Repeat for DSS, RFS, and MFS (similar code with plot_dss, plot_rfs, plot_mfs)
    
    # Kaplan-Meier Curve for Disease-Specific Survival (DSS)
    DSS_surv <- Surv(time = as.numeric(disc.ERneg$OS.Months.), event = as.numeric(disc.ERneg$DSS))  
    DSS_fit <- survfit(DSS_surv ~ PDS_group, data = disc.ERneg, na.action = na.exclude)  
    
    plot_dss <- ggsurvplot(DSS_fit, data = disc.ERneg, pval = TRUE,   
                           title = paste("Kaplan-Meier Curve for DSS -", PDS_name),  
                           risk.table = TRUE, conf.int = FALSE)  
    
    print(ggdraw(plot_dss$plot) + theme(plot.margin = margin(0, 0, 0, 0)))  
    plot_created <- TRUE
    
    # Kaplan-Meier Curve for Relapse-Free Survival (RFS)
    RFS_surv <- Surv(time = as.numeric(disc.ERneg$RFS..Months.), event = as.numeric(disc.ERneg$RFS))  
    RFS_fit <- survfit(RFS_surv ~ PDS_group, data = disc.ERneg, na.action = na.exclude)  
    
    plot_rfs <- ggsurvplot(RFS_fit, data = disc.ERneg, pval = TRUE,   
                           title = paste("Kaplan-Meier Curve for RFS -", PDS_name),  
                           risk.table = TRUE, conf.int = FALSE)  
    
    print(ggdraw(plot_rfs$plot) + theme(plot.margin = margin(0, 0, 0, 0)))  
    plot_created <- TRUE
    
    # Kaplan-Meier Curve for Metastasis-Free Survival (MFS)
    MFS_surv <- Surv(time = as.numeric(disc.ERneg$OS.Months.), event = as.numeric(disc.ERneg$MFS))  
    MFS_fit <- survfit(MFS_surv ~ PDS_group, data = disc.ERneg, na.action = na.exclude)  
    
    plot_mfs <- ggsurvplot(MFS_fit, data = disc.ERneg, pval = TRUE,   
                           title = paste("Kaplan-Meier Curve for MFS -", PDS_name),  
                           risk.table = TRUE, conf.int = FALSE)  
    
    print(ggdraw(plot_mfs$plot) + theme(plot.margin = margin(0, 0, 0, 0)))  
    plot_created <- TRUE
  }  
}  

# Close the PDF file after generating all the plots  
dev.off()

# If no plots were created, print a message
if (!plot_created) {
  cat("No valid plots were generated.\n")
}
