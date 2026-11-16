# Load packages
library(tools)
library(readxl)
library(dplyr)
library(tidyverse)
library(drc)
library(writexl)

# Define analysis parameters
data_path <- "C:/Users/path-to-data"                          # Replace path with your owb path. Use / not \
pulldown_file <- "DNDI-6690_pulldown_DR_Rep1.csv"             # Include extension (.xlsx/.csv)
pulldown_sheet <- ""                                          # Leave as "" if only one sheet
selected_proteins <- "A0A3Q8IRJ0"                             # Proteins to plot (Protein_IDs/topx_log2FC/topx_IC50)
remove_contaminants <- TRUE                                   # TRUE = yes, FALSE = no
remove_low_peptides <- TRUE                                   # TRUE = yes, FALSE = no
add_imputation <- TRUE                                        # TRUE = yes, FALSE = no
export_IC50s <- TRUE                                          # TRUE = yes, FALSE = no
IC50s_filename <- "LdD1LP_pulldown_DR_Rep1_IC50s.xlsx"        # Include .xlsx extension

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          PART 1: FORMAT DATASET                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load data (supports .xlsx and .csv)
setwd(data_path)
ext <- file_ext(pulldown_file)
pulldown <- switch(
  tolower(ext),
  "xlsx" = if (pulldown_sheet == "") {read_excel(pulldown_file)} else 
    {read_excel(pulldown_file, sheet = pulldown_sheet)},
  "csv"  = read.csv(pulldown_file, stringsAsFactors = FALSE))

# Rename columns
pulldown <- pulldown %>%
  rename(any_of(c(
    Name = "Name",
    Protein_ID = "Protein_ID",
    Gene_ID = "Gene_ID",               # Use "" if no gene information
    Peptides = "Peptides",
    DMSO = "DMSO_L",
    Drug_0.01 = "Comp_0.02_L",   # Use "" if you didn't include this concentration
    Drug_0.1 = "Comp_0.2_L",     # Use "" if you didn't include this concentration
    Drug_1 = "Comp_2_L",         # Use "" if you didn't include this concentration
    Drug_10 = "Comp_20_L",       # Use "" if you didn't include this concentration
    Drug_100 = "Comp_200_L"
  )))

# Clean and prepare data
pulldown <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO", "Drug_0.01", "Drug_0.1", "Drug_1", "Drug_10", "Drug_100", "Intensity", "Kinase", "Host_kinase")))
pulldown <- pulldown %>% dplyr::mutate(across(any_of(c("DMSO", "Drug_0.01", "Drug_0.1", "Drug_1", "Drug_10", "Drug_100")), as.numeric))

# Determine data type
cols <- intersect(c("DMSO", "Drug_0.01", "Drug_0.1", "Drug_1", "Drug_10", "Drug_100"), colnames(pulldown))
if (length(cols) > 0) {
  vals <- unlist(pulldown[cols])
  max_val <- max(vals, na.rm = TRUE)
  data_type <- if (max_val < 1e5) "normalised" else "raw"
}

# Create function to handle NA/NaN values
is.any.na <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  is.na(x) | (is.character(x) & tolower(x) %in% c("na", "nan", "NA", "NaN"))
}

# Remove any previous imputation
remove_imputation <- TRUE  # TRUE = yes, FALSE = no
if (remove_imputation && length(cols) > 0 && all(!is.any.na(pulldown[cols]))) {
  val_remove <- case_when(
    any(vals == 11111, na.rm = TRUE) ~ 11111,
    data_type == "normalised" & max_val > 2000 ~ max_val,
    TRUE ~ min(vals, na.rm = TRUE)
  )
  pulldown <- pulldown %>%
    dplyr::mutate(across(all_of(cols), ~ na_if(., val_remove)))
}

# Calculate max intensity
if (!"Intensity" %in% colnames(pulldown) && length(cols) > 0) {
  pulldown$Intensity <- do.call(pmax, c(pulldown[cols], na.rm = TRUE))
}

# Remove contaminants
if (remove_contaminants) {
  pulldown <- dplyr::filter(pulldown, !grepl("cRAP", Protein_ID))
}

# Remove proteins with low peptide counts (≤ 2)
if (remove_low_peptides) {
  pulldown <- pulldown %>% dplyr::filter(Peptides > 2)
}

# Impute missing values
impute_val <- FALSE
if (add_imputation) {
  if (length(cols) > 0) {
    vals <- unlist(pulldown[cols])
    impute_val <- min(vals, na.rm = TRUE) * 0.999999
    pulldown <- pulldown %>%
      mutate(across(all_of(cols), ~ replace(., is.any.na(.), impute_val)))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                      PART 2: CALCULATE IC50 VALUES                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Calculate log2 FC for highest drug concentration
pulldown$log2FC <- log2(pulldown$DMSO / pulldown$Drug_100)

# Extract drug concentrations
conc <- c(Drug_0.01 = 0.02, Drug_0.1 = 0.2, Drug_1 = 2, Drug_10 = 20, Drug_100 = 200)
conc <- conc[names(conc) %in% colnames(pulldown)] # Keep only concentrations that exist in the dataset
conc <- c(DMSO = min(conc) / 10, conc) # Add DMSO as 10-fold lower than the lowest drug concentration

# Calculate IC50 values
calculate_IC50 <- function(response) {
  fit <- tryCatch({
    model <- drm(response ~ conc, fct = LL.4(fixed = c(NA, NA, NA, NA)))
    ed <- ED(model, 50)
    return(c(IC50 = ed[1], sd = ed[2]))
  }, error = function(e) return(c(NA, NA)))
}

IC50_matrix <- t(apply(pulldown[, cols], 1, calculate_IC50))
pulldown$IC50 <- IC50_matrix[, 1]
pulldown$SD <- IC50_matrix[, 2]

# Save IC50 results to Excel file
if (export_IC50s) {
  write_xlsx(pulldown, IC50s_filename)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                   PART 3: GENERATE DOSE-RESPONSE CURVES                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Define function to get top Protein IDs
get_top_proteins <- function(keyword, data) {
  
  # Extract number and method
  n <- as.numeric(str_extract(keyword, "\\d+"))
  method <- str_extract(keyword, "(log2FC|IC50)")
  
  # Filter and sort direction
  data <- data %>%
    filter(!is.na(.data[[method]])) %>%
    arrange(if (method == "log2FC") desc(.data[[method]]) else .data[[method]])
  
  # Return top n Protein_IDs
  head(data$Protein_ID, n)
}

# Expand protein list
selected_proteins <- unlist(strsplit(selected_proteins, ",\\s*"))
full_protein_list <- unique(unlist(lapply(selected_proteins, function(p) {
  if (grepl("^top\\d+_(log2FC|IC50)$", p)) get_top_proteins(p, pulldown) else p
})))

# Define plot function
plot_dose_response <- function(protein_id, protein_row, response_values, imputed, model) {
  par(mgp = c(6.5, 1, 0), mar = c(8, 8, 4, 2))
  plot(model, type = "n", bp = "",
       main = paste0(protein_id, "\n", protein_row$Name),
       ylab = "Intensity",
       xlab = "",
       xlim = c(0.001, 100),
       cex.main = 1.5,
       cex.sub = 1.5,
       cex.lab = 1.5,
       cex.axis = 1.5,
       lwd = 2.5,
       axes = FALSE
  )
  
  mtext("Drug concentration (μM)", side = 1, line = 3, cex = 1.5)
  mtext(sprintf("IC50 = %.3f ± %.3f μM", protein_row$IC50, protein_row$SD), side = 1, line = 6, cex = 1.5, font = 2)
  
  axis(1, at = conc, labels = sub("Drug_", "", names(conc)), cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  
  points(conc[!imputed], response_values[!imputed], pch = 16, cex = 2)
  points(conc[imputed], response_values[imputed], pch = 18, cex = 2)
}

# Loop through each Protein ID
for (protein_id in full_protein_list) {
  protein_row <- pulldown[pulldown$Protein_ID == protein_id, ]
  if (nrow(protein_row) == 0) next  # Skip if protein not found
  response_values <- as.numeric(protein_row[1, cols])
  imputed <- if (add_imputation) response_values == impute_val else rep(FALSE, length(response_values))
  
  # Fit dose-response model
  model <- drm(response_values ~ conc, fct = LL.4(fixed = c(NA, NA, NA, NA)), na.action = na.omit)
  
  # PNG
  png(file = file.path(data_path, paste0(protein_id, ".png")), width = 550, height = 500)
  plot_dose_response(protein_id, protein_row, response_values, imputed, model)
  dev.off()
    
  # SVG
  svg(file = file.path(data_path, paste0(protein_id, ".svg")), width = 7, height = 6.5)
  plot_dose_response(protein_id, protein_row, response_values, imputed, model)
  dev.off()
}

