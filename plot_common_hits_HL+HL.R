# Load packages
library(tools)
library(readxl)
library(dplyr)
library(purrr)
library(tidyr)
library(plotly)
library(drc)
library(writexl)

# Define analysis parameters
data_path <- "C:/Users/path-to-data"                              # Replace path with your owb path. Use / not \
pulldown_file_1 <- "DNDI-6690_pulldown_Low_Rep1-enrichment.xlsx"  # Include extension (.xlsx/.csv)
pulldown_sheet_1 <- ""                                            # Leave as "" if only one sheet
pulldown_file_2 <- "DNDI-6690_pulldown_Low_Rep2-enrichment.xlsx"  # Include extension (.xlsx/.csv)
pulldown_sheet_2 <- ""                                            # Leave as "" if only one sheet
load_level <- "L"                                                 # Choose "H" or "L"
highlight_proteins_red <- "A0A3Q8IRJ0"                            # Protein IDs or "(host) kinases" (pre-processed data only)
highlight_proteins_blue <- ""                                     # Protein IDs or "(host) kinases" (pre-processed data only)
highlight_interactors <- FALSE                                    # TRUE = yes, FALSE = no (pre-processed data only)
remove_contaminants <- FALSE                                      # TRUE = yes, FALSE = no
remove_low_peptides <- FALSE                                      # TRUE = yes, FALSE = no
add_imputation <- FALSE                                           # TRUE = yes, FALSE = no
add_hit_threshold <- TRUE                                         # TRUE = yes, FALSE = no
hit_threshold_file_1 <- "median+6MAD"                             # Choose format: "2", "mean+2SD", "median+3MAD"
hit_threshold_file_2 <- "median+6MAD"                             # Choose format: "2", "mean+2SD", "median+3MAD"
export_results <- TRUE                                            # TRUE = yes, FALSE = no
results_filename <- "DNDI-6690_pulldown_Low_common_hits.xlsx"        # Include .xlsx extension

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                          PART 1: FORMAT DATASET                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Load data (supports .xlsx and .csv)
setwd(data_path)
ext_1 <- file_ext(pulldown_file_1)
ext_2 <- file_ext(pulldown_file_2)
pulldown_1 <- switch(
  tolower(ext_1),
  "xlsx" = if (pulldown_sheet_1 == "") {read_excel(pulldown_file_1)} else 
  {read_excel(pulldown_file_1, sheet = pulldown_sheet_1)},
  "csv"  = read.csv(pulldown_file_1, stringsAsFactors = FALSE))
pulldown_2 <- switch(
  tolower(ext_2),
  "xlsx" = if (pulldown_sheet_2 == "") {read_excel(pulldown_file_2)} else 
  {read_excel(pulldown_file_2, sheet = pulldown_sheet_2)},
  "csv"  = read.csv(pulldown_file_2, stringsAsFactors = FALSE))

# Rename columns
pulldown_1 <- pulldown_1 %>%
  rename(any_of(c(
    Name = "Name",
    Protein_ID = "Protein_ID",
    Gene_ID = "Gene_ID",               # Use "" if no gene information
    Peptides = "Peptides",
    Intensity = "Intensity",           # Use "" if no intensity column - it will be calculated automatically
    DMSO_H = "",                       # Use "" if you don't have H data
    Drug_H = "",                       # Use "" if you don't have H data
    DMSO_L = "DMSO_L",                 # Use "" if you don't have L data
    Drug_L = "Comp_L"                  # Use "" if you don't have L data
  )))
pulldown_2 <- pulldown_2 %>%
  rename(any_of(c(
    Name = "Name",
    Protein_ID = "Protein_ID",
    Gene_ID = "Gene_ID",               # Use "" if no gene information
    Peptides = "Peptides",
    Intensity = "Intensity",           # Use "" if no intensity column - it will be calculated automatically
    DMSO_H = "",                       # Use "" if you don't have H data
    Drug_H = "",                       # Use "" if you don't have H data
    DMSO_L = "DMSO_L",                 # Use "" if you don't have L data
    Drug_L = "Comp_L"                  # Use "" if you don't have L data
  )))

# Detect available data
has_H_1 <- all(c("DMSO_H", "Drug_H") %in% colnames(pulldown_1))
has_L_1 <- all(c("DMSO_L", "Drug_L") %in% colnames(pulldown_1))
has_H_2 <- all(c("DMSO_H", "Drug_H") %in% colnames(pulldown_2))
has_L_2 <- all(c("DMSO_L", "Drug_L") %in% colnames(pulldown_2))

# Clean and prepare data
pulldown_1 <- pulldown_1 %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO_H", "Drug_H", "DMSO_L", "Drug_L", "Intensity", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
pulldown_1 <- pulldown_1 %>% dplyr::mutate(across(any_of(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L")), as.numeric))
pulldown_2 <- pulldown_2 %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO_H", "Drug_H", "DMSO_L", "Drug_L", "Intensity", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
pulldown_2 <- pulldown_2 %>% dplyr::mutate(across(any_of(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L")), as.numeric))

# Determine data type
cols_1 <- intersect(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L"), colnames(pulldown_1))
if (length(cols_1) > 0) {
  vals_1 <- unlist(pulldown_1[cols_1])
  max_val_1 <- max(vals_1, na.rm = TRUE)
  data_type_1 <- if (max_val_1 < 1e5) "normalised" else "raw"
}
cols_2 <- intersect(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L"), colnames(pulldown_2))
if (length(cols_2) > 0) {
  vals_2 <- unlist(pulldown_2[cols_2])
  max_val_2 <- max(vals_2, na.rm = TRUE)
  data_type_2 <- if (max_val_2 < 1e5) "normalised" else "raw"
}

# Create function to handle NA/NaN values
is.any.na <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  is.na(x) | (is.character(x) & tolower(x) %in% c("na", "nan", "NA", "NaN"))
}

# Remove any previous imputation
remove_imputation <- TRUE  # TRUE = yes, FALSE = no
if (remove_imputation && length(cols_1) > 0 && all(!is.any.na(pulldown_1[cols_1]))) {
  val_remove_1 <- case_when(
    any(vals_1 == 11111, na.rm = TRUE) ~ 11111,
    data_type_1 == "normalised" & max_val_1 > 2000 ~ max_val_1,
    TRUE ~ min(vals_1, na.rm = TRUE)
  )
  pulldown_1 <- pulldown_1 %>%
    dplyr::mutate(across(all_of(cols_1), ~ na_if(., val_remove_1)))
}
if (remove_imputation && length(cols_2) > 0 && all(!is.any.na(pulldown_2[cols_2]))) {
  val_remove_2 <- case_when(
    any(vals_2 == 11111, na.rm = TRUE) ~ 11111,
    data_type_2 == "normalised" & max_val_2 > 2000 ~ max_val_2,
    TRUE ~ min(vals_2, na.rm = TRUE)
  )
  pulldown_2 <- pulldown_2 %>%
    dplyr::mutate(across(all_of(cols_2), ~ na_if(., val_remove_2)))
}

# Calculate max intensity
if (!"Intensity" %in% colnames(pulldown_1) && length(cols_1) > 0) { pulldown_1$Intensity <- do.call(pmax, c(pulldown_1[cols_1], na.rm = TRUE)) }
if (!"Intensity" %in% colnames(pulldown_2) && length(cols_2) > 0) { pulldown_2$Intensity <- do.call(pmax, c(pulldown_2[cols_2], na.rm = TRUE)) }

# Freeze raw data for hit threshold calculations
raw_pulldown_1 <- pulldown_1
raw_pulldown_2 <- pulldown_2

# Remove contaminants
if (remove_contaminants) {
  pulldown_1 <- dplyr::filter(pulldown_1, !grepl("cRAP", Protein_ID))
  raw_pulldown_1 <- dplyr::filter(raw_pulldown_1, !grepl("cRAP", Protein_ID))
  pulldown_2 <- dplyr::filter(pulldown_2, !grepl("cRAP", Protein_ID))
  raw_pulldown_2 <- dplyr::filter(raw_pulldown_2, !grepl("cRAP", Protein_ID))
}

# Remove proteins with low peptide counts (≤ 2)
if (remove_low_peptides) {
  pulldown_1 <- pulldown_1 %>% dplyr::filter(Peptides > 2)
  raw_pulldown_1 <- raw_pulldown_1 %>% dplyr::filter(Peptides > 2)
  pulldown_2 <- pulldown_2 %>% dplyr::filter(Peptides > 2)
  raw_pulldown_2 <- raw_pulldown_2 %>% dplyr::filter(Peptides > 2)
}

# Impute missing values
impute_val_1 <- FALSE
impute_val_2 <- FALSE
if (add_imputation) {
  if (length(cols_1) > 0) {
    vals_1 <- unlist(pulldown_1[cols_1])
    impute_val_1 <- min(vals_1, na.rm = TRUE) * 0.999999
    pulldown_1 <- pulldown_1 %>%
      mutate(across(all_of(cols_1), ~ replace(., is.any.na(.), impute_val_1)))
  }
  if (length(cols_2) > 0) {
    vals_2 <- unlist(pulldown_2[cols_2])
    impute_val_2 <- min(vals_2, na.rm = TRUE) * 0.999999
    pulldown_2 <- pulldown_2 %>%
      mutate(across(all_of(cols_2), ~ replace(., is.any.na(.), impute_val_2)))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                        PART 2: MERGE DATASETS                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Split entries with multiple Protein_IDs
pulldown_1 <- pulldown_1 %>%
  mutate(Split_ID = strsplit(Protein_ID, ";"))
pulldown_2 <- pulldown_2 %>%
  mutate(Split_ID = strsplit(Protein_ID, ";"))

# Flatten the list of IDs into a single data frame
pulldown_1_flattened <- pulldown_1 %>%
  unnest(cols = Split_ID) %>%
  rename(Individual_ID = Split_ID)
pulldown_2_flattened <- pulldown_2 %>%
  unnest(cols = Split_ID) %>%
  rename(Individual_ID = Split_ID)

# Rename columns
pulldown_1_flattened <- pulldown_1_flattened %>%
  rename_with(~ paste0(., "_1"), -Individual_ID)
pulldown_2_flattened <- pulldown_2_flattened %>%
  rename_with(~ paste0(., "_2"), -Individual_ID)

# Merge the two data frames by Protein_ID
pulldown <- full_join(pulldown_1_flattened, pulldown_2_flattened, by = "Individual_ID")

# Remove duplicated rows
pulldown <- pulldown %>% 
  group_by(across(c(-Individual_ID, -Protein_ID_1, -Protein_ID_2))) %>%
  summarise(Protein_ID = paste(Individual_ID, collapse = ";"), .groups = "drop")

# Merge duplicated columns
pulldown <- pulldown %>%
  dplyr::mutate(Name = purrr::reduce(dplyr::select(., dplyr::matches("^Name_[0-9]+$")), dplyr::coalesce),
                Gene_ID = purrr::reduce(dplyr::select(., dplyr::matches("^Gene_ID_[0-9]+$")), dplyr::coalesce)) %>%
  dplyr::select(-dplyr::matches("^Name_[0-9]+$"), -dplyr::matches("^Gene_ID_[0-9]+$"))
if (highlight_proteins_red %in% c("kinases", "kinase")) {
  pulldown <- pulldown %>%
    dplyr::mutate(Kinase_interactor = purrr::reduce(dplyr::select(., dplyr::matches("^Kinase_interactor_[0-9]+$")), dplyr::coalesce),
                  Kinase = purrr::reduce(dplyr::select(., dplyr::matches("^Kinase_[0-9]+$")), dplyr::coalesce)) %>%
    dplyr::select(-dplyr::matches("^Kinase_interactor_[0-9]+$"), -dplyr::matches("^Kinase_[0-9]+$"))
}
if (highlight_proteins_blue %in% c("host kinases", "host kinase")) {
  pulldown <- pulldown %>%
    dplyr::mutate(Host_kinase_interactor = purrr::reduce(dplyr::select(., dplyr::matches("^Host_kinase_interactor_[0-9]+$")), dplyr::coalesce),
                  Host_kinase = purrr::reduce(dplyr::select(., dplyr::matches("^Host_kinase_[0-9]+$")), dplyr::coalesce)) %>%
    dplyr::select(-dplyr::matches("^Host_kinase_interactor_[0-9]+$"), -dplyr::matches("^Host_kinase_[0-9]+$"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                           PART 3: PLOT ENRICHMENT                            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Calculate log2 FC
if (has_H_1) { pulldown$log2FC_H_1 <- log2(pulldown$DMSO_H_1 / pulldown$Drug_H_1) }
if (has_L_1) { pulldown$log2FC_L_1 <- log2(pulldown$DMSO_L_1 / pulldown$Drug_L_1) }
if (has_H_1 & has_L_1) { pulldown <- pulldown %>% filter(.data[[paste0("log2FC_", load_level, "_1")]] != 0) } # Filter zero values
if (has_H_2) { pulldown$log2FC_H_2 <- log2(pulldown$DMSO_H_2 / pulldown$Drug_H_2) }
if (has_L_2) { pulldown$log2FC_L_2 <- log2(pulldown$DMSO_L_2 / pulldown$Drug_L_2) }
if (has_H_2 & has_L_2) { pulldown <- pulldown %>% filter(.data[[paste0("log2FC_", load_level, "_2")]] != 0) } # Filter zero values

# Calculate log2 FC values from observed values only
if (has_H_1) { raw_pulldown_1$log2FC_H_1_raw <- log2(raw_pulldown_1$DMSO_H / raw_pulldown_1$Drug_H) }
if (has_L_1) { raw_pulldown_1$log2FC_L_1_raw <- log2(raw_pulldown_1$DMSO_L / raw_pulldown_1$Drug_L) }
if (has_H_1 & has_L_1) { raw_pulldown_1 <- raw_pulldown_1 %>% filter(.data[[paste0("log2FC_", load_level, "_1_raw")]] != 0) } # Filter zero values
if (has_H_2) { raw_pulldown_2$log2FC_H_2_raw <- log2(raw_pulldown_2$DMSO_H / raw_pulldown_2$Drug_H) }
if (has_L_2) { raw_pulldown_2$log2FC_L_2_raw <- log2(raw_pulldown_2$DMSO_L / raw_pulldown_2$Drug_L) }
if (has_H_2 & has_L_2) { raw_pulldown_2 <- raw_pulldown_2 %>% filter(.data[[paste0("log2FC_", load_level, "_2_raw")]] != 0) } # Filter zero values

# Reformat highlighted protein list
highlight_kinases_red <- highlight_proteins_red %in% c("kinases", "kinase")
highlight_host_kinases_blue <- highlight_proteins_blue %in% c("host kinases", "host kinase")
highlight_IDs_red <- if (highlight_kinases_red) { character(0) } else { unlist(strsplit(highlight_proteins_red, "\\s*,\\s*")) }
highlight_IDs_blue <- if (highlight_host_kinases_blue) { character(0) } else { unlist(strsplit(highlight_proteins_blue, "\\s*,\\s*")) }

# Define visual styling
pulldown <- pulldown %>%
  plotly::mutate(Imputed_H_1 = if (add_imputation && all(c("DMSO_H_1","Drug_H_1") %in% colnames(.))) { (pulldown$Drug_H_1 == impute_val_1) | (pulldown$DMSO_H_1 == impute_val_1) } else { FALSE },
                 Imputed_L_1 = if (add_imputation && all(c("DMSO_L_1","Drug_L_1") %in% colnames(.))) { (pulldown$Drug_L_1 == impute_val_1) | (pulldown$DMSO_L_1 == impute_val_1) } else { FALSE },
                 Imputed_H_2 = if (add_imputation && all(c("DMSO_H_2","Drug_H_2") %in% colnames(.))) { (pulldown$Drug_H_2 == impute_val_2) | (pulldown$DMSO_H_2 == impute_val_2) } else { FALSE },
                 Imputed_L_2 = if (add_imputation && all(c("DMSO_L_2","Drug_L_2") %in% colnames(.))) { (pulldown$Drug_L_2 == impute_val_2) | (pulldown$DMSO_L_2 == impute_val_2) } else { FALSE },
                 Contaminant = grepl("cRAP", Protein_ID),
                 LowPeptide = (Peptides_1 <= 2 | Peptides_2 <= 2),
                 Shape = case_when(
                   (Imputed_H_1 | Imputed_L_1 | Imputed_H_2 | Imputed_L_2) & LowPeptide ~ "diamond-open",
                   (Imputed_H_1 | Imputed_L_1 | Imputed_H_2 | Imputed_L_2) & !LowPeptide ~ "diamond",
                   !(Imputed_H_1 | Imputed_L_1 | Imputed_H_2 | Imputed_L_2) & LowPeptide ~ "circle-open",
                   TRUE ~ "circle"
                 ),
                 Colour = case_when(
                   if ("Kinase" %in% colnames(pulldown)) { highlight_kinases_red & Kinase ~ "red" },
                   if ("Host_kinase" %in% colnames(pulldown)) { highlight_host_kinases_blue & Host_kinase ~ "blue" },
                   if ("Kinase_interactor" %in% colnames(pulldown) & highlight_interactors) { highlight_kinases_red & Kinase_interactor ~ "goldenrod" },
                   if ("Host_kinase_interactor" %in% colnames(pulldown) & highlight_interactors) { highlight_host_kinases_blue & Host_kinase_interactor ~ "green" },
                   Protein_ID %in% highlight_IDs_red ~ "red",
                   Protein_ID %in% highlight_IDs_blue ~ "blue",
                   Contaminant ~ "grey",
                   TRUE ~ "black"
                 ),
                 StyleGroup = paste(Shape, Colour, sep = "_")
  )

# Create a function to define the hit threshold
parse_hit_threshold <- function(threshold_input, values) {
  if (grepl("^median([+-])\\d*\\.?\\d*MAD$", threshold_input)) {
    matches <- regmatches(threshold_input, regexec("^median([+-])(\\d*\\.?\\d*)MAD$", threshold_input))[[1]]
    operator <- matches[2]
    multiplier <- as.numeric(matches[3])
    med_val <- median(values, na.rm = TRUE)
    mad_val <- mad(values, constant = 1, na.rm = TRUE)
    if (operator == "+") {
      return(med_val + multiplier * mad_val)
    } else {
      return(med_val - multiplier * mad_val)
    }
  } else if (grepl("^mean([+-])\\d*\\.?\\d*SD$", threshold_input)) {
    matches <- regmatches(threshold_input, regexec("^mean([+-])(\\d*\\.?\\d*)SD$", threshold_input))[[1]]
    operator <- matches[2]
    multiplier <- as.numeric(matches[3])
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    if (operator == "+") {
      return(mean_val + multiplier * sd_val)
    } else {
      return(mean_val - multiplier * sd_val)
    }
  } else if (!is.na(suppressWarnings(as.numeric(threshold_input)))) {
    return(as.numeric(threshold_input))
  }
}

# Calculate hit threshold
log2FC_threshold_H_1 <- FALSE
log2FC_threshold_L_1 <- FALSE
log2FC_threshold_H_2 <- FALSE
log2FC_threshold_L_2 <- FALSE

if (add_hit_threshold) {
  if (load_level == "H") { 
    log2FC_threshold_H_1 <- parse_hit_threshold(hit_threshold_file_1, raw_pulldown_1$log2FC_H_1_raw)
    log2FC_threshold_H_2 <- parse_hit_threshold(hit_threshold_file_2, raw_pulldown_2$log2FC_H_2_raw)
    }
  if (load_level == "L") { 
    log2FC_threshold_L_1 <- parse_hit_threshold(hit_threshold_file_1, raw_pulldown_1$log2FC_L_1_raw)
    log2FC_threshold_L_2 <- parse_hit_threshold(hit_threshold_file_2, raw_pulldown_2$log2FC_L_2_raw)
    }
}

# Initialise plot
fig <- plot_ly()

# Define legend labels
legend_labels <- c(
  "diamond_black" = "imputed",
  "circle_black" = "> 2 peptides",
  "circle-open_black" = "≤ 2 peptides",
  "circle_grey" = "contaminant",
  "circle_goldenrod" = "kinase interactor",
  "circle_green" = "host kinase interactor"
)

if ("circle_blue" %in% pulldown$StyleGroup) {
  legend_labels["circle_blue"] <- "highlighted"
} else if ("diamond_blue" %in% pulldown$StyleGroup) {
  legend_labels["diamond_blue"] <- "highlighted"
}

if ("circle_red" %in% pulldown$StyleGroup) {
  legend_labels["circle_red"] <- "highlighted"
} else if ("diamond_red" %in% pulldown$StyleGroup) {
  legend_labels["diamond_red"] <- "highlighted"
}

# Define display order (top = foreground, bottom = background)
style_order <- c(
  "diamond-open_red",        # highlighted_red, ≤ 2 peptides, imputed
  "diamond_red",             # highlighted_red, > 2 peptides, imputed
  "circle-open_red",         # highlighted_red, ≤ 2 peptides
  "circle_red",              # highlighted_red, ≤ 2 peptides
  "diamond-open_blue",       # highlighted_blue, ≤ 2 peptides, imputed
  "diamond_blue",            # highlighted_blue, > 2 peptides, imputed
  "circle-open_blue",        # highlighted_blue, ≤ 2 peptides
  "circle_blue",             # highlighted_blue, ≤ 2 peptides
  "diamond-open_goldenrod",  # highlighted_goldenrod, ≤ 2 peptides, imputed
  "diamond_goldenrod",       # highlighted_goldenrod, > 2 peptides, imputed
  "circle-open_goldenrod",   # highlighted_goldenrod, ≤ 2 peptides
  "circle_goldenrod",        # highlighted_goldenrod, ≤ 2 peptides
  "diamond-open_green",      # highlighted_green, ≤ 2 peptides, imputed
  "diamond_green",           # highlighted_green, > 2 peptides, imputed
  "circle-open_green",       # highlighted_green, ≤ 2 peptides
  "circle_green",            # highlighted_green, ≤ 2 peptides
  "diamond-open_grey",       # contaminant, ≤ 2 peptides, imputed
  "diamond_grey",            # contaminant, > 2 peptides, imputed
  "circle-open_grey",        # contaminant, ≤ 2 peptides
  "circle_grey",             # contaminant, > 2 peptides
  "diamond-open_black",      # ≤ 2 peptides, imputed
  "diamond_black",           # > 2 peptides, imputed
  "circle-open_black",       # ≤ 2 peptides
  "circle_black"             # > 2 peptides
)

# Add hit threshold lines
if (add_hit_threshold) {
  fig <- fig %>%
    add_trace(
      y = c(get(paste0("log2FC_threshold_", load_level, "_2")), max(pulldown[[paste0("log2FC_", load_level, "_2")]], na.rm = TRUE)),
      x = c(get(paste0("log2FC_threshold_", load_level, "_1")), get(paste0("log2FC_threshold_", load_level, "_1"))),
      type = "scatter",
      mode = "lines",
      line = list(color = "grey", dash = "dot", width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    ) %>%
    add_trace(
      y = c(get(paste0("log2FC_threshold_", load_level, "_2")), get(paste0("log2FC_threshold_", load_level, "_2"))),
      x = c(get(paste0("log2FC_threshold_", load_level, "_1")), max(pulldown[[paste0("log2FC_", load_level, "_1")]], na.rm = TRUE)),
      type = "scatter",
      mode = "lines",
      line = list(color = "grey", dash = "dot", width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    )
}

# Add data points
for (group in rev(style_order)) {
  group_data <- pulldown %>% dplyr::filter(StyleGroup == group)
  if (nrow(group_data) == 0) next
  
  parts <- strsplit(group, "_")[[1]]
  Shape <- parts[1]
  Color <- parts[2]
  
  marker_size <- ifelse(
    (Shape == "circle" && Color == "black") || 
      (Shape == "diamond" && Color == "black") || 
      (Shape == "circle" && Color == "red") ||
      (Shape == "diamond" && Color == "red") ||
      (Shape == "circle" && Color == "blue") ||
      (Shape == "diamond" && Color == "blue") ||
      (Shape == "circle" && Color == "goldenrod") ||
      (Shape == "diamond" && Color == "goldenrod") ||
      (Shape == "circle" && Color == "green") ||
      (Shape == "diamond" && Color == "green"),
    6, 4
  )
  
  legend_name <- legend_labels[group]
  show_legend <- !is.na(legend_name)
  
  if ("Gene_ID" %in% colnames(group_data)) {
    hover_text <- paste0(
      "Protein ID: ", group_data$Protein_ID,
      "<br>Gene ID: ", group_data$Gene_ID,
      "<br>Name: ", group_data$Name,
      "<br>Peptides (file 1):", group_data$Peptides_1,
      "<br>Peptides (file 2):", group_data$Peptides_2
    )
  } else {
    hover_text <- paste0(
      "Protein ID: ", group_data$Protein_ID,
      "<br>Name: ", group_data$Name,
      "<br>Peptides (file 1):", group_data$Peptides_1,
      "<br>Peptides (file 2):", group_data$Peptides_2
    )
  }
  if (has_H_1 & has_H_2) { hover_text <- paste0(hover_text, 
                                                "<br>log2FC_H (file 1): ", signif(group_data$log2FC_H_1, 3),
                                                "<br>log2FC_H (file 2): ", signif(group_data$log2FC_H_2, 3))}
  if (has_L_1 & has_L_2) { hover_text <- paste0(hover_text, 
                                                "<br>log2FC_L (file 1): ", signif(group_data$log2FC_L_1, 3),
                                                "<br>log2FC_L (file 2): ", signif(group_data$log2FC_L_2, 3)) }
  
  fig <- fig %>%
    add_trace(
      data = group_data,
      x = ~get(paste0("log2FC_", load_level, "_1")),
      y = ~get(paste0("log2FC_", load_level, "_2")),
      type = 'scatter',
      mode = 'markers',
      name = ifelse(show_legend, legend_name, ""),
      showlegend = show_legend,
      marker = list(
        symbol = Shape,
        size = marker_size,
        color = Color,
        line = list(color = Color, width = 1)
      ),
      text = hover_text,
      hoverinfo = 'text'
    )
}

# Finalise layout
fig <- fig %>%
  layout(
    xaxis = list(
      title = paste0("log2 FC DMSO/100 µM competitor ", load_level, " (H/L #1)"),
      #range = c(-12, 12), # <- Remove '#' sign before 'range' to specify x-axis limits
      showline = TRUE,
      linecolor = "black",
      ticks = "outside",
      tick0 = 0,
      dtick = 2,
      ticklen = 5,
      tickwidth = 1,
      tickcolor = "black",
      showgrid = FALSE
    ),
    yaxis = list(
      title = paste0("log2 FC DMSO/100 µM competitor ", load_level, " (H/L #2)"),
      #range = c(-12, 12), # <- Remove '#' sign before 'range' to specify y-axis limits
      showline = TRUE,
      linecolor = "black",
      ticks = "outside",
      tick0 = 0,
      dtick = 2,
      ticklen = 5,
      tickwidth = 1,
      tickcolor = "black",
      showgrid = FALSE
    ),
    legend = list(
      title = list(text = "Legend"),
      x = 1.05,
      y = 1,
      bordercolor = "grey",
      borderwidth = 1
    ),
    margin = list(l = 50, r = 50, t = 10, b = 60)
  )

# Display plot
print(fig)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                         PART 4: EXPORT RESULTS                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Format final dataset and create parameter dataframe
final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides_1", "DMSO_H_1", "Drug_H_1", "DMSO_L_1", "Drug_L_1", "Peptides_2", "DMSO_H_2", "Drug_H_2", "DMSO_L_2", "Drug_L_2", "log2FC_H_1", "log2FC_L_1", "log2FC_H_2", "log2FC_L_2", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
parameters <- data.frame(
  Parameter = c("data path", "data file 1", "data sheet 1", "date file 2", "data sheet 2", "script", "remove contaminants", "remove low peptides", "add imputation", "imputed value (LOD) file 1", "imputed value (LOD) file 2", "add hit threshold"),
  Value = c(data_path, pulldown_file_1, pulldown_sheet_1, pulldown_file_2, pulldown_sheet_2, "plot_common_hits_HL+DR.R", remove_contaminants, remove_low_peptides, add_imputation, impute_val_1, impute_val_2, add_hit_threshold),
  stringsAsFactors = FALSE
)

# Apply hit thresholds
if (add_hit_threshold) {
  if (load_level == "H") { 
    pulldown$Hit_1 <- pulldown$log2FC_H_1 > log2FC_threshold_H_1
    pulldown$Hit_2 <- pulldown$log2FC_H_2 > log2FC_threshold_H_2
  }
  if (load_level == "L") { 
    pulldown$Hit_1 <- pulldown$log2FC_L_1 > log2FC_threshold_L_1
    pulldown$Hit_2 <- pulldown$log2FC_L_2 > log2FC_threshold_L_2
  }
  final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides_1", "DMSO_H_1", "Drug_H_1", "DMSO_L_1", "Drug_L_1", "Peptides_2", "DMSO_H_2", "Drug_H_2", "DMSO_L_2", "Drug_L_2", "log2FC_H_1", "log2FC_L_1", "log2FC_H_2", "log2FC_L_2", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor", "Hit_1", "Hit_2")))
  parameters <- data.frame(
    Parameter = c("data path", "data file 1", "data sheet 1", "date file 2", "data sheet 2", "script", "remove contaminants", "remove low peptides", "add imputation", "imputed value (LOD) file 1", "imputed value (LOD) file 2", "add hit threshold", "hit threshold method file 1", "hit threshold H file 1", "hit threshold L file 1", "hit threshold method file 2", "hit threshold H file 2", "hit threshold L file 2"),
    Value = c(data_path, pulldown_file_1, pulldown_sheet_1, pulldown_file_2, pulldown_sheet_2, "plot_common_hits_HL+HL.R", remove_contaminants, remove_low_peptides, add_imputation, impute_val_1, impute_val_2, add_hit_threshold, hit_threshold_file_1, log2FC_threshold_H_1, log2FC_threshold_L_1, hit_threshold_file_2, log2FC_threshold_H_2, log2FC_threshold_L_2),
    stringsAsFactors = FALSE
  )
}

# Export to file
if (export_results) {
  write_xlsx(
    list("Results" = final_dataset, "Parameters" = parameters), results_filename)
}