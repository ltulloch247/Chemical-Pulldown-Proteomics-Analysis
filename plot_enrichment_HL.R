# Load packages
library(tools)
library(readxl)
library(dplyr)
library(plotly)
library(writexl)

# Define analysis parameters
data_path <- "C:/Users/path-to-data"                              # Replace path with your owb path. Use / not \
pulldown_file <- "DNDI-6690_pulldown_Low_Rep1.xlsx"               # Enter filename and include extension (.xlsx/.csv)
pulldown_sheet <- ""                                              # Leave as "" if only one sheet
load_level <- "L"                                                 # Choose "H" or "L" 
highlight_proteins_red <- "A0A3Q8IRJ0"                            # Protein IDs or "kinases" (pre-processed data only)
highlight_proteins_blue <- ""                                     # Protein IDs or "host kinases" (pre-processed data only)
highlight_interactors <- FALSE                                    # TRUE = yes, FALSE = no (pre-processed data only)
remove_contaminants <- TRUE                                       # TRUE = yes, FALSE = no
remove_low_peptides <- TRUE                                       # TRUE = yes, FALSE = no
add_imputation <- FALSE                                           # TRUE = yes, FALSE = no
add_hit_threshold <- TRUE                                         # TRUE = yes, FALSE = no
hit_threshold <- "median+6MAD"                                    # Choose format: "2", "mean+2SD", "median+3MAD"
export_results <- TRUE                                            # TRUE = yes, FALSE = no
results_filename <- "DNDI-6690_pulldown_Low_Rep1-enrichment.xlsx"    # Include .xlsx extension

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
    Gene_ID = "Gene_ID",          # Use "" if no gene information
    Peptides = "Peptides",
    Intensity = "Intensity",      # Use "" if no intensity column - it will be calculated automatically
    DMSO_H = "",                  # Use "" if no H data
    Drug_H = "",                  # Use "" if no H data
    DMSO_L = "DMSO_L",            # Use "" if no L data
    Drug_L = "Comp_L"             # Use "" if no L data
  )))

# Detect available data
has_H <- all(c("DMSO_H", "Drug_H") %in% colnames(pulldown))
has_L <- all(c("DMSO_L", "Drug_L") %in% colnames(pulldown))

# Clean and prepare data
pulldown <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO_H", "Drug_H", "DMSO_L", "Drug_L", "Intensity", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
pulldown <- pulldown %>% dplyr::mutate(across(any_of(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L")), as.numeric))

# Determine data type
cols <- intersect(c("DMSO_H", "Drug_H", "DMSO_L", "Drug_L"), colnames(pulldown))
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

# Freeze raw data for hit threshold calculations
raw_pulldown <- pulldown

# Remove contaminants
if (remove_contaminants) {
  pulldown <- dplyr::filter(pulldown, !grepl("cRAP", Protein_ID))
  raw_pulldown <- dplyr::filter(raw_pulldown, !grepl("cRAP", Protein_ID))
}

# Remove proteins with low peptide counts (≤ 2)
if (remove_low_peptides) {
  pulldown <- pulldown %>% dplyr::filter(Peptides > 2)
  raw_pulldown <- raw_pulldown %>% dplyr::filter(Peptides > 2)
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
#                        PART 2: PLOT ENRICHMENT                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Calculate log2 FC
if (has_H) { pulldown$log2FC_H <- log2(pulldown$DMSO_H / pulldown$Drug_H) }
if (has_L) { pulldown$log2FC_L <- log2(pulldown$DMSO_L / pulldown$Drug_L) }
if (has_H & has_L) { pulldown <- pulldown %>% filter(.data[[paste0("log2FC_", load_level)]] != 0) } # Filter zero values

# Calculate log2 FC values from observed values only
if (has_H) { raw_pulldown$log2FC_H_raw <- log2(raw_pulldown$DMSO_H / raw_pulldown$Drug_H) }
if (has_L) { raw_pulldown$log2FC_L_raw <- log2(raw_pulldown$DMSO_L / raw_pulldown$Drug_L) }
if (has_H & has_L) { raw_pulldown <- raw_pulldown %>% filter(.data[[paste0("log2FC_", load_level, "_raw")]] != 0) } # Filter zero values

# Reformat highlighted protein list
highlight_kinases_red <- highlight_proteins_red %in% c("kinases", "kinase")
highlight_host_kinases_blue <- highlight_proteins_blue %in% c("host kinases", "host kinase")
highlight_IDs_red <- if (highlight_kinases_red) { character(0) } else { unlist(strsplit(highlight_proteins_red, "\\s*,\\s*")) }
highlight_IDs_blue <- if (highlight_host_kinases_blue) { character(0) } else { unlist(strsplit(highlight_proteins_blue, "\\s*,\\s*")) }

# Define visual styling
pulldown <- pulldown %>%
  plotly::mutate(Imputed_H = if (add_imputation && all(c("DMSO_H","Drug_H") %in% colnames(.))) {(Drug_H == impute_val) | (DMSO_H == impute_val)} else FALSE,
                 Imputed_L = if (add_imputation && all(c("DMSO_L","Drug_L") %in% colnames(.))) {(Drug_L == impute_val) | (DMSO_L == impute_val)} else FALSE,
                 Contaminant = grepl("^cRAP", Protein_ID),
                 LowPeptide = Peptides <= 2,
                 Shape = case_when(
                   .data[[paste0("Imputed_", load_level)]] & LowPeptide ~ "diamond-open",
                   .data[[paste0("Imputed_", load_level)]] & !LowPeptide ~ "diamond",
                   !.data[[paste0("Imputed_", load_level)]] & LowPeptide ~ "circle-open",
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
log2FC_threshold_H <- FALSE
log2FC_threshold_L <- FALSE

if (add_hit_threshold) {
  if (load_level == "H") { log2FC_threshold_H <- parse_hit_threshold(hit_threshold, raw_pulldown$log2FC_H_raw) }
  if (load_level == "L") { log2FC_threshold_L <- parse_hit_threshold(hit_threshold, raw_pulldown$log2FC_L_raw) }
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

# Add hit threshold line
if (add_hit_threshold) {
  fig <- fig %>%
    layout(
      shapes = list(
        list(
          type = "line",
          x0 = get(paste0("log2FC_threshold_", load_level)),
          x1 = get(paste0("log2FC_threshold_", load_level)),
          y0 = 0,
          y1 = 1,
          yref = "paper",
          line = list(color = "black", dash = "dot", width = 0.5)
        )
      )
    )
}

# Add data points
for(group in rev(style_order)) {
  group_data <- pulldown %>% filter(StyleGroup == group)
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
      "<br>Peptides: ", group_data$Peptides
    )
  } else {
    hover_text <- paste0(
      "Protein ID: ", group_data$Protein_ID,
      "<br>Name: ", group_data$Name,
      "<br>Peptides: ", group_data$Peptides
    )
  }
  if (has_H) { hover_text <- paste0(hover_text, "<br>log2FC_H: ", signif(group_data$log2FC_H, 3)) }
  if (has_L) { hover_text <- paste0(hover_text, "<br>log2FC_L: ", signif(group_data$log2FC_L, 3)) }
  
  fig <- fig %>%
    add_trace(
      data = group_data,
      x = ~get(paste0("log2FC_", load_level)),
      y = ~Intensity,
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
      title = paste0("log2 FC DMSO/100 µM competitor ", load_level),
      #range = c(-14, 14), # <- Remove '#' sign before 'range' to set x-axis limits
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
      title = "Intensity", type = "log",
      #range = c(3, 11), # <- Remove '#' sign before 'range' to set y-axis limits (log scale)
      tickvals = 10^(3:11),
      tickformat = ".1e",
      showline = TRUE,
      linecolor = "black",
      ticks = "outside",
      ticklen = 5,
      tickwidth = 1,
      tickcolor = "black",
      showticksuffix = "none",
      showexponent = "all",
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
#                         PART 3: EXPORT RESULTS                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Format final dataset and create parameter dataframe
final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO_H", "Drug_H", "DMSO_L", "Drug_L", "Intensity", "log2FC_H", "log2FC_L", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor")))
parameters <- data.frame(
  Parameter = c("data path", "data file", "data sheet", "script", "remove contaminants", "remove low peptides", "add imputation", "imputed value (LOD)", "add hit threshold"),
  Value = c(data_path, pulldown_file, pulldown_sheet, "plot_enrichment_HL.R", remove_contaminants, remove_low_peptides, add_imputation, impute_val, add_hit_threshold),
  stringsAsFactors = FALSE
)

# Apply hit thresholds
if (add_hit_threshold) {
  if (load_level == "H") { pulldown$Hit <- pulldown$log2FC_H > log2FC_threshold_H }
  if (load_level == "L") { pulldown$Hit <- pulldown$log2FC_L > log2FC_threshold_L }
  
  final_dataset <- pulldown %>% dplyr::select(any_of(c("Name", "Protein_ID", "Gene_ID", "Peptides", "DMSO_H", "Drug_H", "DMSO_L", "Drug_L", "Intensity", "log2FC_H", "log2FC_L", "Kinase", "Kinase_interactor", "Host_kinase", "Host_kinase_interactor", "Hit")))
  
  parameters <- data.frame(
    Parameter = c("data path", "data file", "data sheet", "script", "remove contaminants", "remove low peptides", "add imputation", "imputed value (LOD)", "add hit threshold", "hit threshold method", "hit threshold H", "hit threshold L"),
    Value = c(data_path, pulldown_file, pulldown_sheet, "plot_enrichment_HL.R", remove_contaminants, remove_low_peptides, add_imputation, impute_val, add_hit_threshold, hit_threshold,
              if (load_level == "H") log2FC_threshold_H else NA,
              if (load_level == "L") log2FC_threshold_L else NA),
    stringsAsFactors = FALSE
  )
}

# Export to file
if (export_results) {
  write_xlsx(
    list("Results" = final_dataset, "Parameters" = parameters), results_filename)
}
