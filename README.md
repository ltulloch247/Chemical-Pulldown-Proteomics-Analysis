[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24286730.svg)](https://doi.org/10.5281/zenodo.24286730)

# Chemical Pulldown Proteomics Analysis Using R

## Overview

This repository contains R scripts for the analysis and visualisation of mass spectrometry-based chemical pulldown proteomics data.

The scripts are designed to identify proteins that interact with a bead-immobilised drug (affinity-based chemical probe) by quantitatively comparing proteins captured in the absence of competitor (DMSO vehicle control) with those captured in the presence of free drug, which competes with the immobilised compound for protein binding.

The `plot_enrichment_HL.R` script compares protein intensities between samples with and without competitor and calculates enrichment as the log₂ fold-change (DMSO/competitor). High values indicate strong competition-dependent enrichment, while values close to zero indicate little or no detectable competition. The script can identify significantly enriched proteins using user-defined thresholds, such as a specified number of median absolute deviations (MAD) above the median enrichment.

Two biological replicates can be compared using the `plot_common_hits_HL+HL.R` script, which plots enrichment values from replicate 1 against those from replicate 2, highlighting proteins that are consistently enriched above the specified thresholds in both experiments.

Dose-response pulldowns, in which a range of competitor concentrations is used, can be analysed using the `plot_DRC.R` script. This fits four-parameter dose-response curves to protein intensity measurements across the concentration series and calculates IC₅₀ values.

The scripts are written in R and can be run using RStudio. The `install_packages.R` script should be run first to install the required R packages.

## Usage

Each script contains a user-configurable section near the beginning, where input files, data columns and analysis parameters can be specified.

To perform an analysis:

1. Open the relevant `.R` script in RStudio.
2. Set `data_path` to the directory containing the input data.
3. Specify the input filename (`.xlsx` or `.csv`) and, where applicable, the Excel worksheet.
4. Update the column-name mappings to match the input dataset.
5. Adjust the analysis parameters as required, following the comments within the script.
6. Run the script using **Source** in RStudio.

Interactive enrichment plots can be viewed in RStudio and exported as images or HTML files. Processed results can be exported to Excel, while dose-response curves are automatically saved as PNG and SVG files.

## Enrichment analysis

Protein enrichment is calculated as:

`log2FC = log2(DMSO intensity / competitor intensity)`

Higher values indicate greater competition-dependent enrichment, consistent with selective capture of the protein by the immobilised compound.

Enrichment thresholds can be specified as a fixed log₂ fold-change, mean plus a multiple of the standard deviation (SD), or median plus a selected number of median absolute deviations (MAD). We recommend using six MAD above the median (`median+6MAD`).

In some cases, competition may be sufficiently effective that a protein is no longer detected in the competitor-treated sample. Without imputation, the calculated enrichment may be undefined or infinite, preventing the protein from being displayed appropriately. To include these proteins in the analysis, missing intensity values can optionally be imputed using a value approximating the limit of detection (LOD).

Results are presented as enrichment plots and can also be exported to Excel.

## Dose-response analysis

Dose-response curves are fitted using four-parameter log-logistic regression implemented in the R package `drc`.

**Drug concentrations:** The competitor concentrations can be modified on line 112 of `plot_DRC.R` (under `# Extract drug concentrations`) to match the experimental design.

**DMSO control:** For plotting and curve fitting on a logarithmic concentration axis, the DMSO control is assigned an artificial concentration equal to one-tenth of the lowest non-zero competitor concentration (line 114). This is a plotting and modelling convention and does not represent an actual drug concentration in the DMSO control.

**Imputed data:** Where a protein is not detected in a sample, its intensity can optionally be imputed using a value approximating the LOD to facilitate curve fitting.

## Data interpretation

Protein enrichment reflects the reduction in protein recovery following competition with an unlabelled compound and does not necessarily establish direct compound–protein binding.

Missing intensity values can optionally be imputed using a value slightly below the lowest detected intensity in the dataset. Enrichment values involving imputed measurements should be interpreted cautiously, as highly abundant proteins may exhibit greater calculated enrichment simply because of their higher initial intensity.

Reproducible enrichment across independent experiments and concentration-dependent competition provide additional evidence supporting potential compound–protein interactions.

## Authors

- **Laurine Brouck**
- **Susan Wyllie**

University of Dundee, UK.

## Citation

If these scripts are used in published research, please cite the associated software release.

**Zenodo DOI:** To be added following deposition.
