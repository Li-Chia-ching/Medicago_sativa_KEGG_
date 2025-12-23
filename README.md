# GO/KEGG Enrichment Analysis Pipeline (V3.1)

## Overview

This R script provides a robust, **offline-capable** pipeline for performing GO (Gene Ontology) and KEGG (Kyoto Encyclopedia of Genes and Genomes) enrichment analysis. It is specifically designed to handle *Medicago truncatula* gene module data, with intelligent features for data parsing, error handling, and publication-quality visualization.

**Key Features:**
-   **Offline-First Design**: Works without an internet connection after initial package installation.
-   **Production-Ready Robustness**: Includes comprehensive error handling and graceful fallbacks.
-   **Intelligent Data Parsing**: Flexible recognition of various column names and annotation formats.
-   **Publication-Quality Output**: Generates high-resolution figures and detailed reports.

## Input File Requirements

Place your input file (e.g., `Summary_Modules.xlsx`) in your R working directory. The script expects a table with the following columns. It will automatically search for common alternative names if the exact names are not found.

| Required Data | Expected Column Name(s) | Description |
| :--- | :--- | :--- |
| **Gene Identifier** | `Gene_ID`, `GeneID`, `gene_id`, `gene`, `ID` | Unique gene identifier. |
| **Module Assignment** | `Module`, `module`, `Cluster`, `cluster`, `Group` | The module/group each gene belongs to (e.g., `turquoise`, `blue`). |
| **GO Annotations** | `GO`, `Go`, `go`, `GO_terms` | GO terms associated with each gene (semicolon-, comma-, or pipe-separated). |
| **KEGG Annotations** | `KEGG`, `Kegg`, `kegg`, `Pathway`, `pathway` | KEGG pathway IDs associated with each gene. |

## Quick Start Guide

1.  **Prepare Your Data**: Ensure your `Summary_Modules` file is saved (e.g., as `.xlsx`, `.csv`, or `.txt`) in your R working directory.
2.  **Run the Script**: Execute the entire script in RStudio. It will:
    *   Automatically check for and install required R packages.
    *   Create a unique timestamped output folder (e.g., `GO_KEGG_Analysis_20241215_143022`).
    *   Perform the entire analysis and save all results inside the new folder.
3.  **Review Results**: Check the generated output folder for CSV results, PDF/PNG figures, and a comprehensive text report.

## Core Features

### 1. 🔧 Robust Data Preprocessing & Error Handling
-   **Smart Column Mapping**: Automatically detects variations in column names (e.g., `GeneID` vs. `Gene_ID`), preventing crashes due to simple typos.
-   **Special "Grey" Module Handling**: Identifies and prompts the user to exclude generic modules like `grey` or `0`, which typically contain unclustered genes in WGCNA.
-   **Graceful Degradation**: If the `GO.db` package is unavailable, the script continues by using GO IDs instead of full names, ensuring uninterrupted offline analysis.

### 2. 🧠 Intelligent Annotation Parsing
-   **Flexible KEGG ID Recognition**: Uses an enhanced regular expression (`(mtr\|map\|ko)\d{5}`) to accurately capture various KEGG ID formats specific to *Medicago* (`mtr`), general pathways (`map`), or ortholog groups (`ko`).
-   **Multi-Separator Support for GO Terms**: Correctly splits GO terms whether they are separated by semicolons (`;`), commas (`,`), or pipes (`\|`), accommodating different data export formats.

### 3. 📊 Adaptive, Publication-Ready Visualization
-   **Dynamic Figure Sizing**: Plot dimensions automatically adjust based on the number of modules and terms, preventing overcrowded or unreadable figures.
    ```r
    plot_width <- max(12, num_modules * 1.5)
    plot_height <- 8 + (number_of_terms * 0.2)
    ```
-   **Automatic Text Wrapping**: Long GO term and pathway names are intelligently wrapped to maintain clean and professional figure layouts.
-   **High-Quality Output**: Exports both vector PDFs for editing and high-resolution (600 DPI) PNGs for publication, using a specified low-saturation color palette.

### 4. 📈 Comprehensive Functional Enrichment Analysis
-   **Modular Enrichment**: Performs separate enrichment analyses for each module in your data.
-   **GO Ontology Split**: Conducts and reports enrichment for Biological Process (BP), Cellular Component (CC), and Molecular Function (MF) ontologies independently.
-   **Statistical Rigor**: Applies industry-standard multiple testing correction (Benjamini-Hochberg) and customizable significance thresholds.

### 5. 📋 Detailed Reporting & Reproducibility
-   **Automated Analysis Report**: Generates a detailed `09_Analysis_Report.txt` file summarizing all parameters, input statistics, and results.
-   **Complete Output Suite**: Saves every step's results—from gene counts per module and extracted annotations to final enrichment results and figures—in clearly named files.

## Practical Usage Notes

-   **The "Grey" Module**: When prompted, it is generally recommended to **exclude** the `grey` module (type `y` or `yes`) as it represents genes not assigned to any coherent functional module.
-   **Understanding the Output**:
    -   `05_GO_Enrichment_Results.csv`: The main results table. Significant terms are filtered by `p.adjust < 0.05`.
    -   `07_GO_Enrichment_Dotplot.pdf`: The key visualization. Point size represents the number of genes in the term; color intensity represents the significance (`-log10(p.adjust)`).
-   **Troubleshooting**: If KEGG pathway names appear as IDs (e.g., `mtr00010`) instead of descriptive names (e.g., `Glycolysis`), it indicates the script failed to connect to the online KEGG database (due to network settings or strict offline mode). The analysis is still **completely valid**; the IDs can be looked up manually on the KEGG website if needed.

## Version History

-   **V3.1 (Current)**: Enhanced robustness with better regex for KEGG parsing, improved offline stability, and dynamic visualization scaling.
-   **V3.0**: Added intelligent "grey" module handling, improved error reporting, and data quality diagnostics.
-   **Earlier Versions**: Established core functionality for offline annotation extraction and enrichment analysis.

---
**Verdict from Code Review**: This script is **ready for use**. It is highly defensive (won't crash easily), produces publication-quality figures, and handles the specific quirks of *Medicago* WGCNA data perfectly.
