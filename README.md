# 🔬 Crosslink Analyzer for XiNET / Output: Proteome Discoverer & pLink

A flexible R-based pipeline for processing, filtering, and organizing crosslinking mass spectrometry (XL-MS) data for downstream visualization and analysis (e.g., XiNET-compatible outputs).

---

## 📌 Overview

This tool is designed to standardize and streamline the analysis of crosslinked peptide datasets generated from:

* Proteome Discoverer (PD)
* pLink
* XiNET-compatible workflows

It enables:

* 📊 Score-based filtering of crosslinks
* 🔁 Replicate consistency filtering
* 🧪 Grouping by experimental conditions
* 📁 Structured output generation for downstream analysis

---

## ⚙️ Features

* Supports multiple input formats (PD / pLink)
* Customizable filtering thresholds
* Flexible replicate handling (per treatment or global)
* Modular structure for easy adaptation
* Clean, reproducible workflow in R

---

## 📂 Project Structure

```
project/
│── data/                  # Input files (crosslink datasets)
│── database/              # FASTA database
│── output/                # Generated results
│── script.R               # Main analysis script
│── README.md              # Documentation
```

---

## 🚀 Getting Started

### 1. Install Requirements

Make sure you have R installed, then install the required packages:

```r
install.packages(c(
  "tidyverse", "dplyr", "data.table",
  "stringr", "readr", "purrr", "seqinr"
))

# Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
```

---

### 2. Configure the Script

Open `script.R` and modify the configuration section:

* Set your working directory
* Provide the path to your FASTA database
* Define input/output folders
* Adjust filtering parameters
* Define treatment groups and replicates

Example:

```r
database_path <- "path/to/database.fasta"
input_dir <- "path/to/input_files"
output_dir <- "path/to/output"

min_score <- 40
min_replicas <- 2
replica_mode <- "Treatment"
```

---

### 3. Define Treatment Groups

Group your replicates by experimental condition:

```r
treatment_groups <- list(
  "Condition_A" = c("rep1", "rep2", "rep3"),
  "Condition_B" = c("rep4", "rep5", "rep6")
)
```

---

### 4. Run the Analysis

Simply run the script in R:

```r
source("script.R")
```

---

## 📊 Input Requirements

* Crosslink data files in `.csv` or compatible format
* Consistent file naming across replicates
* Matching identifiers between files and `treatment_groups`
* FASTA database for sequence reference

---

## 📈 Output

The pipeline generates:

* Filtered crosslink datasets
* Replicate-consistent interactions
* Organized files ready for visualization tools (e.g., XiNET)

Outputs are saved in the specified output directory.

---

## 🧠 Methodology

The workflow applies:

1. **Initial filtering** (optional) based on score threshold
2. **Replication filtering**

   * By treatment group OR across all samples
3. **Data restructuring** for compatibility with downstream tools

---

## 🔧 Customization

You can easily extend this pipeline to:

* Add new filtering criteria
* Integrate statistical analysis
* Export to additional formats
* Automate batch processing
* Contributions are welcome!
* ## 📖 Reference
  

If you use this pipeline in your work, please cite:

Frances, N., Giustini, C., Finazzi, G., Ferro, M., & Albanese, P.  
**Deciphering Photosynthetic Protein Networks: A Crosslinking-MS Strategy for Studying Functional Thylakoid Membranes**  
*bioRxiv* (2025)  
https://doi.org/10.1101/2025.10.07.681025  

🔗 https://www.biorxiv.org/content/10.1101/2025.10.07.681025v1
