# Pesticide Analysis Pipeline

## Overview

This script processes a multi-sheet Excel workbook containing pesticide observation data and generates a comprehensive set of analytical outputs, including:

- Binary presence/absence matrices
- Relative frequency summaries
- t-SNE dimensionality reduction visualizations
- Pesticide co-occurrence pairs and triplets
- Association rule mining metrics (support and confidence)
- Pair-support heatmaps
- Hamming distance analyses
- Excel exports and LaTeX-ready tables

The workflow is designed to run end-to-end from a single input workbook.

## Input

The script expects an Excel file named:

```text
data.xlsx
```

Each worksheet is treated as an individual dataset/locality and is loaded into a pandas DataFrame for analysis.

## Output

All generated results are written to:

```text
output/
```

The directory is created automatically if it does not already exist.

## Main Processing Steps

### 1. Data Preparation

- Load all workbook sheets.
- Clean worksheet names and filenames.
- Convert pesticide data into binary presence/absence matrices.
- Export consolidated binary datasets.

### 2. Dimensionality Reduction

- Generate t-SNE embeddings from binary data.
- Explore similarity patterns between locality-month observations.

### 3. Relative Frequency Analysis

- Calculate pesticide occurrence frequencies.
- Produce locality-based frequency summaries.

### 4. Itemset Mining

- Identify frequently occurring pesticide pairs.
- Identify frequently occurring pesticide triplets.

### 5. Association Rule Analysis

- Compute support and confidence metrics.
- Generate pairwise association statistics.
- Generate triplet-based association rules.
- Filter rules using configurable thresholds.

Default thresholds:

```python
MIN_SUPPORT = 0.1
MIN_CONFIDENCE = 0.5
```

### 6. Visualization

- Pair-support heatmaps
- Clustering outputs
- Distance-based analyses

### 7. LaTeX Export

- Generate LaTeX tables for reports and publications.
- Export top association rules by locality.

## Configuration

Global parameters:

```python
TOP_N_PAIRS = 50
TOP_N_PAIR_HEATMAP = 20
TOP_N_PESTICIDES = 10
TOP_N_RULES = 50
MIN_SUPPORT = 0.1
MIN_CONFIDENCE = 0.5
```

## Dependencies

The script uses the following Python packages:

```text
numpy
pandas
matplotlib
seaborn
scipy
scikit-learn
openpyxl
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Running the Script

Place `data.xlsx` in the same directory as the script and run:

```bash
python main.py
```

## Generated Artifacts

Examples of generated artifacts include:

- Binary datasets
- Frequency tables
- Association rule spreadsheets
- Heatmap images
- Hamming distance outputs
- LaTeX tables for publication-ready reporting

## Entry Point

The pipeline is executed through:

```python
if __name__ == '__main__':
    main()
```

## Notes

- Excel worksheet names are automatically sanitized to comply with Excel naming limits.
- Output directories are created automatically.
- The pipeline is intended for exploratory analysis of pesticide occurrence and co-occurrence patterns across multiple localities and sampling periods.

## Original file
The attached file
```text
data.xlsx
```
represents a test set of artificially generated data. The original dataset can be obtained on demand. However, we provide a file with pair rules generated from the original data used in the paper: 
```text
original_pair_rules_metrics_by_locality.xlsx
```

## Licence
Creative Commons Attribution 4.0 International, Copyright (C) 2026 VSB - Technical University of Ostrava