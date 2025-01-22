# Omint: Omics Integration

Omint is a Shiny application designed for integrating omics data. It allows users to analyze relationships between different omic types: Lipid, Protein, and Metabolite. The application provides a ranked list of entities (proteins, lipids, or metabolites) based on their proximity to other omics, calculated using hyperbolic distance metrics. Users can input their data (e.g., Swisslipid IDs, UniProt IDs, PubChem IDs) and obtain ranked results along with downloadable CSV files for further analysis.

## Features
- **Omic Type Selection**: Choose between "Lipid", "Protein", or "Metabolite" to analyze the omic connections.
- **Distance Ranking**: Rank related omic entities based on proximity (e.g., Lipids vs Proteins, Metabolites vs Lipids).
- **Subset Processing**: Upload lists of omic IDs (e.g., lipid, protein, metabolite) to analyze.
- **Downloadable Results**: Download ranked results as CSV files.
- **Interactive UI**: User-friendly interface with progress bars, results tables, and valid subsets visualization.

## Requirements
- R version 4.0.0 or higher
- Required libraries:
  - `shiny`
  - `dplyr`

You can install these libraries with the following commands:
```r
install.packages("shiny")
install.packages("dplyr")
```

## How to Use
Download the "App.R" software and the file "nodes_multi_omics_category2.csv" in Git. Change the directory in the code as follows:
      data <- read.csv('Directory for nodes_multi_omics_category2.csv')
      
Then click "Run App".
1. **Input Data**:
   - Select the omic type (Lipid, Protein, or Metabolite).
   - Input comma-separated omic IDs for lipids, proteins, or metabolites.
   - Specify the number of smallest distances to consider.

2. **Processing**:
   - Click the "Process" button to start processing the data. The app will calculate the distances and rank the entities based on their connection to the selected omic type.

3. **Results**:
   - The ranked list of entities will be displayed.
   - Valid subsets (entities that are valid for the selected omic type) will also be shown.
   - You can download the results as CSV files for each entity type (e.g., "Lipid (Proteins)", "Metabolites (Proteins)").

## Data Input Example
- Lipids with Swisslipid IDs: `SLM:000000347,SLM:000000505,SLM:000000795`
- Proteins with UniProt IDs: `A0A0B4J2D5,A0AV02,A0AV96`
- Metabolites with PubChem IDs: `pubchem:863,pubchem:87642,pubchem:171548`

## Code Overview
- **UI Components**: Includes radio buttons, numeric inputs, text inputs for omic IDs, and action buttons.
- **Server Logic**:
  - The server processes the input data, calculates distance rankings using preloaded CSV distance matrices (`lipid_protein_hyperbolicDistance2.csv`, `protein_metabolite_hyperbolicDistance2.csv`, `metabolite_lipid_hyperbolicDistance2.csv`).
  - Ranks are computed based on the smallest distance values and are displayed in sorted order.
  - Downloadable buttons allow users to save results.

