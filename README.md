# Omint: Omics Integration

Omint is a Shiny application designed for integrating omics data. It allows users to analyze relationships between different omic types: Lipid, Protein, and Metabolite. The application provides a ranked list of entities (proteins, lipids, or metabolites) based on their proximity to other omics, calculated using hyperbolic distance metrics. Users can input their data (e.g., Swisslipid IDs, UniProt IDs, PubChem IDs) and obtain ranked results along with downloadable CSV files for further analysis.

## Features
- **Tissue Locations Selection**: Choose Tissue location(s) from a drop-down list of Tissues, to indicate the cell.
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
Download the "App.R" software and the files "nodes_multi_omics_category2.csv" and "protein_tissueLocation_clean.csv" in Git. Change the directory in the code as follows:
      
      protein_tissueLocation_map <- read.csv('Directory for protein_tissueLocation_clean.csv')
      data <- read.csv('Directory for nodes_multi_omics_category2.csv')
      
Then click "Run App".
1. **Input Data**:
   - Select the Tissue Locations (e.g. kidney, liver, lung).
   - Select the omic type (Lipid, Protein, or Metabolite).
   - Input comma-separated omic IDs for lipids, proteins, or metabolites.
   - Specify the number of smallest distances to consider. NOTE: it should be less than or equal to the number of omic IDs above.

2. **Processing**:
   - Click the "Process" button to start processing the data. The app will calculate the distances and rank the entities based on their connection to the selected omic type.

3. **Results**:
   - The ranked list of entities will be displayed.
   - Valid subsets (entities that are valid for the selected omic type) will also be shown.
   - You can download the results as CSV files for each entity type (e.g., "Lipid (Proteins)", "Metabolites (Proteins)").

## Example Usage: Cardiovascular Disease Associations
This example demonstrates how users can analyze proteins in the heart muscle associated with cardiovascular diseases and obtain a ranked list of related lipids and metabolites.

### Step-by-Step Instructions
1. **Launch the Application**:  
   Open RStudio, set your working directory to the folder containing "App.R" and "nodes_multi_omics_category2.csv", and click **Run App**.

2. **Select Tissue Locations**:  
   In the user interface, choose **heart muscle** as your Tissue Location.

3. **Select Input Omic Type**:  
   In the user interface, choose **Protein** as your input omic type.

4. **Input Protein IDs**:  
   Enter the following comma-separated list of protein identifiers (UniProt IDs) associated with cardiovascular diseases:
   ```
   P51569,Q60612,P35951,Q61009,P08226,O14649,Q9QUK6,P06804,P16671,P14901
   ```

5. **Specify Parameters**:  
   Set the number of smallest distances to consider (default is 3). Ensure that this number does not exceed the count of input proteins.

6. **Process the Data**:  
   Click the **Process** button. The application will:
   - Validate the provided protein IDs.
   - Compute the hyperbolic distances between these proteins and all lipids/metabolites in the dataset.
   - Rank the lipids and metabolites based on the computed association scores.

7. **View the Results**:  
   Once processing is complete, you will see:
   - A ranked list of lipids and metabolites that are most strongly associated with the input proteins.
   - Evidence details for the top associations, showing the contributing distances from each protein.

8. **Download the Output**:  
   Use the download buttons to export the ranked lists as CSV files for further analysis.


## Code Overview
- **UI Components**: Includes radio buttons, numeric inputs, text inputs for omic IDs, and action buttons.
- **Server Logic**:
  - The server processes the input data and calculates distance rankings.
  - Ranks are computed based on the smallest distance values and are displayed in sorted order.
  - Downloadable buttons allow users to save results.

For citations:
Anyaegbunam, U. A., Vagiona, A.-C., ten Cate, V., Bauer, K., Schmidlin, T., Distler, U., Tenzer, S., Araldi, E., Bindila, L., Wild, P., & Andrade-Navarro, M. A. (2025). A Map of the Lipid–Metabolite–Protein Network to Aid Multi-Omics Integration. Biomolecules, 15(4), 484. https://doi.org/10.3390/biom15040484
