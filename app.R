library(shiny)
library(dplyr)

# Define UI for the application
ui <- fluidPage(
  titlePanel("Omint: Omics Integration"),
  p("Get a ranked list of other omics based on its connection to one omic list."),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons(
        "omic_type",
        "Select Omic Type:",
        choices = list("Lipid" = "lipid", "Protein" = "protein", "Metabolite" = "metabolite"),
        selected = "lipid"
      ),
      numericInput("n", "Number of smallest distances to consider:", value = 3, min = 1),
      textInput("lipids_input", "Lipids with Swisslipid id (comma-separated):", "SLM:000000347,SLM:000000505,SLM:000000795"),
      textInput("proteins_input", "Proteins with UniProt id (comma-separated):", "A0A0B4J2D5,A0AV02,A0AV96"),
      textInput("metabolites_input", "Metabolites with PubChem id (comma-separated):", "pubchem:863,pubchem:87642,pubchem:171548"),
      actionButton("submit", "Process"),
      uiOutput("downloadButtons")
    ),
    
    mainPanel(
      uiOutput("resultsPanel"),
      uiOutput("validSubsetPanel")  # New output to display the valid_subset
    )
  )
)

# Define server logic required to rank entities
server <- function(input, output) {
  
  result_data <- reactiveVal(NULL)
  valid_subsets_data <- reactiveVal(NULL)  # Store the valid subsets
  
  # Function to clean, deduplicate, and split the user input
  process_input <- function(input_string) {
    unique(trimws(unlist(strsplit(input_string, ","))))
  }
  
  # Function to get ranked list of entities from a given subset
  get_ranked_list <- function(subset, df, n, key, match_on = "index") {
    if (match_on == "index") {
      valid_subset <- subset[subset %in% rownames(df)]
    } else {
      valid_subset <- subset[subset %in% colnames(df)]
    }
    
    if (length(valid_subset) == 0) {
      return(list(valid_subset = valid_subset, data = data.frame(entity = character(), rank = numeric(), evidence = character(), stringsAsFactors = FALSE)))
    }
    
    subset_smallest_values <- function(col) {
      subset_values <- col[valid_subset]
      smallest_indices <- order(subset_values)[1:n]
      smallest_values <- subset_values[smallest_indices]
      smallest_names <- names(subset_values)[smallest_indices]
      evidence <- paste(smallest_names, round(smallest_values, 5), sep = ":", collapse = ", ")
      rank <- round(1 / sum(smallest_values), 5)
      return(list(rank = rank, evidence = evidence))
    }
    
    results <- apply(df, ifelse(match_on == "index", 2, 1), subset_smallest_values)
    ranks <- sapply(results, function(x) x$rank)
    evidences <- sapply(results, function(x) x$evidence)
    
    df_weighted_sorted <- data.frame(
      entity = names(ranks),
      score = ranks,
      evidence = evidences,
      stringsAsFactors = FALSE
    )
    df_weighted_sorted <- df_weighted_sorted[order(df_weighted_sorted$score, decreasing = TRUE), ]
    colnames(df_weighted_sorted)[1] <- key
    return(list(valid_subset = valid_subset, data = df_weighted_sorted))
  }
  
  observeEvent(input$submit, {
    withProgress(message = 'Processing...', value = 0, {
      incProgress(0.1, detail = "Loading data...")
      lipid_protein_df <- read.csv('/home/uchenna/Documents/Rstudio/Omint/lipid_protein_hyperbolicDistance2.csv', row.names = 1, check.names = FALSE)
      protein_metabolite_df <- read.csv('/home/uchenna/Documents/Rstudio/Omint/protein_metabolite_hyperbolicDistance2.csv', row.names = 1, check.names = FALSE)
      metabolite_lipid_df <- read.csv('/home/uchenna/Documents/Rstudio/Omint/metabolite_lipid_hyperbolicDistance2.csv', row.names = 1, check.names = FALSE)
      
      lipids_subset <- process_input(input$lipids_input)
      proteins_subset <- process_input(input$proteins_input)
      metabolites_subset <- process_input(input$metabolites_input)
      
      n <- input$n
      output_data <- list()
      valid_subsets <- list()
      
      if (input$omic_type == "lipid") {
        incProgress(0.2, detail = "Processing lipids...")
        lipid_protein_result <- get_ranked_list(lipids_subset, lipid_protein_df, n, "Protein", match_on = "columns")
        output_data[["Protein (Lipids)"]] <- lipid_protein_result$data
        valid_subsets[["Lipids (Valid for Proteins)"]] <- lipid_protein_result$valid_subset
        
        lipid_metabolite_result <- get_ranked_list(lipids_subset, metabolite_lipid_df, n, "Metabolite", match_on = "columns")
        output_data[["Metabolites (Lipids)"]] <- lipid_metabolite_result$data
        valid_subsets[["Lipids (Valid for Metabolites)"]] <- lipid_metabolite_result$valid_subset
      }
      
      if (input$omic_type == "protein") {
        incProgress(0.4, detail = "Processing proteins...")
        protein_lipid_result <- get_ranked_list(proteins_subset, lipid_protein_df, n, "Lipid", match_on = "index")
        output_data[["Lipids (Proteins)"]] <- protein_lipid_result$data
        valid_subsets[["Proteins (Valid for Lipids)"]] <- protein_lipid_result$valid_subset
        
        protein_metabolite_result <- get_ranked_list(proteins_subset, protein_metabolite_df, n, "Metabolite", match_on = "index")
        output_data[["Metabolites (Proteins)"]] <- protein_metabolite_result$data
        valid_subsets[["Proteins (Valid for Metabolites)"]] <- protein_metabolite_result$valid_subset
      }
      
      if (input$omic_type == "metabolite") {
        incProgress(0.6, detail = "Processing metabolites...")
        metabolite_protein_result <- get_ranked_list(metabolites_subset, protein_metabolite_df, n, "Protein", match_on = "columns")
        output_data[["Protein (Metabolites)"]] <- metabolite_protein_result$data
        valid_subsets[["Metabolites (Valid for Proteins)"]] <- metabolite_protein_result$valid_subset
        
        metabolite_lipid_result <- get_ranked_list(metabolites_subset, metabolite_lipid_df, n, "Lipid", match_on = "index")
        output_data[["Lipids (Metabolites)"]] <- metabolite_lipid_result$data
        valid_subsets[["Metabolites (Valid for Lipids)"]] <- metabolite_lipid_result$valid_subset
      }
      
      incProgress(0.8, detail = "Finalizing...")
      result_data(output_data)
      valid_subsets_data(valid_subsets)
      
      output$downloadButtons <- renderUI({
        if (length(output_data) > 0) {
          downloadButtons <- lapply(names(output_data), function(key) {
            downloadButton(outputId = paste0("download_", key), label = paste("Download", key))
          })
          do.call(tagList, downloadButtons)
        }
      })
      
      lapply(names(output_data), function(key) {
        output[[paste0("download_", key)]] <- downloadHandler(
          filename = function() {
            paste(key, Sys.Date(), ".csv", sep = "_")
          },
          content = function(file) {
            df <- result_data()[[key]]
            df$score <- formatC(df$score, format = "e", digits = 3)
            write.csv(df, file, row.names = FALSE)
          }
        )
      })
      
      output$validSubsetPanel <- renderUI({
        if (length(valid_subsets_data()) > 0) {
          fluidRow(
            column(6, tableOutput("leftValidSubset")),
            column(6, tableOutput("rightValidSubset"))
          )
        }
      })
      
      output$leftValidSubset <- renderTable({
        data.frame(ValidSubset = valid_subsets_data()[[1]])
      })
      
      output$rightValidSubset <- renderTable({
        data.frame(ValidSubset = valid_subsets_data()[[2]])
      })
      
      output$resultsPanel <- renderUI({
        if (length(output_data) > 1) {
          fluidRow(
            column(6, tableOutput("leftPanel")),
            column(6, tableOutput("rightPanel"))
          )
        } else {
          tableOutput("singlePanel")
        }
      })
      
      if (length(output_data) > 1) {
        output$leftPanel <- renderTable({
          df <- output_data[[1]]
          df$score <- formatC(df$score, format = "e", digits = 3)
          df
        })
        output$rightPanel <- renderTable({
          df <- output_data[[2]]
          df$score <- formatC(df$score, format = "e", digits = 3)
          df
        })
      } else if (length(output_data) == 1) {
        output$singlePanel <- renderTable({
          df <- output_data[[1]]
          df$score <- formatC(df$score, format = "e", digits = 3)
          df
        })
      }
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
