library(shiny)
library(dplyr)

# Load mapping data
protein_tissueLocation_map <- read.csv('/home/uchenna/Documents/python/protein_tissueLocation_clean.csv')


# Define UI for the application
ui <- fluidPage(
  titlePanel("Omint: Omics Integration"),
  p("Get a ranked list of other omics based on its connection to one omic list."),
  
  sidebarLayout(
    sidebarPanel(
      # New multi-select dropdown for selecting tissue locations, sorted alphabetically
      selectInput(
        inputId = "tissue_location",
        label = "Select Tissue Location:",
        choices = sort(unique(protein_tissueLocation_map$tissue_location)),
        selected = sort(unique(protein_tissueLocation_map$tissue_location))[1],
        multiple = TRUE
      ),
      
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
      uiOutput("validSubsetPanel")  # output to display the valid_subset
    )
  )
)



# Define server logic required to rank entities
server <- function(input, output) {
  #Create variables with reactive values 
  result_data <- reactiveVal(NULL)
  valid_subsets_data <- reactiveVal(NULL)  # Store the valid subsets
  
  # Function to clean, de-duplicate, and split the user input
  process_input <- function(input_string) {
    unique(trimws(unlist(strsplit(input_string, ","))))
  }
  
  # Function to get ranked list of entities from a given subset
  get_ranked_list <- function(subset, df, n, key, match_on = "index") {
    
    #Extract valid subset. That is, subsets in map
    if (match_on == "index") {
      valid_subset <- rownames(df)
    } else {
      valid_subset <- colnames(df)
    }
    
    if (length(valid_subset) == 0) {
      return(list(valid_subset = valid_subset, data = data.frame(entity = character(), rank = numeric(), evidence = character(), stringsAsFactors = FALSE)))
    }
    
    #Function to return the rank and evidence based on n least distances among all subset distance-values
    subset_smallest_values <- function(col) {
      subset_values <- col[valid_subset]
      smallest_indices <- order(subset_values)[1:n]
      smallest_values <- subset_values[smallest_indices]
      smallest_names <- names(subset_values)[smallest_indices]
      evidence <- paste(smallest_names, round(smallest_values, 5), sep = ":", collapse = ", ")
      rank <- round(1 / sum(smallest_values), 5)
      return(list(rank = rank, evidence = evidence))
    }
    
    #Process dataframe df by applying subset_smallest_values function on the row or col
    results <- apply(df, ifelse(match_on == "index", 2, 1), subset_smallest_values)
    ranks <- sapply(results, function(x) x$rank) #Extracts the rank component of each element of results
    evidences <- sapply(results, function(x) x$evidence) #Extracts the evidence component of each element of results
    
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
      
      # Load the data
      data <- read.csv('/home/uchenna/Documents/python/nodes_multi_omics_category2.csv')
      #lipids_subset <- process_input(input$lipids_input)
      #proteins_subset <- process_input(input$proteins_input)
      #metabolites_subset <- process_input(input$metabolites_input)
      
      # Separate the data by categories
      if (input$omic_type == "protein") {
        proteins_subset <- process_input(input$proteins_input)
        
        proteins1 <- data %>% filter(category == "protein")
        proteins <- proteins1 %>% filter(id %in% proteins_subset) #Extract only intersections of proteins_subset and proteins
        if (!is.null(input$tissue_location) && length(input$tissue_location) > 0) { #if the user selected tissue locations
          #further filter the protein to include only the ones in the selected tissue locations
          # First, extract the proteins from the mapping data corresponding to the selected tissue locations
          selected_proteins <- protein_tissueLocation_map %>%
            filter(tissue_location %in% input$tissue_location) %>%
            pull(protein)
          
          # Then, filter proteins1 to keep only those whose 'id' appears in the selected proteins list
          proteins <- proteins %>%
            filter(id %in% selected_proteins)
        }  
        lipids <- data %>% filter(category == "lipid")
        metabolites <- data %>% filter(category == "metabolite")
        
        # Initialize distance matrices
        lipid_protein_df <- matrix(NA, nrow = nrow(proteins), ncol = nrow(lipids),
                                   dimnames = list(proteins$id, lipids$id))
        protein_metabolite_df <- matrix(NA, nrow = nrow(proteins), ncol = nrow(metabolites),
                                        dimnames = list(proteins$id, metabolites$id))
      }
      
      if (input$omic_type == "lipid") {
        lipids_subset <- process_input(input$lipids_input)
        
        lipids1 <- data %>% filter(category == "lipid")
        lipids <- lipids1 %>% filter(id %in% lipids_subset) #Extract only intersections of lipids_subset and lipids
        
        proteins <- data %>% filter(category == "protein")
        if (!is.null(input$tissue_location) && length(input$tissue_location) > 0) { #if the user selected tissue locations
          #further filter the protein to include only the ones in the selected tissue locations
          # First, extract the proteins from the mapping data corresponding to the selected tissue locations
          selected_proteins <- protein_tissueLocation_map %>%
            filter(tissue_location %in% input$tissue_location) %>%
            pull(protein)
          
          # Then, filter proteins1 to keep only those whose 'id' appears in the selected proteins list
          proteins <- proteins %>%
            filter(id %in% selected_proteins)
        }  
        metabolites <- data %>% filter(category == "metabolite")
        
        # Initialize distance matrices
        lipid_protein_df <- matrix(NA, nrow = nrow(proteins), ncol = nrow(lipids),
                                   dimnames = list(proteins$id, lipids$id))
        metabolite_lipid_df <- matrix(NA, nrow = nrow(metabolites), ncol = nrow(lipids),
                                      dimnames = list(metabolites$id, lipids$id))
      }
      
      if (input$omic_type == "metabolite") {
        metabolites_subset <- process_input(input$metabolites_input)
        
        metabolites1 <- data %>% filter(category == "metabolite")
        metabolites <- metabolites1 %>% filter(id %in% metabolites_subset) #Extract only intersections of metabolites_subset and metabolites
        
        proteins <- data %>% filter(category == "protein")
        if (!is.null(input$tissue_location) && length(input$tissue_location) > 0) { #if the user selected tissue locations
          #further filter the protein to include only the ones in the selected tissue locations
          # First, extract the proteins from the mapping data corresponding to the selected tissue locations
          selected_proteins <- protein_tissueLocation_map %>%
            filter(tissue_location %in% input$tissue_location) %>%
            pull(protein)
          
          # Then, filter proteins1 to keep only those whose 'id' appears in the selected proteins list
          proteins <- proteins %>%
            filter(id %in% selected_proteins)
        }  
        
        lipids <- data %>% filter(category == "lipid")
        
        # Initialize distance matrices
        metabolite_lipid_df <- matrix(NA, nrow = nrow(metabolites), ncol = nrow(lipids),
                                      dimnames = list(metabolites$id, lipids$id))
        protein_metabolite_df <- matrix(NA, nrow = nrow(proteins), ncol = nrow(metabolites),
                                        dimnames = list(proteins$id, metabolites$id))
      }
      
      # Function to calculate hyperbolic distance
      hyperbolic_distance <- function(ZI_r, ZI_theta, ZJ_r, ZJ_theta) {
        delta <- pi - abs(pi - abs(ZI_theta - ZJ_theta))
        d <- cosh(ZI_r) * cosh(ZJ_r) - sinh(ZI_r) * sinh(ZJ_r) * cos(delta)
        d[d < 1] <- 1
        d <- acosh(d)
        d[ZI_r == ZJ_r & ZI_theta == ZJ_theta] <- 0
        return(d)
      }
      #compute lipid_protein_df and metabolite_lipid_df if input subset are lipids
      if (input$omic_type == "lipid") {
        # Compute hyperbolic distances between lipids and proteins
        for (i in 1:nrow(proteins)) {
          for (j in 1:nrow(lipids)) {
            lipid_protein_df[i, j] <- hyperbolic_distance(
              ZI_r = proteins$r[i], ZI_theta = proteins$theta[i],
              ZJ_r = lipids$r[j], ZJ_theta = lipids$theta[j]
            )
          }
        }
        
        # Compute hyperbolic distances between metabolites and lipids
        for (i in 1:nrow(metabolites)) {
          for (j in 1:nrow(lipids)) {
            metabolite_lipid_df[i, j] <- hyperbolic_distance(
              ZI_r = metabolites$r[i], ZI_theta = metabolites$theta[i],
              ZJ_r = lipids$r[j], ZJ_theta = lipids$theta[j]
            )
          }
        }
      lipid_protein_df <- as.data.frame(lipid_protein_df, stringsAsFactors = FALSE)
      metabolite_lipid_df <- as.data.frame(metabolite_lipid_df, stringsAsFactors = FALSE)
      }
      
      #compute lipid_protein_df and protein_metabolite_df if input subset are proteins
      if (input$omic_type == "protein") {
        # Compute hyperbolic distances between lipids and proteins
        for (i in 1:nrow(proteins)) {
          for (j in 1:nrow(lipids)) {
            lipid_protein_df[i, j] <- hyperbolic_distance(
              ZI_r = proteins$r[i], ZI_theta = proteins$theta[i],
              ZJ_r = lipids$r[j], ZJ_theta = lipids$theta[j]
            )
          }
        }
        
        # Compute hyperbolic distances between proteins and metabolites
        for (i in 1:nrow(proteins)) {
          for (j in 1:nrow(metabolites)) {
            protein_metabolite_df[i, j] <- hyperbolic_distance(
              ZI_r = proteins$r[i], ZI_theta = proteins$theta[i],
              ZJ_r = metabolites$r[j], ZJ_theta = metabolites$theta[j]
            )
          }
        }
      lipid_protein_df <- as.data.frame(lipid_protein_df, stringsAsFactors = FALSE)
      protein_metabolite_df <- as.data.frame(protein_metabolite_df, stringsAsFactors = FALSE)  
      }  
      #compute protein_metabolite_df and metabolite_lipid_df if input subset are metabolites
      if (input$omic_type == "metabolite") {
        # Compute hyperbolic distances between metabolites and lipids
        for (i in 1:nrow(metabolites)) {
          for (j in 1:nrow(lipids)) {
            metabolite_lipid_df[i, j] <- hyperbolic_distance(
              ZI_r = metabolites$r[i], ZI_theta = metabolites$theta[i],
              ZJ_r = lipids$r[j], ZJ_theta = lipids$theta[j]
            )
          }
        }
        # Compute hyperbolic distances between proteins and metabolites
        for (i in 1:nrow(proteins)) {
          for (j in 1:nrow(metabolites)) {
            protein_metabolite_df[i, j] <- hyperbolic_distance(
              ZI_r = proteins$r[i], ZI_theta = proteins$theta[i],
              ZJ_r = metabolites$r[j], ZJ_theta = metabolites$theta[j]
            )
          }
        }
      metabolite_lipid_df <- as.data.frame(metabolite_lipid_df, stringsAsFactors = FALSE)
      protein_metabolite_df <- as.data.frame(protein_metabolite_df, stringsAsFactors = FALSE)  
      }
      # Convert matrices to data frames for easier handling
      #lipid_protein_df <- as.data.frame(lipid_protein_df, stringsAsFactors = FALSE)
      #metabolite_lipid_df <- as.data.frame(metabolite_lipid_df, stringsAsFactors = FALSE)
      #protein_metabolite_df <- as.data.frame(protein_metabolite_df, stringsAsFactors = FALSE)  
      
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

      
