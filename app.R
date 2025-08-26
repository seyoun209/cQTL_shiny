# Check and install required packages
required_packages <- c("shiny", "plotly", "dplyr", "tidyr", "shinyjs", "viridis")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load libraries
library(shiny)
library(plotly)
library(dplyr)
library(tidyr)
library(shinyjs)
library(viridis)

# Load the data once
load("data/counts_annot.Rdata")

# Filter for significant genes to improve performance
sig_genes <- counts_annot %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  pull(ensgID) %>%
  unique()

# Get unique gene symbols for dropdown and handle NAs
gene_options <- counts_annot %>%
  filter(ensgID %in% sig_genes) %>%
  distinct(SYMBOL, ensgID) %>%
  # Replace NA SYMBOL values with their Ensembl IDs
  mutate(SYMBOL = ifelse(is.na(SYMBOL), ensgID, SYMBOL)) %>%
  arrange(SYMBOL)

# UI definition
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Gene Expression Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      # Server-side selectize for better performance
      selectizeInput("gene", "Select Gene:", 
                    choices = NULL, # Will be populated server-side
                    options = list(
                      placeholder = 'Type to search for a gene',
                      maxOptions = 100
                    )),
      
      radioButtons("plotType", "Plot Type:",
                 choices = list("Boxplot" = "box", "Scatter" = "scatter"),
                 selected = "box"),
      
      selectInput("colorBy", "Color By:",
                choices = list("Sex" = "sex", 
                              "Condition" = "condition",
                              "Age" = "age"),
                selected = "sex")
    ),
    
    mainPanel(
      plotlyOutput("expressionPlot", height = "600px"),
      
      # Gene info panel
      wellPanel(
        h4("Gene Information"),
        verbatimTextOutput("geneInfo")
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Initialize the gene selector with server-side options
  updateSelectizeInput(session, 'gene', 
                      choices = setNames(gene_options$ensgID, gene_options$SYMBOL),
                      server = TRUE)
  
  # Reactive data subset for selected gene
  selected_data <- reactive({
    req(input$gene)
    counts_annot %>%
      filter(ensgID == input$gene) %>%
      mutate(log2_count = log2(count + 1))
  })
  
  # Gene info output
  output$geneInfo <- renderPrint({
    req(nrow(selected_data()) > 0)
    gene_data <- selected_data() %>% 
      select(SYMBOL, ensgID, seqnames, start, end, log2FoldChange, padj) %>%
      distinct() %>%
      head(1)
    
    cat("Symbol:", gene_data$SYMBOL %||% "Not available", "\n")
    cat("Ensembl ID:", gene_data$ensgID, "\n")
    cat("Location:", gene_data$seqnames, ":", gene_data$start, "-", gene_data$end, "\n")
    cat("Log2 Fold Change:", round(gene_data$log2FoldChange, 3), "\n")
    cat("Adjusted p-value:", format(gene_data$padj, scientific = TRUE), "\n")
    cat("Sample size:", nrow(selected_data()), "observations\n")
  })
  
  # Expression plot
  output$expressionPlot <- renderPlotly({
    req(nrow(selected_data()) > 0)
    
    gene_data <- selected_data()
    color_var <- input$colorBy
    
    # Custom colors to match JavaScript version
    my_colors <- list(
      condition = c("CTL" = "#9FCCE4", "FNF" = "#FAB394"),
      sex = c("M" = "#0075B0", "F" = "#F8B7CD"),
      age = viridis::plasma(10) # Age will use a gradient
    )
    
    if(input$plotType == "box") {
      if(color_var == "sex") {
        # Create a combined factor for condition+sex to get 4 separate boxplots
        gene_data$condition_sex <- paste(gene_data$condition, gene_data$sex, sep = "-")
        
        # Set the order of the combined factor
        gene_data$condition_sex <- factor(gene_data$condition_sex, 
                                         levels = c("CTL-M", "CTL-F", "FNF-M", "FNF-F"))
        
        # Create a new color palette for the combined groups
        combined_colors <- c(
          "CTL-M" = "#0075B0",  # Male CTL - blue
          "CTL-F" = "#F8B7CD",  # Female CTL - pink
          "FNF-M" = "#005587",  # Male FNF - darker blue
          "FNF-F" = "#E47A9E"   # Female FNF - darker pink
        )
        
        p <- ggplot(gene_data, aes(x = condition_sex, y = log2_count, fill = condition_sex)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          scale_fill_manual(values = combined_colors) +
          geom_jitter(aes(color = condition_sex), width = 0.2, size = 3, alpha = 0.7) +
          scale_color_manual(values = combined_colors) +
          labs(x = "Condition-Sex", fill = "Group", color = "Group")
        
      } else if(color_var == "age") {
        # For age, just show scatter plot with age on x-axis
        p <- ggplot(gene_data, aes(x = age, y = log2_count, color = condition)) +
          geom_point(size = 3, alpha = 0.7) +
          scale_color_manual(values = my_colors$condition) +
          geom_smooth(aes(group = condition, color = condition), method = "lm", se = FALSE) +
          labs(x = "Age")
        
      } else {
        # Default boxplot by condition
        p <- ggplot(gene_data, aes(x = condition, y = log2_count, group = condition)) +
          geom_boxplot(aes(fill = condition), alpha = 0.7, outlier.shape = NA) +
          scale_fill_manual(values = my_colors$condition) +
          geom_jitter(aes(color = condition), width = 0.2, size = 3, alpha = 0.7) +
          scale_color_manual(values = my_colors$condition)
      }
    } else {
      # For scatter plot, always show age vs expression
      p <- ggplot(gene_data, aes(x = age, y = log2_count))
      
      # Color by selected variable
      if(color_var == "sex") {
        # Create a combined factor for condition+sex
        gene_data$condition_sex <- paste(gene_data$condition, gene_data$sex, sep = "-")
        gene_data$condition_sex <- factor(gene_data$condition_sex, 
                                        levels = c("CTL-M", "CTL-F", "FNF-M", "FNF-F"))
        
        # Combined colors
        combined_colors <- c(
          "CTL-M" = "#0075B0",  # Male CTL - blue
          "CTL-F" = "#F8B7CD",  # Female CTL - pink
          "FNF-M" = "#005587",  # Male FNF - darker blue
          "FNF-F" = "#E47A9E"   # Female FNF - darker pink
        )
        
        p <- p + geom_point(aes(color = condition_sex), size = 3, alpha = 0.7) +
               scale_color_manual(values = combined_colors) +
               geom_smooth(aes(group = condition_sex, color = condition_sex), method = "lm", se = FALSE) +
               labs(color = "Group")
        
      } else if(color_var == "condition") {
        p <- p + geom_point(aes(color = condition), size = 3, alpha = 0.7) +
               scale_color_manual(values = my_colors$condition) +
               geom_smooth(aes(group = condition, color = condition), method = "lm", se = FALSE)
               
      } else {
        # If coloring by age, group by condition for trend lines
        p <- p + geom_point(aes(color = age), size = 3, alpha = 0.7) +
               scale_color_viridis_c(option = "plasma") +
               geom_smooth(aes(group = condition, color = condition), method = "lm", se = FALSE) +
               scale_color_manual(values = my_colors$condition, guide = "none") + 
               # Add a second legend for condition trend lines
               annotate("segment", x = min(gene_data$age), 
                       xend = min(gene_data$age) + 5, 
                       y = max(gene_data$log2_count), 
                       yend = max(gene_data$log2_count), 
                       color = my_colors$condition["CTL"], linetype = 1) +
               annotate("text", x = min(gene_data$age) + 7, 
                       y = max(gene_data$log2_count), 
                       label = "CTL trend", hjust = 0) +
               annotate("segment", x = min(gene_data$age), 
                       xend = min(gene_data$age) + 5, 
                       y = max(gene_data$log2_count) - 0.2, 
                       yend = max(gene_data$log2_count) - 0.2, 
                       color = my_colors$condition["FNF"], linetype = 1) +
               annotate("text", x = min(gene_data$age) + 7, 
                       y = max(gene_data$log2_count) - 0.2, 
                       label = "FNF trend", hjust = 0)
      }
    }
    
    # Common styling for both plots
    p <- p + labs(
      title = paste("Expression of", unique(gene_data$SYMBOL)[1] %||% input$gene),
      y = "log2(count + 1)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
    
    ggplotly(p, tooltip = c("y", "x", "color", "fill", "donor"))
  })
}

# Helper function to handle NULL values
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x

# Run the application
shinyApp(ui = ui, server = server)