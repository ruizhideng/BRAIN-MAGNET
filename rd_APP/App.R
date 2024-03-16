library(rsconnect)
library(dplyr)
library(readxl)
# library(gridExtra)
library(shiny)
library(ggplot2)
library(ggseqlogo)
# library(cowplot)
library(patchwork)
library(DT)
library(tidyverse)
# library(reticulate)
# np <- import("numpy")



# set working path
setwd("/Users/Ruizhi/Work/EMC/Projects/Deep_learning/rd_APP/")
# load enhancer activity which corresponding to the contri score order
enhancer_activity <- read.table("data/Enhancer_activity.txt", sep = "\t", header = TRUE) %>%
  unite("tmp", Chr:Start, sep = ":") %>%
  unite("location", tmp:End, sep = "-")

# load NSC data
# nsc_hypoth_scores <- np$load("data/shap_explanations_NSC.npy")
# nsc_inps <- np$load("data/inp_NSC.npy")
# nsc_contri_scores <- np$multiply(nsc_hypoth_scores, nsc_inps)
# 
# rm(nsc_hypoth_scores, nsc_inps)
# gc()
# 
# saveRDS(nsc_contri_scores,"data/nsc_contri_scores.rds")
nsc_contri_scores <- readRDS("data/nsc_contri_scores.rds")


# load ESC data
# esc_hypoth_scores <- np$load("data/shap_explanations_ESC.npy")
# esc_inps <- np$load("data/inp_ESC.npy")
# esc_contri_scores <- np$multiply(esc_hypoth_scores, esc_inps)
# 
# rm(esc_hypoth_scores, esc_inps)
# gc()
# 
# saveRDS(esc_contri_scores,"data/esc_contri_scores.rds")
esc_contri_scores <- readRDS("data/esc_contri_scores.rds")


# Define UI 
ui <- fluidPage(
  titlePanel("ChIP-STARR-seq Visualization"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput(
        inputId = "cell_type",
        label = "Cell Type:",
        choices = c("NSC", "ESC"),
        selected = "NSC"
      ),
      textInput(
        inputId = "enhancer_ID",
        label = "Enhancer ID",
        value = "chr1:154862304-154862970",
        placeholder = "For example, chr1:1000-2000"
      ),
      actionButton("plot_button", "Plot"),
      downloadButton("download_plot", "Download Plot")
    ),
    
    mainPanel(
      fluidRow(
        column(width = 12, align = "center",
               h4("The plot of contribution scores", style = "font-weight: bold;"),
               p("choose the enhancer regions from the following table"),
               tags$div(
                 id = "plot_container",
                 style = "max-height: 50vh;",
                 plotOutput(outputId = "contri_plot",height = "250px")
               )
        )
      ),
      fluidRow(
        column(width = 12,
               h4("Enhancer Activity"),
               DT::dataTableOutput("enhancer_table")  # Use DT::dataTableOutput for DataTable
        )
      )
    )
  )
)

# Define server 
server <- function(input, output, session) {
  
  selected_data <- reactiveVal()  # Create a reactiveVal to store the selected data
  plot_data <- reactiveVal(NULL)  # Create a reactiveVal to store the plot data
  
  # Observe changes in the data type input and set the appropriate data
  observe({
    if (input$cell_type == "NSC") {
      selected_data(nsc_contri_scores)
    } else {
      selected_data(esc_contri_scores)
    }
  })
  
  # Only generate the plot when the "Plot" button is clicked
  observeEvent(input$plot_button, {
    req(input$enhancer_ID)  # Ensure that the input is not empty
    
    plots <- list()
    
    # Skip empty lines
    if (input$enhancer_ID == "") {
      next
    }
    
    enhancer_order <- as.numeric(row.names(enhancer_activity)[which(enhancer_activity$location == input$enhancer_ID)])
    
    examp <- selected_data()[enhancer_order, 1:4, 1:1000]
    rownames(examp) <- c('A', 'C', 'G', 'T')
    
    p <- ggseqlogo(examp, method='custom', seq_type='dna', font='roboto_bold') + 
      scale_x_continuous(breaks = seq(0, 1000, by = 100), expand = c(0, 0), limits = c(0, 1005)) + 
      ggtitle(input$enhancer_ID) +
      theme_bw() +
      theme(panel.grid=element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.2),
            axis.text=element_text(size=3,colour = "black"),
            axis.ticks = element_line(linewidth = 0.2),
            legend.position="none",
            plot.title = element_text(color="black", size=4))
    
    # Adjust the size of the ggseqlogo plot
    # p <- p + plot_layout(widths = unit(9, "cm"), heights = unit(1.5, "cm"))
    
    # Store the plot data in the reactiveVal
    plot_data(p)
  })
  
  # Render the plot
  output$contri_plot <- renderPlot({
    plot_data()  # Get the plot data from the reactiveVal
  }, res = 300)
  
  # Define a function to generate the PDF plot and save it to a file
  generatePDFPlot <- function(plot, filename) {
    pdf(filename, width = 9, height = 4)
    print(plot)
    dev.off()
  }
  
  # Define a download handler for the plot
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("plot_", gsub("[:/]", "_", input$enhancer_ID), ".pdf", sep = "")
    },
    content = function(file) {
      plot <- plot_data()  # Get the plot data from the reactiveVal
      generatePDFPlot(plot, file)  # Generate the PDF plot and save it to the file
    }
  )
  
  # Render the enhancer_activity table using DT::renderDT
  output$enhancer_table <- DT::renderDT({
    enhancer_activity
  })
}

# Create a Shiny app object 
shinyApp(ui = ui, server = server)
