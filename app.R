#app.R
library(shiny)
library(tidyverse)
library(shinythemes)
library(DT)

js <- c(
  "table.on('key', function(e, datatable, key, cell, originalEvent){",
  "  var targetName = originalEvent.target.localName;",
  "  if(key == 13 && targetName == 'body'){",
  "    $(cell.node()).trigger('dblclick.dt');",
  "  }",
  "});",
  "table.on('keydown', function(e){",
  "  var keys = [9,13,37,38,39,40];",
  "  if(e.target.localName == 'input' && keys.indexOf(e.keyCode) > -1){",
  "    $(e.target).trigger('blur');",
  "  }",
  "});",
  "table.on('key-focus', function(e, datatable, cell, originalEvent){",
  "  var targetName = originalEvent.target.localName;",
  "  var type = originalEvent.type;",
  "  if(type == 'keydown' && targetName == 'input'){",
  "    if([9,37,38,39,40].indexOf(originalEvent.keyCode) > -1){",
  "      $(cell.node()).trigger('dblclick.dt');",
  "    }",
  "  }",
  "});"
)
# Define UI
ui <- fluidPage(
  titlePanel("Variant Analysis Report"),
  br(),
  
  
    sidebarLayout(
      sidebarPanel(
        width = 5,
        selectInput(inputId = "Gene", label = "Select a gene", choices = c("ALK","BRAF","BRCA2","CDKN2A","DDR2","EGFR","ERBB2","FGFR1","KRAS","MET","NRAS","PDGFRA","PIK3CA","ROS1","TP53")),
        selectInput(inputId = "Patient", label = "Select a Patient", choices = c("DNA001r1","DNA001r2","DNA002r1","DNA002r2","DNA003r1","DNA003r2","DNA004r1","DNA004r2","DNA005r1","DNA005r2","DNA006r1","DNA007r1","DNA008r1","DNA009r1","DNA010r1","DNA011r1","DNA012r1","DNA013r1","DNA014r1","DNA015r1","DNA016r1","DNA017r1","DNA018r1","DNA019r1","DNA020r1","DNA022r1","DNA023r1","DNA027r1","DNA029r1","DNA030r1","DNA032r1","DNA033r1","DNA034r1","DNA035r1","DNA036r1","DNA037r1","DNA038r1","DNA041r1","DNA042r1","DNA043r1","DNA045r1","DNA047r1","DNA048r1","DNA049r1","DNA051r1","DNA052r1","DNA053r1","DNA054r1","DNA055r1","DNA056r1","DNA057r1","DNA059r1","DNA061r1","DNA062r1","DNA063r1","DNA064r1","DNA065r1","DNA066r1","DNA068r1","DNA069r1","DNA070r1","DNA071r1","DNA072r1","DNA074r1","DNA075r1","DNA076r1","DNA077r1","DNA078r1","DNA085r1","DNA086r1","DNA087r1","DNA088r1","DNA089r1","DNA092r1","DNA093r1")),
        tags$style(type='text/css', ".selectize-input { font-size: 24px; line-height: 24px;} .selectize-dropdown { font-size: 20px; line-height: 20px; }"),
        
        titlePanel(windowTitle = T, "Sample information:"),
        DTOutput("Metadata"),
        
        titlePanel(windowTitle = T, "Sequencing data Quality Control:"),
        
        tags$li(uiOutput("ui_open_tab1")),
        tags$li(uiOutput("ui_open_tab2"))
        
      ),
      mainPanel(
        titlePanel(windowTitle = T, "Lung Cancer Targeted Sequencing Panel Results"),
        DTOutput("table")
      )
    ),
  br(),
  titlePanel("*supported by TUBITAK TEYDEB 7200066")

  )

# Define server function
server <- function(input, output) {
  
  snpfile <- reactive({
    paste0("data/", input$Patient, "_snpEff_genes.txt")
    })
  
  output$table <- renderDT({
    
    snpeff_data <- read_delim(snpfile(), delim = "\t")
    
    datatable(
      snpeff_data[which(snpeff_data$GeneName == input$Gene), c(rep(TRUE, 4), colSums(snpeff_data[which(snpeff_data$GeneName == input$Gene), -c(1:4)])!=0)],
      selection = "none",
      editable = TRUE, 
      callback = JS(js),
      extensions = "KeyTable",
      options = list(
        keys = TRUE
      )
    )
  })
  
  output$ui_open_tab1 <- renderUI({
    a("FastQC Results - Read 1", target="_blank", href=paste0(input$Patient, "_R1_fastqc.html"))
  })
  output$ui_open_tab2 <- renderUI({
    a("FastQC Results - Read 2", target="_blank", href=paste0(input$Patient, "_R2_fastqc.html"))
  })
  
  sampleSheet <- readxl::read_xlsx("data/sampleSheet.xlsx",col_names = T)
  
  output$Metadata <- renderDT({
    sampleSheet %>% filter(SampleID == input$Patient) %>% datatable(options = list(dom = 't'))
  })
  
  
  }

# Create Shiny object
shinyApp(ui = ui, server = server)


