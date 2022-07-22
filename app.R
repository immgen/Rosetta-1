### Load in required libraries
library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(shinyWidgets)
library(ComplexHeatmap)
library(ggplot2)
library(ggExtra)
library(circlize)
library(Matrix)
library(shinycssloaders)
library(shinyhelper)
library(magrittr)
library(dplyr)
library(shinyjs)
library(viridis)
library(egg)
library(readxl)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(Seurat)

### Comment this line out when deploying
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))


### Read in the ADT Information table
ADT_table <-
  as.data.frame(read.csv("./adt_hash_seq_file_2022_panel.csv"))

### Master file for experiments
sample_metadata <- read.csv("./sample_metadata.csv")

### Displaying webpage
ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  ### Header
  titlePanel(fluidRow(column(
    4, img(src = "ImmGenIcon.png")
  ),
  column(
    4,
    h1(
      id = "heading",
      "Rosetta",
      align = "center",
      style = 'padding:0px;'
    )
  )),
  windowTitle = "Protein/RNA Mosaic Viewer"),
  tags$head(tags$style(HTML(
    "#heading{color: darkred}"
  )),
  tags$title("Protein/RNA Mosaic Viewer")),
  
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      actionButton(
        'help',
        img(src = "What.png",
            style = "padding-bottom:0px;border-radius: 0px;border-width: 0px")
      ),
      HTML("<br> <br>"),
      selectInput(
        "selected_igt_id",
        "IGT ID",
        choices = unique(sample_metadata$igt_description),
        selected = '1_ImmGen Pilot4 Panel 1'
      ),
      selectInput(
        "selected_igt_samples",
        "Samples",
        choices = NULL,
        multiple = TRUE
      ),
      conditionalPanel(condition = "input.tabs == 'Heatmap'",
                       uiOutput("clusterSelector")),
      textOutput("description"),
      downloadLink('downloadData', 'Ab list is here'),
      width = 2
    ),
    mainPanel(tabsetPanel(
      id = 'tabs',
      tabPanel(
        'Heatmap',
        value = 'Heatmap',
        splitLayout(withSpinner(
          plotOutput(
            "RNA_plot",
            height = 400,
            width = 500,
            brush = brushOpts(
              id = "plot4_brush",
              resetOnNew =
                FALSE,
              delay = 10000,
              delayType =
                "debounce"
            )
          )
        ),
        withSpinner(
          plotOutput(
            "ADT_plot",
            height = 400,
            width = 500,
            brush = brushOpts(
              id = "ADT_brush",
              resetOnNew =
                FALSE,
              delay = 10000,
              delayType =
                "debounce"
            )
          )
        ),
        align = "center"),
        withSpinner(plotOutput("heatmap"))
      ),
      tabPanel(
        'Flow Cytometry Style Plots',
        value = 'facs',
        fluidRow(
          align = "center",
          column(
            id = "first",
            4,
            style = 'padding:10px;',
            withSpinner(
              plotOutput(
                "prot_plot1",
                height = 250,
                width = 250,
                brush = brushOpts(id =
                                    "facs1_brush",
                                  resetOnNew =
                                    TRUE)
              )
            ),
            fixedRow(
              column(
                6,
                align = "center",
                selectizeInput("plot1_X",
                               "X axis:",
                               as.list(c("")),
                               selected =
                                 "TCRB")
              ),
              column(
                6,
                align = "center",
                selectizeInput("plot1_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD19")
              )
            ),
            align = "center"
          ),
          column(
            id = "second",
            4,
            fixedRow(
              column(
                1,
                radioButtons(
                  "highlight_filter",
                  "",
                  choices = c("Highlight", "Filter"),
                  width = "10px"
                )
              ),
              column(
                style = 'padding:10px;',
                3,
                offset = 2,
                withSpinner(
                  plotOutput(
                    "prot_plot2",
                    height = 250,
                    width = 250,
                    brush = brushOpts(id =
                                        "facs2_brush", resetOnNew = FALSE)
                  )
                )
              )
            ),
            fixedRow(
              column(
                4,
                align = "center",
                selectizeInput("plot2_X", "X axis:", as.list(c("")),
                               selected =
                                 "CD4")
              ),
              column(
                4,
                align = "center",
                selectizeInput("plot2_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD8a")
              ),
              column(4, align =
                       "center", selectizeInput("plot2_gate", "Gate 2:", c("A", "none")))
            )
          ),
          column(
            id = "third",
            4,
            fixedRow(
              column(
                1,
                radioButtons(
                  "highlight_filter2",
                  "",
                  choices = c("Highlight", "Filter"),
                  width = "10px"
                )
              ),
              column(
                style = 'padding:10px;',
                3,
                offset = 2,
                withSpinner(plotOutput(
                  "prot_plot3", height = 250, width = 250
                ))
              )
            ),
            fixedRow(
              column(
                4,
                align = "center",
                selectizeInput("plot3_X", "X axis:", as.list(c("")),
                               selected =
                                 "CD44")
              ),
              column(
                4,
                align = "center",
                selectizeInput("plot3_Y", "Y axis:", as.list(c("")),
                               selected =
                                 "CD62L")
              ),
              column(
                4,
                align = "center",
                selectizeInput("plot3_gate", "Gate 3:", c("A", "B", "A+B", "none"), selected =
                                 "A+B")
              )
            ),
            align = "right"
          )
        ),
        column(
          style = 'padding-top:70px;padding-bottom:0px',
          width = 12,
          selectizeInput(
            "UMAP_gate",
            "UMAP gate:",
            c("A", "B", "A+B", "none"),
            selected = "A+B",
            width = "200px"
          ),
          align = "center"
        ),
        tags$style(
          type = 'text/css',
          ".selectize-input { font-size: 12px} .selectize-dropdown { text-align:center !important; font-size: 12px } .container-fluid {  min-width: 1800px; }",
          HTML(
            "
      #first {
          border: 1px solid black;height=350px;width=350px;
      }
      #second {
          border: 1px solid black;height=350px;width=350px;
      }
      #third {
          border: 1px solid black;height=350px;width=350px;
      }
    "
          )
        ),
        withSpinner(plotOutput("rnaumap"))
      )
    ))
  )
)

server <- function(input, output, session) {
  observeEvent(input$help, {
    showModal(
      modalDialog(
        align = "center",
        title = "Help",
        HTML(
          "This is an interactive tool to help compare single cell Protein and RNA landscapes.
      <br> <br>
      On the Heatmap tab (default), select cells on either the Protein or RNA UMAP.
           The web tool will output the most differential proteins and most differential genes for your selected
           cell population as two heatmaps. It will also highlight the cells you selected on the other unselected UMAP.
           The outputed heatmaps show the same cell orders between the two heatmaps. <br> <br> Note: It may take
           10-30 seconds for the sample list to appear for subsetting when a dataset is selected.
           <br> <br> To clear,
           click anywhere on the plot you highlighted. <br> <br>"
        ),
        img(
          src = "HeatmapStep1.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br>
      Select the flow cytometry style plots tab to plot protein markers on XY axes. Each subsequent plot will output
      the selected cells according to the selected proteins for the XY axes. Gates can be selected and
      RNA and ADT umaps will be outputed after each plot depending on the highlighting or filtering used. <br> <br>
      On the flow cytometry style tab, select 2 protein markers for the first XY axis and select the population
      you are interested in.<br> <br>"
        ),
        img(
          src = "FlowStep1.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br> The selected populations will show on the plot plotted as the X axis and Y axis. <br> <br>"
        ),
        img(
          src = "Gatechange.png",
          height = "25%",
          width = "25%"
        ),
        HTML(
          "<br> <br> Next, select the cells of interest in the second plot in order to plot in the third and final plot.
      This final plot will show the cells plotted on the XY axis selected."
        ),
        img(
          src = "Flowrevisited.png",
          height = "50%",
          width = "50%"
        ),
        HTML(
          "<br> <br> To change what cells are shown in the umap (either Gate A, Gate B, Gate A and B,
      or all cells in dataset), select the appropriate option under UMAP gate."
        )
      )
    )
  })
  
  
  
  #according to the dataset selected, update the adt protein list
  observeEvent(input$selected_igt_id, {
    updateSelectizeInput(
      session,
      inputId = "plot1_X",
      label = "X axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "TCRB"
    )
    updateSelectizeInput(
      session,
      inputId = "plot1_Y",
      label = "Y axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "TCRGD"
    )
    updateSelectizeInput(
      session,
      inputId = "plot2_X",
      label = "X axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "CD4"
    )
    updateSelectizeInput(
      session,
      inputId = "plot2_Y",
      label = "Y axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "CD8A"
    )
    updateSelectizeInput(
      session,
      inputId = "plot3_X",
      label = "X axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "CD44"
    )
    updateSelectizeInput(
      session,
      inputId = "plot3_Y",
      label = "Y axis",
      choices = as.list(ADT_table$Protein_Symbol),
      selected = "CD62L"
    )
  })
  
  
  #for the link to download the antibody list used for the experiment
  output$downloadData <- downloadHandler(
    filename <- function() {
    "adt_hash_seq_file_2022_panel.csv"
  },
    content <- function(file) {
      write.csv(ADT_table, file)
  })
  
  #depending on the dataset, select the clusters needed to show on the main screen
  output$clusterSelector <- renderUI({
    selectizeInput("cluster", "Cluster color", as.list(c("seurat_clusters", "hash.ID", 'sample_name')),
                   selected = "sample_name")
    
  })
  
  output$description <- renderText({
    description <- read.csv("descriptions.csv")
    desc <-
      description$Description[description$Dataset == input$selected_igt_id]
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  #OBSERVE EVENT ----
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  # -------------------------------------------------------------------
  # Linked plots (middle and right)
  facs1 <- reactiveValues(x = NULL, y = NULL)
  ranges4 <- reactiveValues(x = NULL, y = NULL)
  ADTrange <- reactiveValues(x = NULL, y = NULL)
  facs2 <- reactiveValues(x = NULL, y = NULL)
  
  selectedIGT <- reactive({
    filter(sample_metadata, igt_id == strsplit(input$selected_igt_id, "_")[[1]][1])
    
  })
  
  observeEvent(selectedIGT(), {
    choices = unique(selectedIGT()[, 2])
    
    updateSelectInput(inputId = "selected_igt_samples", choices = choices)
  })
  
  
  datasetInput <- reactive({
    sc <-
      readRDS(paste0("./SCRNA", strsplit(input$selected_igt_id, "_")[[1]][1], "/dataset.Rds"))
    
    samples_to_subset <- input$selected_igt_samples
    
    Idents(sc) <- sc$sample_name
    sc <- subset(sc, idents = samples_to_subset)
    
    adtumap1 <- sc@reductions$adt.umap@cell.embeddings[, 1]
    adtumap2 <- sc@reductions$adt.umap@cell.embeddings[, 2]
    rnaumap1 <- sc@reductions$rnaumap@cell.embeddings[, 1]
    rnaumap2 <- sc@reductions$rnaumap@cell.embeddings[, 2]
    
    data <- as.data.frame(t(as.matrix(sc@assays$ADT@counts)))
    data <- log2(data + 1)
    data <- data[, order(colnames(data))]
    
    rna_data <- sc@assays$RNA@data
    return(list(
      data,
      adtumap1,
      adtumap2,
      rnaumap1,
      rnaumap2,
      sc@meta.data,
      rna_data
    ))
  })
  
  datasetInput_d <- datasetInput %>% debounce(100)
  
  output$prot_plot1 <- renderPlot({
    
    data = datasetInput_d()

    data = data[[1]]
    
    ggplot(data, aes(data[, input$plot1_X], data[, input$plot1_Y])) + geom_point(col =
                                                                                   rgb(0, 0, 0, 0.1)) +
      xlab(input$plot1_X) + ggtitle(paste0(input$plot1_X, " vs ", input$plot1_Y, " (select gate A)")) +
      ylab(input$plot1_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
  })
  
  output$prot_plot2 <- renderPlot({
    data <- datasetInput_d()
    data = data[[1]]
    
    if (!is.null(facs1$x)) {
      if (input$plot2_gate == "none") {
        data %>% ggplot(aes(data[, input$plot2_X], data[, input$plot2_Y])) + geom_point(col =
                                                                                          rgb(0, 0, 0, 0.1)) +
          xlab(input$plot2_X) + ggtitle(paste0(input$plot2_X, " vs ", input$plot2_Y, " (select gate B)")) +
          ylab(input$plot2_Y) +
          xlim(0, 10) + ylim(0, 10) + theme_bw()
      } else{
        if (input$plot2_gate == "A") {
          data2 <-
            data[which(data[, input$plot1_X] < facs1$x[2] &
                         data[, input$plot1_X] > facs1$x[1]),]
          data2 <-
            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                          data2[, input$plot1_Y] > facs1$y[1]),]
          
          if (input$highlight_filter == "Filter") {
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) +
              geom_point(data = data2,
                         aes(data2[, input$plot2_X], data2[, input$plot2_Y]),
                         col = rgb(1, 0, 0, 0.1)) +
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(data =
                                                                                                 data,
                                                                                               aes(data[, input$plot2_X], data[, input$plot2_Y]),
                                                                                               col = "grey89") +
              geom_point(data = data2,
                         aes(data2[, input$plot2_X], data2[, input$plot2_Y]),
                         col = rgb(1, 0, 0, 0.1)) +
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
          }
        } else{
          if (input$plot2_gate == "B" & !is.null(facs2$x)) {
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[, input$plot2_Y] > facs2$y[1]),]
            
            data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(col =
                                                                                                 rgb(0, 0, 0, 0.1)) +
              #geom_point(data=data2,aes(data2[,input$ADTmarker],data2[,input$ADTmarker2[2]]),col="red")+
              xlab(input$plot2_X) + ggtitle(paste0(
                input$plot2_X,
                " vs ",
                input$plot2_Y,
                " (select gate B)"
              )) + ylab(input$plot2_Y) +
              xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            if (input$plot2_gate == "A+B" & !is.null(facs2$x)) {
              data2 <-
                data[which(data[, input$plot1_X] < facs1$x[2] &
                             data[, input$plot1_X] > facs1$x[1]),]
              data2 <-
                data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                              data2[, input$plot1_Y] > facs1$y[1]),]
              
              data2 <-
                data2[which(data2[, input$plot2_X] < facs2$x[2] &
                              data2[, input$plot2_X] > facs2$x[1]),]
              data2 <-
                data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                              data2[, input$plot2_Y] > facs2$y[1]),]
              
              data2 %>% ggplot(aes(data2[, input$plot2_X], data2[, input$plot2_Y])) + geom_point(col =
                                                                                                   rgb(0, 0, 0, 0.1)) +
                #geom_point(data=data2,aes(data2[,input$ADTmarker],data2[,input$ADTmarker2[2]]),col="red")+
                xlab(input$plot2_X) + ggtitle(paste0(
                  input$plot2_X,
                  " vs ",
                  input$plot2_Y,
                  " (select gate B)"
                )) + ylab(input$plot2_Y) +
                xlim(0, 10) + ylim(0, 10) + theme_bw()
            }
          }
        }
      }
    } else{
      if (input$plot2_X == "") {
        ggplot() + geom_point() + theme_bw()
      } else{
        ggplot(data, aes(data[, input$plot2_X], data[, input$plot2_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot2_X) + ggtitle(paste0(input$plot2_X, " vs ", input$plot2_Y, " (select gate B)")) +
          ylab(input$plot2_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      }
    }
    
  })
  
  output$prot_plot3 <- renderPlot({
    data <- datasetInput_d() 
    data = data[[1]]
    
    if (!is.null(facs2$x)) {
      if (input$plot3_gate == "none") {
        ggplot(data, aes(data[, input$plot3_X], data[, input$plot3_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
          ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      } else{
        if (input$plot3_gate == "A") {
          #subset data to gate A
          data2 <-
            data[which(data[, input$plot1_X] < facs1$x[2] &
                         data[, input$plot1_X] > facs1$x[1]),]
          data2 <-
            data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                          data2[input$plot1_Y] > facs1$y[1]),]
          if (input$highlight_filter2 == "Filter") {
            ggplot() +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
              ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          } else{
            ggplot() + geom_point(data = data,
                                  aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                  col = "grey89") +
              geom_point(data = data[rownames(data2),],
                         aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                         col = rgb(0, 0, 1, 0.1)) +
              xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
              ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
          }
        } else{
          if (input$plot3_gate == "B" && input$plot2_gate != "A") {
            #subset data to gate B
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[input$plot2_Y] > facs2$y[1]),]
            if (input$highlight_filter2 == "Filter") {
              ggplot() +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            } else{
              ggplot() +  geom_point(data = data,
                                     aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                     col = "grey89") +
                geom_point(data = data[rownames(data2),],
                           aes(data[rownames(data2), input$plot3_X], data[rownames(data2), input$plot3_Y]),
                           col = rgb(0, 0, 1, 0.1)) +
                xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
            }
          } else{
            if ((input$plot3_gate == "A+B") ||
                (input$plot3_gate == "B" &&
                 input$plot2_gate == "A")) {
              if (!is.null(facs1$x)) {
                data2 <-
                  data[which(data[, input$plot1_X] < facs1$x[2] &
                               data[, input$plot1_X] > facs1$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                data2[input$plot1_Y] > facs1$y[1]),]
                
                data3 <-
                  data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                data2[, input$plot2_X] > facs2$x[1]),]
                data3 <-
                  data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                data3[input$plot2_Y] > facs2$y[1]),]
              } else{
                data3 <-
                  data[which(data[, input$plot2_X] < facs2$x[2] &
                               data[, input$plot2_X] > facs2$x[1]),]
                data3 <-
                  data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                data3[input$plot2_Y] > facs2$y[1]),]
              }
              
              if (input$highlight_filter2 == "Filter") {
                ggplot() +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot3_X],
                                 data[rownames(data3), input$plot3_Y]),
                             col = rgb(1, 0, 1, 0.1)) +
                  geom_point(col = rgb(1, 0, 1, 0.5)) +
                  xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                  ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) +
                  theme_bw()
              } else{
                ggplot() + geom_point(data = data,
                                      aes(data[, input$plot3_X], data[, input$plot3_Y]),
                                      col = "grey89") +
                  geom_point(data = data[rownames(data3),],
                             aes(data[rownames(data3), input$plot3_X],
                                 data[rownames(data3), input$plot3_Y]),
                             col = rgb(1, 0, 1, 0.1)) +
                  geom_point(col = rgb(1, 0, 1, 0.5)) +
                  xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
                  ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) +
                  theme_bw()
              }
            }
          }
        }
      }
    } else{
      if (input$plot3_X == "") {
        ggplot() + geom_point() + theme_bw()
      } else{
        ggplot(data, aes(data[, input$plot3_X], data[, input$plot3_Y])) + geom_point(col =
                                                                                       rgb(0, 0, 0, 0.1)) +
          xlab(input$plot3_X) + ggtitle(paste0(input$plot3_X, " vs ", input$plot3_Y)) +
          ylab(input$plot3_Y) + xlim(0, 10) + ylim(0, 10) + theme_bw()
      }
    }
    
    
  })
  
  output$rnaumap <- renderPlot({
    data_list <- datasetInput_d() 
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    
    
    if (input$tabs == "facs") {
      if (is.null(facs1$x) & is.null(facs2$x)) {
        g1 <-
          ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                  "grey89") +
          ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                               20)) + removeGrid()
        g2 <-
          ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                  "grey89") +
          ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                   20)) + removeGrid()
        g2 + g1
      } else{
        if (xor(!is.null(facs1$x),!is.null(facs2$x))) {
          if (is.null(facs1$x)) {
            data2 <-
              data[which(data[, input$plot2_X] < facs2$x[2] &
                           data[, input$plot2_X] > facs2$x[1]),]
            data2 <-
              data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                            data2[, input$plot2_Y] > facs2$y[1]),]
          } else{
            data2 <-
              data[which(data[, input$plot1_X] < facs1$x[2] &
                           data[, input$plot1_X] > facs1$x[1]),]
            data2 <-
              data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                            data2[, input$plot1_Y] > facs1$y[1]),]
          }
          
          g1 <-
            ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                    "grey89") +
            geom_point(
              data = umaps[rownames(data2),],
              aes(rnaumap_1, rnaumap_2),
              col = rgb(1, 0, 0),
              size = 0.5
            ) +
            ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                 20)) + removeGrid()
          g2 <-
            ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                    "grey89") +
            geom_point(
              data = umaps[rownames(data2),],
              aes(adtumap_1, adtumap_2),
              col = rgb(1, 0, 0),
              size = 0.5
            ) +
            ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                     20)) + removeGrid()
          g2 + g1
        } else{
          if (!is.null(facs1$x) & !is.null(facs2$x)) {
            data2 <-
              data[which(data[, input$plot1_X] < facs1$x[2] &
                           data[, input$plot1_X] > facs1$x[1]),]
            data2 <-
              data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                            data2[, input$plot1_Y] > facs1$y[1]),]
            
            data3 <-
              data2[which(data2[, input$plot2_X] < facs2$x[2] &
                            data2[, input$plot2_X] > facs2$x[1]),]
            data3 <-
              data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                            data3[input$plot2_Y] > facs2$y[1]),]
            if (input$UMAP_gate == "none") {
              g1 <-
                ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                        "grey89") +
                ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                     20)) + removeGrid()
              g2 <-
                ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                        "grey89") +
                ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                         20)) + removeGrid()
              g2 + g1
            } else{
              if (input$UMAP_gate == "A") {
                data2 <-
                  data[which(data()[, input$plot1_X] < facs1$x[2] &
                               data[, input$plot1_X] > facs1$x[1]),]
                data2 <-
                  data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                data2[input$plot1_Y] > facs1$y[1]),]
                
                g1 <-
                  ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  geom_point(
                    data = umaps[rownames(data2),],
                    aes(rnaumap_1, rnaumap_2),
                    col = rgb(0, 0, 1),
                    size = 0.5
                  ) +
                  ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                       20)) + removeGrid()
                g2 <-
                  ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                          "grey89") +
                  geom_point(
                    data = umaps[rownames(data2),],
                    aes(adtumap_1, adtumap_2),
                    col = rgb(0, 0, 1),
                    size = 0.5
                  ) +
                  ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                           20)) + removeGrid()
                g2 + g1
              } else{
                if (input$UMAP_gate == "B") {
                  data2 <-
                    data[which(data[, input$plot2_X] < facs2$x[2] &
                                 data[, input$plot2_X] > facs2$x[1]),]
                  data2 <-
                    data2[which(data2[, input$plot2_Y] < facs2$y[2] &
                                  data2[input$plot2_Y] > facs2$y[1]),]
                  
                  g1 <-
                    ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(rnaumap_1, rnaumap_2),
                      col = rgb(1, 0, 0),
                      size = 0.5
                    ) +
                    ggtitle("RNA UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                         20)) + removeGrid()
                  g2 <-
                    ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                            "grey89") +
                    geom_point(
                      data = umaps[rownames(data2),],
                      aes(adtumap_1, adtumap_2),
                      col = rgb(1, 0, 0),
                      size = 0.5
                    ) +
                    ggtitle("Protein UMAP") + theme_bw() + theme(plot.title = element_text(hjust = 0.5, size =
                                                                                             20)) + removeGrid()
                  g2 + g1
                } else{
                  if (input$UMAP_gate == "A+B") {
                    if (!is.null(facs1$x)) {
                      data2 <-
                        data[which(data[, input$plot1_X] < facs1$x[2] &
                                     data[, input$plot1_X] > facs1$x[1]),]
                      data2 <-
                        data2[which(data2[, input$plot1_Y] < facs1$y[2] &
                                      data2[input$plot1_Y] > facs1$y[1]),]
                      
                      data3 <-
                        data2[which(data2[, input$plot2_X] < facs2$x[2] &
                                      data2[, input$plot2_X] > facs2$x[1]),]
                      data3 <-
                        data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                      data3[input$plot2_Y] > facs2$y[1]),]
                    } else{
                      data3 <-
                        data[which(data[, input$plot2_X] < facs2$x[2] &
                                     data[, input$plot2_X] > facs2$x[1]),]
                      data3 <-
                        data3[which(data3[, input$plot2_Y] < facs2$y[2] &
                                      data3[input$plot2_Y] > facs2$y[1]),]
                    }
                    
                    g1 <-
                      ggplot(umaps, aes(rnaumap_1, rnaumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data3),],
                        aes(rnaumap_1, rnaumap_2),
                        col = rgb(1, 0, 1, 0.5),
                        size = 0.5
                      ) +
                      theme_bw() + removeGrid() + ggtitle("RNA UMAP") +
                      theme(plot.title = element_text(hjust = 0.5, size = 20))
                    g2 <-
                      ggplot(umaps, aes(adtumap_1, adtumap_2)) + geom_point(size = 0.5, col =
                                                                              "grey89") +
                      geom_point(
                        data = umaps[rownames(data3),],
                        aes(adtumap_1, adtumap_2),
                        col = rgb(1, 0, 1, 0.5),
                        size = 0.5
                      ) +
                      theme_bw() + removeGrid() + ggtitle("Protein UMAP") +
                      theme(plot.title = element_text(hjust = 0.5, size = 20))
                    g2 + g1
                  }
                }
              }
            }
          }
        }
      }
    }
    
  })
  
  autoInvalidate <- reactiveTimer(2000)
  
  output$RNA_plot <- renderPlot({
    data_list <- datasetInput_d()
    
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    metadata = data_list[[6]]
    
    
    umaps2 <-
      umaps[which(
        umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
          umaps$rnaumap_2 < ranges4$y[2] &
          umaps$rnaumap_2 > ranges4$y[1]
      ),]
    
    cluster <- input$cluster
    umaps$cluster <- metadata[, cluster]
    if (is.null(input$ADT_brush) & !is.null(input$plot4_brush)) {
      umaps %>% ggplot(aes(rnaumap_1, rnaumap_2, col = cluster)) + geom_point(size =
                                                                         0.5) +
        labs(colour = cluster) +
        ggtitle("RNA UMAP (select cells to show top differential heatmaps)") +
        theme_bw() + removeGrid()
    } else{
      if (!is.null(input$ADT_brush) & is.null(input$plot4_brush)) {
        umaps %>% ggplot(aes(rnaumap_1, rnaumap_2)) + geom_point(col = "grey89") +
          geom_point(
            data = umaps2,
            aes(rnaumap_1, rnaumap_2),
            col = "red",
            size = 0.5
          ) + labs(colour = cluster) +
          ggtitle("RNA UMAP (select cells to show top differential heatmaps)") +
          theme_bw() + removeGrid()
      } else{
        if (is.null(input$ADT_brush) & is.null(input$plot4_brush)) {
          umaps %>% ggplot(aes(rnaumap_1, rnaumap_2, col = cluster)) +
            geom_point(size = 0.5) +
            labs(colour = cluster) +
            ggtitle("RNA UMAP (select cells to show top differential heatmaps)") +
            theme_bw() + removeGrid()
        }
      }
    }
  })
  
  output$ADT_plot <- renderPlot({
    data_list <- datasetInput_d() 
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    metadata = data_list[[6]]
    
    
    umaps2 <-
      umaps[which(
        umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
          umaps$rnaumap_2 < ranges4$y[2] &
          umaps$rnaumap_2 > ranges4$y[1]
      ),]
    
    cluster <- input$cluster
    umaps$cluster <- metadata[, cluster]
    
    if (!is.null(input$plot4_brush)) {
      umaps %>% ggplot(aes(adtumap_1, adtumap_2)) + geom_point(col = "grey89") +
        ggtitle("Protein UMAP (select cells to show top differential heatmaps)") +
        geom_point(
          data = umaps2,
          aes(adtumap_1, adtumap_2),
          col = "red",
          size = 0.5
        ) +
        labs(colour = cluster) +
        theme_bw() + removeGrid()
    } else{
      umaps %>% ggplot(aes(adtumap_1, adtumap_2, col = cluster)) + geom_point(size =
                                                                                0.5) +
        ggtitle("Protein UMAP (select cells to show top differential heatmaps)") +
        geom_point(
          data = umaps2,
          aes(adtumap_1, adtumap_2),
          col = "red",
          size = 0.5
        ) +
        labs(colour = cluster) +
        theme_bw() + removeGrid()
    }
  })
  
  output$heatmap <- renderPlot({
    data_list <- datasetInput_d() 
    data = data_list[[1]]
    umaps = data.frame(
      adtumap_1 = data_list[[2]],
      adtumap_2 = data_list[[3]],
      rnaumap_1 = data_list[[4]],
      rnaumap_2 = data_list[[5]]
    )
    metadata = data_list[[6]]
    
    rna_data = data_list[[7]]
    
    cluster <- input$cluster
    umaps$cluster <- metadata[, cluster]
    
    if (!is.null(input$plot4_brush) & is.null(input$ADT_brush)) {
      umaps2 <-
        umaps[which(
          umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
            umaps$rnaumap_2 < ranges4$y[2] &
            umaps$rnaumap_2 > ranges4$y[1]
        ),]
      data2 <-
        data[which(
          umaps$rnaumap_1 < ranges4$x[2]  & umaps$rnaumap_1 > ranges4$x[1] &
            umaps$rnaumap_2 < ranges4$y[2] &
            umaps$rnaumap_2 > ranges4$y[1]
        ),]
    } else{
      if (!is.null(input$ADT_brush) & is.null(input$plot4_brush)) {
        umaps2 <-
          umaps[which(
            umaps$adtumap_1 < ADTrange$x[2]  & umaps$adtumap_1 > ADTrange$x[1] &
              umaps$adtumap_2 < ADTrange$y[2] &
              umaps$adtumap_2 > ADTrange$y[1]
          ),]
        data2 <-
          data[which(
            umaps$adtumap_1 < ADTrange$x[2]  & umaps$adtumap_1 > ADTrange$x[1] &
              umaps$adtumap_2 < ADTrange$y[2] &
              umaps$adtumap_2 > ADTrange$y[1]
          ),]
        
      }
    }
    
    if (!is.null(input$plot4_brush) |
        !is.null(input$ADT_brush)) {
      FC_ADT <-
        apply(data2, 2, mean) / apply(data[-which(rownames(data) %in% rownames(data2)),], 2, mean)
      FC_ADT <- names(FC_ADT[order(FC_ADT, decreasing = T)])[1:30]
      
      print(rna_data)
      FC_RNA <-
        rowMeans(rna_data[, rownames(data2)]) / rowMeans(rna_data[,-which(colnames(rna_data) %in%
                                                                            rownames(data2))])
      FC_RNA <-
        FC_RNA[names(which((
          rowSums(rna_data > 0) / dim(data2)[2] > 0.05 * dim(data2)[2]
        ) == T))]
      FC_RNA <- names(FC_RNA[order(FC_RNA, decreasing = T)])[1:30]
      FC_RNA <- rowMeans(rna_data[FC_RNA,])
      FC_RNA <- names(FC_RNA)
      
      if (dim(data2)[1] > 2000) {
        s <- sample(1:dim(data2)[1], 2000, replace = F)
        data2 <- data2[s,]
      }
      rnadata <-
        as.data.frame(as.matrix(rna_data[rownames(rna_data) %in% FC_RNA,
                                         rownames(data2)]))
      
      
      column_order = names(colSums(data2[, FC_ADT])[order(colSums(data2[, FC_ADT]), decreasing =
                                                            T)])
      data2 <- t(data2[, column_order])
      h <-
        Heatmap(
          as.matrix(data2),
          cluster_rows = F,
          cluster_columns = T,
          show_row_names = T,
          row_names_gp = gpar(fontsize = 15),
          heatmap_legend_param = list(title = ""),
          show_column_dend = F,
          show_column_names = F,
          column_title = "Protein Heatmap (top protein markers differentiating selected cells from rest of cells)"
        )
      
      ht <- draw(h)
      h2 <-
        Heatmap(
          as.matrix(rnadata[, column_order(ht)]),
          cluster_rows = T,
          cluster_columns = F,
          row_names_gp = gpar(fontsize = 15),
          heatmap_legend_param = list(title = ""),
          show_column_names = F,
          column_title = "RNA Heatmap (top gene markers differentiating selected cells from rest of cells)"
        )
      
      draw(h2 + h, auto_adjust = F)
    }
  })
  
  output$sample_table <- renderTable({
    data_list <- datasetInput_d() 
    metadata = data_list[[6]] 
    
    samples <- factor(metadata$sample_name)
    samples
  })
  
  # OBSERVE ----
  
  observe({
    brush_facs1 <- input$facs1_brush
    if (!is.null(brush_facs1)) {
      facs1$x <- c(brush_facs1$xmin, brush_facs1$xmax)
      facs1$y <- c(brush_facs1$ymin, brush_facs1$ymax)
      
    } else {
      facs1$x <- NULL
      facs1$y <- NULL
    }
  })
  
  observe({
    brush <- input$plot4_brush
    if (!is.null(brush)) {
      ranges4$x <- c(brush$xmin, brush$xmax)
      ranges4$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges4$x <- NULL
      ranges4$y <- NULL
    }
  })
  
  observe({
    brush_ADT <- input$ADT_brush
    if (!is.null(brush_ADT)) {
      ADTrange$x <- c(brush_ADT$xmin, brush_ADT$xmax)
      ADTrange$y <- c(brush_ADT$ymin, brush_ADT$ymax)
      
    } else {
      ADTrange$x <- NULL
      ADTrange$y <- NULL
    }
  })
  
  observe({
    brush_facs2 <- input$facs2_brush
    if (!is.null(brush_facs2)) {
      facs2$x <- c(brush_facs2$xmin, brush_facs2$xmax)
      facs2$y <- c(brush_facs2$ymin, brush_facs2$ymax)
      
    } else {
      facs2$x <- NULL
      facs2$y <- NULL
    }
  })
  
}

shinyApp(ui, server)
