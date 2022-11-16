library(shiny)
library(shinythemes)
library(Seurat)
library(uwot)
library(DT)
library(ggplot2)
#library(plotly)
library(scCustomize)
library(shinyalert)
library(shinydisconnect)
library(shinycssloaders)
#Read .csv files
#data <-read.csv("data/table_shiny.csv",dec = ",")
datos <-suppressWarnings(read.csv("data/ICLN98_HumanMeta.csv",dec=","))
data <-read.csv("data/BM_gene.csv",dec = ",")
dat <-read.csv("data/SP_gene.csv",dec = ",")
da <-read.csv("data/TH_gene.csv",dec = ",")
d <-read.csv("data/LN_gene.csv",dec = ",")

# Read the .rds seurat file
set_5 <- readRDS("data/ImmuneTissuesAllWithNewReduction.rds")


#pulling out data from the integrated object
Idents(set_5) <- "MainCluster"
set_5$MainCluster<-factor(set_5$MainCluster, levels=c("SP_25","BM_9","TH_39","BM_20","LN_44","BM_31","SP_6","BM_17","LN_4","LN_21","LN_1","LN_35","TH_9","LN_18","TH_20","TH_42","LN_42","BM_37","SP_24","BM_4","SP_26","TH_41","SP_27","BM_6","BM_32","SP_16","SP_19","SP_11","SP_2","BM_18","SP_18","TH_31","LN_33","SP_7","BM_22","SP_5","BM_27","TH_27","LN_41","SP_14","SP_10","BM_21","SP_12","TH_21","SP_9","LN_14","SP_17","TH_8","LN_9","LN_25","LN_28","LN_13","TH_26","TH_16","LN_29","TH_43","LN_40","SP_13","BM_14","BM_1","SP_20","LN_15","LN_6","LN_20","TH_30","LN_8","BM_8","SP_3","SP_21","BM_34","SP_22","TH_36","LN_32","LN_19","LN_10","TH_1","LN_43","BM_2","TH_25","SP_8","LN_2","TH_2","BM_25","LN_3","LN_26","LN_24","LN_11","LN_5","LN_17","LN_12","TH_38","BM_16","TH_14","BM_5","TH_28","TH_18","BM_24","BM_15","BM_10","BM_26","BM_30","TH_33","TH_4","BM_28","TH_35","LN_39","BM_39","LN_30","LN_7","BM_3","BM_35","TH_11","BM_19","TH_22","TH_7","TH_5","BM_33","TH_17","TH_3","LN_27","BM_7","TH_6","LN_36","SP_4","SP_1","TH_37","LN_37","BM_23","SP_23","LN_38","LN_31","BM_29","TH_15","SP_15","TH_40","BM_38","TH_10","BM_11","LN_16","LN_22","BM_12","TH_12","TH_24","TH_23","TH_29","TH_19","LN_34","BM_36","TH_34","LN_23","BM_13","TH_13","TH_32"))

lev<-c("SP_25","BM_9","TH_39","BM_20","LN_44","BM_31","SP_6","BM_17","LN_4","LN_21","LN_1","LN_35","TH_9","LN_18","TH_20","TH_42","LN_42","BM_37","SP_24","BM_4","SP_26","TH_41","SP_27","BM_6","BM_32","SP_16","SP_19","SP_11","SP_2","BM_18","SP_18","TH_31","LN_33","SP_7","BM_22","SP_5","BM_27","TH_27","LN_41","SP_14","SP_10","BM_21","SP_12","TH_21","SP_9","LN_14","SP_17","TH_8","LN_9","LN_25","LN_28","LN_13","TH_26","TH_16","LN_29","TH_43","LN_40","SP_13","BM_14","BM_1","SP_20","LN_15","LN_6","LN_20","TH_30","LN_8","BM_8","SP_3","SP_21","BM_34","SP_22","TH_36","LN_32","LN_19","LN_10","TH_1","LN_43","BM_2","TH_25","SP_8","LN_2","TH_2","BM_25","LN_3","LN_26","LN_24","LN_11","LN_5","LN_17","LN_12","TH_38","BM_16","TH_14","BM_5","TH_28","TH_18","BM_24","BM_15","BM_10","BM_26","BM_30","TH_33","TH_4","BM_28","TH_35","LN_39","BM_39","LN_30","LN_7","BM_3","BM_35","TH_11","BM_19","TH_22","TH_7","TH_5","BM_33","TH_17","TH_3","LN_27","BM_7","TH_6","LN_36","SP_4","SP_1","TH_37","LN_37","BM_23","SP_23","LN_38","LN_31","BM_29","TH_15","SP_15","TH_40","BM_38","TH_10","BM_11","LN_16","LN_22","BM_12","TH_12","TH_24","TH_23","TH_29","TH_19","LN_34","BM_36","TH_34","LN_23","BM_13","TH_13","TH_32")
levels(set_5)<-lev

data_list = list(set_5=set_5)

# ui.R

ui <- fluidPage(theme = shinytheme("cerulean"),
                
                
                titlePanel(h1("Visualization of sc-RNAseq Data")),
                navbarPage(" ",
                           tabPanel(icon("home"),
                                    
                                    fluidRow(column(tags$img(src="Pig.png",width="250px",height="200px"),width=2),
                                             column(
                                               
                                               br(),
                                               p("
This project aims at understanding pig immune system for food production and translation research. This will provide an immune cell atlas as a basis for future research.
Moreover, it will improve cell type and tissues specific gene expression data for genetic selection.
", 
                                                 style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                                               br(),
                                               
                                               p("
                                                 Immune tissues were collected from two 6 month old healthy pigs.
                                                 Created clusters of single cell data and proved they are unique and distinguishable. Identified gene expression patterns and markers for different immune cell types.
                                                 Identified tissue specific vs. peripheral immune cell types by comparing against porcine PBMCs
                                                 Identified tissue-specific differences between porcine and human cell types.
                                                 Used canonical markers,"
                                                 , a(href="https://www.frontiersin.org/articles/10.3389/fgene.2021.689406/full", "porcine PBMC data",target="_blank") ,
                                                 "and human tissue-specific data to annotate the porcine immune cell atlas."
                                                 ,style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                                               
                                               width=8),
                                             column(
                                               br(),
                                               tags$img(src="FAANG_logo_RGBc.png",width="200px",height="130px"),
                                               br(),
                                               br(),
                                               p("For more information please check the",em("FAANG official website"),"page clicking",
                                                 br(),
                                                 a(href="https://www.faang.org", "Here",target="_blank"),style="text-align:center;color:black"),
                                               
                                               width=2)),
                                    useShinyalert(),
                                    
                                    
                                    hr(),
                                    tags$style(".fa-database {color:#e87722}"),
                                    h3(p(em("Gene Expression for Single Cell Immune Tissues"),icon("database",lib = "font-awesome"),style="color:black;text-align:center")),
                                    tabsetPanel(
                                      tabPanel("Immune_Tissues_Meta",
                                               fluidRow(column(DT::dataTableOutput("Immune_Tissues_Meta"),
                                                               width = 12))),
                                    tabPanel("BM_Gene_List",
                                             fluidRow(column(DT::dataTableOutput("BM_Gene_List"),
                                                             width = 12))),
                                    tabPanel("SP_Gene_List",
                                             fluidRow(column(DT::dataTableOutput("SP_Gene_List"),
                                                             width = 12))),
                                    tabPanel("TH_Gene_List",
                                             fluidRow(column(DT::dataTableOutput("TH_Gene_List"),
                                                             width = 12))),
                           
                                    tabPanel("LN_Gene_List",
                                            fluidRow(column(DT::dataTableOutput("LN_Gene_List"),
                                                             width = 12)))),
                                    
                                    
                                    
                                    hr(),
                                    p(em("Developed by"),br("Tuggle Lab"),style="text-align:center; font-family: times")
                           ),
                           tabPanel("GENES",icon = icon("user"), 
                                    tags$style(".glyphicon-signal {color:#E87722}"),
                                    h3(p(em("VISUALIZATION OF GENES"),icon("signal", lib="glyphicon"),style="color:black;text-align:center")),
                                    br(),
                                    
                                    column(width = 4,
                                           textInput("gene2", label = "Gene Symbol/Ensembl ID", value = "CD3E")),
                                    submitButton("Update Now", icon=(""),useShinyalert()),
                                    
                                    column(width = 12,
                                           withSpinner(plotOutput("genePlot2",height=700,width="1200px"))),
                                    
                                    
                                    
                                    
                                    
                                    column(width = 12,
                                           selectInput("dataset5", label = h3("Gene Expression among Clusters"),
                                                       choices = list("Integrated"="set_5"),
                                                       selected = "set_5"),
                                           
                                           withSpinner(plotOutput("vPlot5",height=700,width="2900px"
                                                                  
                                           ))),
                                    column(width = 12,
                                           selectInput("dataset5", label = h3("ALL TISSUES COMBINED"),
                                                       choices = list("Integrated"="set_5"),
                                                       selected = "set_5"),
                                           
                                           
                                           #radioButtons("prev","The Integrated Seurat object of Four Tissues:", c("Preview","No Preview"),selected="Preview", inline=TRUE),
                                           radioButtons(inputId = "run","The Integrated Seurat object of Four Tissues:", c("Preview","No Preview"),selected="Preview", inline=TRUE),
                                           
                                           hr(),
                                           p("The Integrated Seurat object of Four Tissues:", style = "color:#888888;"),
                                           
                                           withSpinner(plotOutput("genePlot6",height=700,width="1200px",
                                                                  dblclick="plot_1_dblclick",
                                                                  brush=brushOpts(
                                                                    id="plot1_brush",
                                                                    resetOnNew = TRUE
                                                                  ))),
                                           
                                           textOutput("select")
                                           
                                           
                                    )
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                           )
                           
                           
                           
                )
)



#server.R


server<- function(input,output, session){
  
  
  output$Immune_Tissues_Meta <- DT::renderDataTable(
    DT::datatable({
      datos
    },
    options=list(lengthMenu=list(c(5,15,20),c('5','15','20')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c("Tissues","Number of cells post QC","Number of Features","Number of Clusters")
    ))
  output$BM_Gene_List <- DT::renderDataTable(
    DT::datatable({
      data
    },
    options=list(lengthMenu=list(c(100,400,800),c('100','400','800')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c("Cluster","Gene")
    ))
  output$SP_Gene_List <- DT::renderDataTable(
    DT::datatable({
      dat
    },
    options=list(lengthMenu=list(c(100,400,800),c('100','400','800')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c("Cluster","Gene")
    ))
  output$TH_Gene_List <- DT::renderDataTable(
    DT::datatable({
      da
    },
    options=list(lengthMenu=list(c(100,400,800),c('100','400','800')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c("Cluster","Gene")
    ))
  output$LN_Gene_List <- DT::renderDataTable(
    DT::datatable({
      d
    },
    options=list(lengthMenu=list(c(100,400,800),c('100','400','800')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c("Cluster","Gene")
    ))
  
  
  
  shinyalert(
    title = "Help",
    text = "For visualization of individual gene expression for each cluster and tissue, as well as the integrated clusters, please use the “genes” webpage above ",
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
 

  
  
  
  dataset5Input <- reactive({
    infile <- data_list[[input$dataset5]]
    
  })
  
  
  #output plot
  
  output$genePlot2 <- renderPlot({
    Sys.sleep(5)
    
    shinyalert(
      title = "Note",
      text = "This page might take 20 seconds to load.",
      size = "m", 
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = FALSE,
      type = "success",
      showConfirmButton = TRUE,
      showCancelButton = FALSE,
      confirmButtonText = "OK",
      confirmButtonCol = "#AEDEF4",
      timer = 0,
      imageUrl = "",
      animation = TRUE
    )
    FeaturePlot_scCustom(dataset5Input(),split.by="Tissue",num_columns = 2,reduction="Tissue_UMAP",features=toupper(input$gene2))
    
    
  })

  
  #  output$genePlot6<-renderPlot({
  #Sys.sleep(5)
  #   switch(input$prev, 
  #                 "Preview"=FeaturePlot_scCustom(dataset5Input(),reduction="umap",features= (input$gene2)),
  #          "No Preview"= "nothing"
  #                  )
  #   Sys.sleep(10)
  #  })
  
  observe({
    observeEvent(input$run,{
      if (input$run=="Preview"){
        output$genePlot6<- renderPlot({
          FeaturePlot_scCustom(dataset5Input(),reduction="umap",features= toupper(input$gene2))
        })
      }
      else if (input$run=="No Preview"){
        output$select<- renderText({
          #input$gene2
          "No Preview Selected"
          
        })
      }
      
    })
  })
  
  
  #  output$genePlot6 <- renderPlot({
  #    Sys.sleep(10)
  #   switch(input$preview,
  #    "preview"=FeaturePlot_scCustom(dataset5Input(),reduction="umap",features= (input$gene3))
  #    )
  #  })
  
  
  output$vPlot5 <- renderPlot({
    Sys.sleep(5)
    VlnPlot(dataset5Input(),split.by="MainCluster",features =toupper (input$gene2),split.plot=TRUE)
    
    
    
  })
  
}

shinyApp(ui, server)





