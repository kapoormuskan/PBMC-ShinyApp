library(shiny)
library(shinythemes)
library(Seurat)
library(uwot)
library(DT)
library(ggplot2)
library(plotly)
library(scCustomize)

datos <-read.csv("data/ICLN98_HumanMeta.csv",dec = ",")
data <-read.csv("data/ICLN00_HumanMeta.csv",dec = ",")
dat <-read.csv("data/THY00_HumanMeta.csv",dec = ",")
da <-read.csv("data/THY98_HumanMeta.csv",dec = ",")

set_1 <- readRDS("data/pbmc_10k_v3.rds")
set_2 <- readRDS("data/Spleen_Combined.rds")
set_3 <- readRDS("data/pbmc_10k_v3.rds")

set_4 <- readRDS("data/Spleen_Combined.rds")
set_5 <- readRDS("data/ImmuneTissuesAllWithNewReduction.rds")
Idents(set_5) <- "MainCluster"
#set_5$MainCluster<-factor(set_5$MainCluster, levels=c("SP_25","BM_9","TH_39","BM_20","LN_44","BM_31","SP_6","BM_17","LN_4","LN_21","LN_1","LN_35","TH_9","LN_18","TH_20","TH_42","LN_42","BM_37","SP_24","BM_4","SP_26","TH_41","SP_27","BM_6","BM_32","SP_16","SP_19","SP_11","SP_2","BM_18","SP_18","TH_31","LN_33","SP_7","BM_22","SP_5","BM_27","TH_27","LN_41","SP_14","SP_10","BM_21","SP_12","TH_21","SP_9","LN_14","SP_17","TH_8","LN_9","LN_25","LN_28","LN_13","TH_26","TH_16","LN_29","TH_43","LN_40","SP_13","BM_14","BM_1","SP_20","LN_15","LN_6","LN_20","TH_30","LN_8","BM_8","SP_3","SP_21","BM_34","SP_22","TH_36","LN_32","LN_19","LN_10","TH_1","LN_43","BM_2","TH_25","SP_8","LN_2","TH_2","BM_25","LN_3","LN_26","LN_24","LN_11","LN_5","LN_17","LN_12","TH_38","BM_16","TH_14","BM_5","TH_28","TH_18","BM_24","BM_15","BM_10","BM_26","BM_30","TH_33","TH_4","BM_28","TH_35","LN_39","BM_39","LN_30","LN_7","BM_3","BM_35","TH_11","BM_19","TH_22","TH_7","TH_5","BM_33","TH_17","TH_3","LN_27","BM_7","TH_6","LN_36","SP_4","SP_1","TH_37","LN_37","BM_23","SP_23","LN_38","LN_31","BM_29","TH_15","SP_15","TH_40","BM_38","TH_10","BM_11","LN_16","LN_22","BM_12","TH_12","TH_24","TH_23","TH_29","TH_19","LN_34","BM_36","TH_34","LN_23","BM_13","TH_13","TH_32"))

#lev<-c("SP_25","BM_9","TH_39","BM_20","LN_44","BM_31","SP_6","BM_17","LN_4","LN_21","LN_1","LN_35","TH_9","LN_18","TH_20","TH_42","LN_42","BM_37","SP_24","BM_4","SP_26","TH_41","SP_27","BM_6","BM_32","SP_16","SP_19","SP_11","SP_2","BM_18","SP_18","TH_31","LN_33","SP_7","BM_22","SP_5","BM_27","TH_27","LN_41","SP_14","SP_10","BM_21","SP_12","TH_21","SP_9","LN_14","SP_17","TH_8","LN_9","LN_25","LN_28","LN_13","TH_26","TH_16","LN_29","TH_43","LN_40","SP_13","BM_14","BM_1","SP_20","LN_15","LN_6","LN_20","TH_30","LN_8","BM_8","SP_3","SP_21","BM_34","SP_22","TH_36","LN_32","LN_19","LN_10","TH_1","LN_43","BM_2","TH_25","SP_8","LN_2","TH_2","BM_25","LN_3","LN_26","LN_24","LN_11","LN_5","LN_17","LN_12","TH_38","BM_16","TH_14","BM_5","TH_28","TH_18","BM_24","BM_15","BM_10","BM_26","BM_30","TH_33","TH_4","BM_28","TH_35","LN_39","BM_39","LN_30","LN_7","BM_3","BM_35","TH_11","BM_19","TH_22","TH_7","TH_5","BM_33","TH_17","TH_3","LN_27","BM_7","TH_6","LN_36","SP_4","SP_1","TH_37","LN_37","BM_23","SP_23","LN_38","LN_31","BM_29","TH_15","SP_15","TH_40","BM_38","TH_10","BM_11","LN_16","LN_22","BM_12","TH_12","TH_24","TH_23","TH_29","TH_19","LN_34","BM_36","TH_34","LN_23","BM_13","TH_13","TH_32")
#levels(set_5)<-lev
#sub_seurat <- subset(set_4, idents = "1")
#cells <- Cells(sub_seurat)
#features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4","MS4A1")
#inFileName<-readRDS("data/Spleen_Combined.rds")
#data_list = list(set_1=set_1,set_2=set_2,set_3=set_3)
data_list = list(set_1=set_1,set_2=set_2,set_3=set_3,set_4=set_4,set_5=set_5)

cd_genes<- c("EBF1","CD19","MS4A1","CD79B","HLA-DOB","CD5","GZMA","CCL5","KLRK1","KLRB1","CD244","NLRP3","CD14","STEAP4","CD163","DEFB1")

genes<- c("ACTB","GAPDH")





ui <- fluidPage(theme = shinytheme("cerulean"),
                
                titlePanel("Visualization of sc-RNAseq Data"),
                navbarPage("FAANG",
                           tabPanel(icon("home"),
                                    
                                    fluidRow(column(tags$img(src="Pig.png",width="250px",height="200px"),width=2),
                                             column(
                                               
                                               br(),
                                               p("
This project aims at understanding pig immune system for food production and translation research. This will provide an immune cell atlas as a basis for future research.
Moreover, it will imporve cell type and tissues specific gene expression data for genetic selection.
", 
                                                 style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                                               br(),
                                               
                                               p("
                                                 Immune tissues were collected from two 6 month old healthy pigs.
                                                 Created clusters of single cell data and proved they are unique and distinguishable. Identified gene expression patterns and markers for different immune cell types.
                                                 Identified tissue specific vs. peripheral immune cell types by comparing against porcine PBMCs
                                                 Identified tissue-specific differences between porcine and human cell types.
                                                 Used canonical markers, porcine PBMC data, and human tissue-specific data to annotate the porcine immune cell atlas.
                                                 ",style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                                               
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
                                    
                                    hr(),
                                    tags$style(".fa-database {color:#e87722}"),
                                    h3(p(em("Datasets"),icon("database",lib = "font-awesome"),style="color:black;text-align:center")),
                                    tabsetPanel(
                                      tabPanel("ICLN98_HumanMeta",
                                               fluidRow(column(DT::dataTableOutput("ICLN98_HumanMeta"),
                                                               width = 12))),
                                      tabPanel("ICLN00_HumanMeta",
                                               fluidRow(column(DT::dataTableOutput("ICLN00_HumanMeta"),
                                                               width = 12))),
                                      tabPanel("THY00_HumanMeta",
                                               fluidRow(column(DT::dataTableOutput("THY00_HumanMeta"),
                                                               width = 12))),
                                      tabPanel("THY98_HumanMeta",
                                               fluidRow(column(DT::dataTableOutput("THY98_HumanMeta"),
                                                               width = 12)))),
                                    
                                    
                                    
                                    hr(),
                                    p(em("Developed by"),br("Tuggle Lab"),style="text-align:center; font-family: times")
                           ),
                           tabPanel("GENES",icon = icon("user"), 
                                    tags$style(".glyphicon-signal {color:#E87722}"),
                                    h3(p(em("VISUALIZATION OF GENES"),icon("signal", lib="glyphicon"),style="color:black;text-align:center")),
                                    br(),
                                    
                                    column(width = 4,
                                           textInput("gene2", label = "Gene Symbol/Ensemble ID", value = "GAPDH")),
                                           submitButton("Update Now", icon=("refresh")),
                                    
                                    column(width = 12,
                                           plotOutput("genePlot2",height=700,width="1200px")),
                                    
                                  
                                    
                                    
                                    column(width = 12,
                                           selectInput("dataset5", label = h3("ALL TISSUES COMBINED"),
                                                       choices = list("Intergrated"="set_5"),
                                                       selected = "set_5"),
                                           
                                           plotOutput("genePlot6",height=700,width="1200px",
                                                      dblclick="plot_1_dblclick",
                                                      brush=brushOpts(
                                                        id="plot1_brush",
                                                        resetOnNew = TRUE
                                                      ))),
                                    column(width = 12,
                                           selectInput("dataset5", label = h3("INTEGRATED"),
                                                       choices = list("Intergrated"="set_3"),
                                                       selected = "set_5"),
                                          
                                           plotOutput("vPlot5",height=700,width="2900px"
                                                      
                                           ))
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                                    
                           )
                           
                           
                )
)






server<- function(input,output){
  output$ICLN98_HumanMeta <- DT::renderDataTable(
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
    colnames=c(" ","orig.ident","nCount_RNA","nFeature_RNA","SampleID","UmiSums","GenesDetected","prcntTop","prcntMito","DuplicatedBarcodes","PassViability","PassGenesDet","PassLibSize","PassBarcodeFreq","PassAll","Scrublet","PassScrub","nCount_SCT","nFeature_SCT","mapping.score","predicted.cell_ontology_class.score","predicted.cell_ontology_class","predicted.free_annotation.score","predicted.free_annotation")
    ))
  output$ICLN00_HumanMeta <- DT::renderDataTable(
    DT::datatable({
      data
    },
    options=list(lengthMenu=list(c(5,15,20),c('5','15','20')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c(" ","orig.ident","nCount_RNA","nFeature_RNA","SampleID","UmiSums","GenesDetected","prcntTop","prcntMito","DuplicatedBarcodes","PassViability","PassGenesDet","PassLibSize","PassBarcodeFreq","PassAll","Scrublet","PassScrub","nCount_SCT","nFeature_SCT","mapping.score","predicted.cell_ontology_class.score","predicted.cell_ontology_class","predicted.free_annotation.score","predicted.free_annotation")
    ))
  output$THY00_HumanMeta <- DT::renderDataTable(
    DT::datatable({
      dat
    },
    options=list(lengthMenu=list(c(5,15,20),c('5','15','20')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c(" ","orig.ident","nCount_RNA","nFeature_RNA","SampleID","UmiSums","GenesDetected","prcntTop","prcntMito","DuplicatedBarcodes","PassViability","PassGenesDet","PassLibSize","PassBarcodeFreq","PassAll","Scrublet","PassScrub","nCount_SCT","nFeature_SCT","mapping.score","predicted.cell_ontology_class.score","predicted.cell_ontology_class","predicted.free_annotation.score","predicted.free_annotation")
    ))
  output$THY98_HumanMeta <- DT::renderDataTable(
    DT::datatable({
      da
    },
    options=list(lengthMenu=list(c(5,15,20),c('5','15','20')),pagelength=10,
                 
                 columnDefs=list(list(className='dt-center',targets="_all"))
    ),
    filter="top",
    selection='multiple',
    style='bootstrap',
    class= 'cell-border stripe',
    rownames=FALSE,
    colnames=c(" ","orig.ident","nCount_RNA","nFeature_RNA","SampleID","UmiSums","GenesDetected","prcntTop","prcntMito","DuplicatedBarcodes","PassViability","PassGenesDet","PassLibSize","PassBarcodeFreq","PassAll","Scrublet","PassScrub","nCount_SCT","nFeature_SCT","mapping.score","predicted.cell_ontology_class.score","predicted.cell_ontology_class","predicted.free_annotation.score","predicted.free_annotation")
    ))
  
  
  datasetInput <- reactive({
    infile <- data_list[[input$dataset]]
    
  })
  dataset2Input <- reactive({
    infile <- data_list[[input$dataset2]]
    
  })
  dataset3Input <- reactive({
    infile <- data_list[[input$dataset3]]
    
  })
  dataset4Input <- reactive({
    infile <- data_list[[input$dataset4]]
    
  })
  dataset5Input <- reactive({
    infile <- data_list[[input$dataset5]]
    
  })
  
  
  #ranges<- reactiveValues(x=NULL, y=NULL)
  

  #output$cex<-renderText(input$size)
  
  
  
  
  output$genePlot2 <- renderPlot({
    FeaturePlot_scCustom(dataset5Input(),split.by="Tissue",num_columns = 2,reduction="Tissue_UMAP",features=(input$gene2))
  })
  
 
 

  
  
  output$genePlot6 <- renderPlot({
    
    FeaturePlot_scCustom(dataset5Input(),reduction="umap",features= (input$gene2))
  })
  output$vPlot5 <- renderPlot({
    VlnPlot(dataset5Input(),split.by="MainCluster",features = (input$gene2),split.plot=TRUE)
    
    
    
  })

}

shinyApp(ui, server)