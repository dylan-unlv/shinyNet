library("BiocManager")
options(repos = BiocManager::repositories())
getOption("repos")
library(biomaRt)
library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(readr)
library(stringr)
library(tidyr)
library(shinyjs)
library(sortable)
library(tibble)
library(DESeq2)
library(visNetwork)
library(igraph)
library(poweRlaw)
library(rclipboard)
library(dplyr)



options(shiny.maxRequestSize = 10000*1024^2)


# Define UI
ui <- shinyUI(navbarPage(title = "shinyNet",
                         theme=shinytheme('flatly'),
                         
                         tabPanel(title = "Upload",
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput("nodefile", label = "Node Matrix (.tsv of nodes and attributes)"),
                                      fileInput("edgefile", label = "Edge Matrix (.tsv of edges and attributes [e.g. weight])"),
                                      actionButton(inputId="run","Analyze Network"),
                                      hr(),
                                      # plotOutput("hist"), 
                                      width = 3), 
                                    
                                    mainPanel(
                                      
                                      DT::dataTableOutput("nodetable"),
                                      hr(),
                                      DT::dataTableOutput("edgetable"), 
                                      width = 9
                                      
                                    )
                                  ) ),
                      tabPanel(title = "Visualize",
                        sidebarLayout(
                        sidebarPanel(
                             uiOutput("modules"),
                             hr(),
                             radioButtons('whichdeg', 
                                          'Select local or global degree to be mapped to node size',
                                          c('Global','Local'), selected='Global'),
                             rclipboardSetup(),
                             hr(), 
                             uiOutput("copymods"),
                             hr(),
                             actionButton(inputId="graph","Visualize Network"),
                             width = 3), 
                                             
                                             
                                mainPanel(
                                               
                                  fluidPage(
                                    visNetworkOutput('graphnet', width='1200px', height='1200px')
                                    )
                                               
                                 )
                           )
                                           
                     ),
                     tabPanel(title = "Hub Genes",
                              sidebarLayout(
                                sidebarPanel(
                                  actionButton(inputId="hubs","Find Hub Genes"),
                                  hr(),
                                  plotOutput("pfit"),
                                  hr(), 
                                  radioButtons('hublevel', 
                                               'Global or Local Stats',
                                               c('Global','Local'), selected='Global'),
                                  
                                  hr(), 
                                  uiOutput("copyhubs"),
                                  width = 3), 
                                
                                
                                mainPanel(
                                  
                                  fluidPage(
                                    plotOutput("plhist"), 
                                    DT::dataTableOutput("hubnodes")
                                  )
                                  
                                )
                              )
                              
                     ),
                     tabPanel(title = "KEGG Analysis",
                                sidebarLayout(
                                  sidebarPanel(
                                    actionButton(inputId="runkegg","Run KEGG"),
                                    hr(), 
                                    radioButtons('kegglevel', 
                                                 'Global or from selected list',
                                                 c('Global','Local'), selected='Global'),
                                    textInput('kegg_glist',label='Paste gene list (one gene each line)'),
                                    width = 3), 
                                  
                                  
                                  mainPanel(
                                    
                                    fluidPage( 
                                      DT::dataTableOutput("hubtable")
                                    )
                                    
                                  )
                                )
                                
                     )
                                  
                                    
         )
  )
                         
                         
                        
