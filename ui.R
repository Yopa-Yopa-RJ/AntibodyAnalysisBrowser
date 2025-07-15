# ui.R
library(shiny)
library(plotly)
library(DT)
library(NGLVieweR)
source("R/helpers.R")

ui <- navbarPage("Antibody Analysis", id = "Antibody Analysis",
 tabPanel("Home",
  sidebarLayout(
   sidebarPanel(width=3,
    selectInput("chain","Chain Status",choices=c("ALL","single pair","single chain","extra pair"),
                selected="ALL",multiple=FALSE),
    selectInput("donor","Donor",choices=c("ALL","Mouse01","Mouse02","Mouse03"),
                selected="ALL",multiple=FALSE),
    selectInput("isotype","Isotype",choices=c("ALL","IgA","IgG","IgM"),selected="ALL",multiple=FALSE),
    selectInput("prod","Productivity",choices=c("ALL","TRUE","FALSE"),selected="ALL"),
    sliderInput("clone","Clone Size",min=1,max=50,value=c(1,50)),
    selectInput("vgene","Variable Gene",choices=c("ALL","IGHV1","IGHV3","IGKV1"),selected="ALL",multiple=FALSE),
    sliderInput("cdr","CDR-H3 Length",min=0,max=30,value=c(0,30)),
    sliderInput("shm","SHM Rate",min=0,max=20,value=c(0,20)),
    textInput("inspect_id","Inspect ID","seq_1"),
    actionButton("inspect_btn","Inspect")
   ),
   mainPanel(
    fluidRow(
     column(3,wellPanel(h4("Diversity"),textOutput("diversity"))),
     column(3,wellPanel(h4("Switch%"),textOutput("switch"))),
     column(3,wellPanel(h4("AvgClone"),textOutput("avgclone"))),
     #column(2,wellPanel(h4("AffMat"),textOutput("affmat"))),
     column(3,wellPanel(h4("ClinSucc"),textOutput("clinical"))),
     column(3,wellPanel(h4("MarketPot"),textOutput("market")))
    ),
    plotOutput("chain_status",height="500px"),
    plotOutput("clone_id_size_3",height="500px"),
    plotOutput("sampleid",height="500px"),
    plotlyOutput("scatter",height="500px"),
    DTOutput("table")
   )
  )
 ),
 tabPanel("Inspect", value="inspect",
  uiOutput("inspect_ui"),
 ),
 #tabPanel("3D Structure", NGLVieweROutput("structure_viewer", height="400px"))
)