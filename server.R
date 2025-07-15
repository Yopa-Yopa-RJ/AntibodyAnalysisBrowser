# server.R
library(shiny)
library(dplyr)
library(plotly)
library(DT)
library(ape);
library(NGLVieweR)
source("R/helpers.R")

data <- read.csv("data/vdj_data.csv",stringsAsFactors=FALSE)
sc5p_v2_hs_chain_status <- read_csv("data/sc5p_v2_hs_chain_status.csv")
sc5p_v2_hs_clone_size <- read_csv("data/sc5p_v2_hs_clone_size.csv")
sc5p_v2_hs_sampleid <- read_csv("data/sc5p_v2_hs_sampleid.csv")

server <- function(input, output, session) {
 filt <- reactive({
   df <- data
   if(!"ALL"%in%input$chain) df <- df[df$chain_status%in%input$chain,]
   if(!"ALL"%in%input$donor) df <- df[df$donor_sample%in%input$donor,]
   if(!"ALL"%in%input$isotype) df <- df[df$isotype_status%in%input$isotype,]
   if(input$prod!="ALL") df <- df[df$productivity==input$prod,]
   df <- df[df$clone_size>=input$clone[1]&df$clone_size<=input$clone[2],]
   if(!"ALL"%in%input$vgene) df <- df[df$variable_gene%in%input$vgene,]
   df <- df[df$cdr_h3_length>=input$cdr[1]&df$cdr_h3_length<=input$cdr[2],]
   df <- df[df$median_shm_rate>=input$shm[1]&df$median_shm_rate<=input$shm[2],]
   df
 })
 output$diversity <- renderText(round(calculate_clonal_diversity(filt()),2))
 output$switch    <- renderText(paste0(round(calculate_class_switching(filt()),1),"%"))
 output$avgclone  <- renderText(round(calculate_avg_clone_size(filt()),1))
 output$affmat    <- renderText(round(calculate_affinity_maturation(filt()),2))
 output$clinical  <- renderText("N/A"); output$market <- renderText("N/A")
 
 #Cluster Plots
 output$chain_status    <- chain_status_plot()
 output$clone_id_size_3 <- clone_id_size_3_plot()
 output$sampleid        <- sampleid_plot()
 
 output$scatter  <- renderPlotly({p<-ggplot(filt(),aes(clone_size,median_shm_rate,
                              color=isotype_status,text=sequence_id))+geom_point()+theme_minimal()
                              ggplotly(p,tooltip="text")})
 output$table    <- renderDT(filt(),options=list(pageLength=5))
 
 observeEvent(input$inspect_btn, { updateNavbarPage(session, "Antibody Analysis", selected="inspect") })
 
 output$inspect_ui<-renderUI({
   req(input$inspect_id)
   seq=input$inspect_id
   tagList(
    h3(paste("Inspect",seq)),
    plotOutput("phylo"),
    NGLVieweROutput("structure_viewer",height="300px"),
    tableOutput("metrics")
   )
 })
 
 output$structure_viewer <- renderNGLVieweR({ NGLVieweR("7CID") %>% addRepresentation("cartoon") }) ###
 
 output$phylo <- renderPlot({plot.phylo(rtree(5))})
 #output$structure <- renderNGLVieweR(NGLVieweR("1IGT", height="400px"))
 output$metrics <- renderTable({
   seq <- input$inspect_id
   data.frame(Metric=c("Affinity","Patent","GRAVY","Charge","Aggregation","Tm",
                       "HumanSim","TcellEp","NHmotifs","ExprTiter","Glyco","SEC","HalfLife","FcRn",
                       "IC50","CrossReact","Specificity","Toxicity","SHMdist","Branch","GermDiv"),
              Value=c(
                calculate_affinity(seq),calculate_patentability(seq),
                calculate_gravy(seq),calculate_charge(seq),calculate_aggregation(seq),calculate_tm(seq),
                calculate_human_similarity(seq),calculate_tcell_epitopes(seq),calculate_nh_motifs(seq),
                calculate_expression_titer(seq),calculate_glyco_sites(seq),calculate_sec_purity(seq),
                calculate_half_life(seq),calculate_fcrn(seq),calculate_ic50(seq),
                calculate_cross_reactivity(seq),calculate_specificity(seq),calculate_toxicity(seq),
                calculate_shm_distance(seq),calculate_branch_length(seq),calculate_germline_divergence(seq)
              ))
 })
}
#shinyApp(ui,server)