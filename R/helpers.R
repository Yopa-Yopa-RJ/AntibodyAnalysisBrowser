# helpers.R
library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(shiny)
library(dplyr)
library(plotly)
library(DT)
library(ape);
library(NGLVieweR)

#Load Cluster data
sc5p_v2_hs_chain_status <- read_csv("data/sc5p_v2_hs_chain_status.csv")
sc5p_v2_hs_clone_size <- read_csv("data/sc5p_v2_hs_clone_size.csv")
sc5p_v2_hs_sampleid <- read_csv("data/sc5p_v2_hs_sampleid.csv")

#Make scatter Plots
chain_status_plot    <- function(){ renderPlot({
  # Chain status plot
  sc5p_v2_hs_chain_status <- sc5p_v2_hs_chain_status |> mutate( y = ifelse(str_detect(x, ","), str_extract(x, "[\\-\\.0-9]*$"), y), x = ifelse(str_detect(x, ","), str_extract(x, "^[\\-\\.0-9]*"), x), x = as.double(x), y = as.double(y) )
  
  sc5p_v2_hs_chain_status |> ggplot(aes(x, y, color = label)) +
    geom_point() +
    labs(title = "chain_status",
         x = "vdj1",
         y = "vdj2") +
    scale_color_manual(values = c("No_contig" = "orange", "extra contig" = "blue", "single pair" = "brown"))
})
}
clone_id_size_3_plot <- function(){ renderPlot({
  # clone size plot
  sc5p_v2_hs_clone_size <- sc5p_v2_hs_clone_size |> mutate( y = ifelse(str_detect(x, ","), str_extract(x, "[\\-\\.0-9]*$"), y), x = ifelse(str_detect(x, ","), str_extract(x, "^[\\-\\.0-9]*"), x), x = as.double(x), y = as.double(y),label = ifelse(is.na(label), " >=3", as.character(label)) , label = as.factor(label))
  
  sc5p_v2_hs_clone_size |> ggplot(aes(x, y, col = label) ) +
    geom_point() +
    labs(title = "clone_id_size_3",
         x = "vdj1",
         y = "vdj2")+
    scale_color_manual(values = c("1" = "grey", "2" = "blue", " >=3" = "orange"))
})
}
sampleid_plot        <- function(){ renderPlot({
  # sampleid plot
  sc5p_v2_hs_sampleid <- sc5p_v2_hs_sampleid |> mutate( y = ifelse(str_detect(x, ","), str_extract(x, "[\\-\\.0-9]*$"), y), x = ifelse(str_detect(x, ","), str_extract(x, "^[\\-\\.0-9]*"), x), x = as.double(x), y = as.double(y), label = ifelse(str_detect(label, "vdj_nextgem_hs_pbmc"), "vdj_nextgem_hs_pbmc", label) )
  
  sc5p_v2_hs_sampleid |> ggplot(aes(x, y, col = label)) +
    geom_point() +
    labs(title = "sampleid",
         x = "vdj1",
         y = "vdj2") +
    scale_color_manual(values = c("sc5p_v2_hs_PBMC_10k" = "orange", "sc5p_v2_hs_PBMC_1k" = "blue", "vdj_nextgem_hs_pbmc" = "green", "vdj_v1_hs_pbmc" = "brown"))
})
}

# Placeholder metric functions with integration notes

calculate_affinity <- function(id) {10 + runif(1,-5,5)}

calculate_patentability <- function(id) sample(c("High","Medium","Low"),1)

calculate_gravy <- function(seq) NA

calculate_charge <- function(seq) NA

calculate_aggregation <- function(seq) runif(1,0,1)

calculate_tm <- function(seq) 70 + runif(1,-5,5)

calculate_human_similarity <- function(seq) runif(1,80,100)

calculate_tcell_epitopes <- function(seq) sample(0:5,1)

calculate_nh_motifs <- function(seq) sample(0:3,1)

calculate_expression_titer <- function(seq) round(runif(1,1,10),1)

calculate_glyco_sites <- function(seq) sample(0:3,1)

calculate_sec_purity <- function(seq) round(runif(1,90,100),1)

calculate_half_life <- function(seq) round(runif(1,100,300),1)

calculate_fcrn <- function(seq) runif(1,0.1,1)

calculate_ic50 <- function(seq) round(100/runif(1,1,10),1)

calculate_cross_reactivity <- function(seq) runif(1,0,1)

calculate_specificity <- function(seq) runif(1,0,1)

calculate_toxicity <- function(seq) sample(c("Low","Medium","High"),1)

calculate_shm_distance <- function(seq) runif(1,0,0.2)

calculate_branch_length <- function(seq) runif(1,0,0.5)

calculate_germline_divergence <- function(seq) runif(1,0,0.2)

