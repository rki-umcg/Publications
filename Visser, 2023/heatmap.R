# General information ----
# Project: 
# Title: Heatmap, RNA-seq data, Annie Visser
# Author: Rick Wilbrink
# Department: Rheumatology and Clinical Immunology
# Affiliation: University Medical Center Groningen
# Email: r.wilbrink01@umcg.nl
# Collaboration: please ask permission from the author before using this script
# Vignette: https://jokergoo.github.io/ComplexHeatmap-reference/book/

# Importing libraries ----
library(gridtext)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(readr)
library(gplots)
library(matrixStats)
library(xlsx)
library(dplyr)
library(methods)

# Set working directory ----
setwd("D:/Onedrive/PhD/R/Projects/Heatmaps/Gene study AV/output")

# Importing data & Data wrangling and transformation ----
df <- read_excel("D:/Onedrive/PhD/R/Projects/Heatmaps/Gene study AV/data/data 20-04-2023/ht_data.xlsx") 
df <- as.data.frame(df) # create a data.frame
rnames <- df[,1] # create rnames that cotain the gene name information
rownames(df) <- rnames # change the original rownames (none so far) to the variable rnames
df <- data.matrix(df) # a matrix is used for the heatmap
df <- df[,-1] # remove the first column that has redundant information
df <- t(scale(t(df))) # create z-scores

# Metadata ----

# columns: import metadata
column_metadata <- read_excel("D:/Onedrive/PhD/R/Projects/Heatmaps/Gene study AV/data/data 20-04-2023/metadata.xlsx")

# create heatmap annotation variables
column_ann <- data.frame(column_metadata$group,
                         column_metadata$SSA,
                         column_metadata$SSB,
                         #column_metadata$`MX1`,
                         column_metadata$`IFN3`,
                         column_metadata$`IFN12`,
                         column_metadata$`IFNM1.2`,
                         column_metadata$`IFNM3.4`,
                         column_metadata$`IFNM5.12`) 

# change column names(if necessary)
colnames(column_ann) <- c('Group', 
                          'SSA',
                          'SSB',
                          #'MX1',
                          'IFN3',
                          'IFN12',
                          'IFNM1.2',
                          'IFNM3.4',
                          'IFNM5.12') 


# HA colors ----
#col.MX1 <- colorRamp2(c(0, 0.5, 1, 1.5, 2), c("#ffffff", "#c7d9f7", "#8eb3ef", "#568de6", "#1d67de"))
col.IFN3 <- colorRamp2(c(-5, 0, 15, 30, 45), c('#658ff2','#ffffff', '#fdcc8a','#fc8d59','#e34a33'))
col.IFN12 <- colorRamp2(c(-20, 0, 25, 50, 85), c('#658ff2','#ffffff', '#fdcc8a','#fc8d59','#e34a33'))
col.IFN.1.2 <- colorRamp2(c(-10, 0, 20, 40, 60), c('#658ff2','#ffffff', '#fdcc8a','#fc8d59','#e34a33'))
col.IFN.3.4 <- colorRamp2(c(-5, 0, 7, 14, 21), c('#658ff2','#ffffff', '#fdcc8a','#fc8d59','#e34a33'))
col.IFN.5.12 <- colorRamp2(c(-5, 0, 5, 10, 15), c('#658ff2','#ffffff', '#fdcc8a','#fc8d59','#e34a33')) 
col.CD45 <- colorRamp2(c(0, 6, 13), c('#e2e2e2','#687793','#0f3961'))

column_colours <- list('Group' = c('C' = "#0000FF", 'pSS' = "#FF0004"), 
                       'SSA' = c('pos' = '#FF8F00', 'neg' = '#00FF2D'),
                       'SSB' = c('pos' = '#FF8F00', 'neg' = '#00FF2D'),
                      # 'MX1' = col.MX1,
                       'IFN3' = col.IFN3,
                       'IFN12' = col.IFN12,
                       'IFNM1.2' = col.IFN.1.2,
                       'IFNM3.4' = col.IFN.3.4,
                       'IFNM5.12'= col.IFN.5.12)


# Heatmap column annotation ----
HA <- HeatmapAnnotation(df= column_ann,
                        which = 'column',
                        col = column_colours,
                        annotation_width = unit(c(1, 4), 'cm'),
                        gap = unit(1, 'mm'),
                        show_annotation_name = TRUE,
                        show_legend =TRUE,
                        gp = gpar(col = "black"),
                        annotation_name_side = 'left',
                        annotation_name_gp = gpar(col = 'black', fontsize = 12, fontface = 'bold'),   
                        annotation_legend_param = list(name = TRUE,          
                                                       title_gp = gpar(col = "black", fontface = "bold", fontsize = 12),
                                                       grid_height = unit(5, "mm"),
                                                       grid_width = unit(10, "mm"),
                                                       labels_gp = gpar(col = "black", fontface = "bold", fontsize = 12),
                                                       legend_height = unit(5, "cm"))
                        )

# Heatmap color ---- 
col = colorRamp2(c(-2, 0, 2), c("green", "black" , "red"))

# Complexheatmap
HT <- Heatmap(df,
              
              #Color heatmap
              col = col, # add color to your heatmap
              
              #Border around the heatmap
              border = TRUE,
              
              ###### Rows ######
              #Row names
              row_names_gp = gpar(fontsize=12, fontface="bold"),
              show_row_names = TRUE,
              row_names_centered = TRUE,
              row_names_side = "left",
              
              #Dendograms (branches)
              show_row_dend=TRUE,
              row_dend_side = "left",
              # row_names_gp = gpar(cex=fontsize(1)),
              row_dend_width = unit(2, "cm"),
              
              #Clustering
              cluster_rows = TRUE,
              clustering_distance_rows ="euclidean", 
              clustering_method_rows = "ward.D",
              
              ###### Columns ######
              
              #Column names
              column_names_gp = gpar(fontsize=12, fontface="bold"),
              column_names_side = "bottom",
              column_names_rot = 45,
              
              ##Dendograms (branches)
              show_column_dend=TRUE,
              column_dend_side = "top",
              column_dend_gp = gpar(cex=fontsize(1)),
              column_dend_height = unit(2, "cm"),
              
              #Clustering
              cluster_columns = TRUE,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "ward.D",
              
              #Order
              #column_order = column_order,
              top_annotation = HA, 
              
              #Legend
              heatmap_legend_param = list(
                title = gt_render("Z-scores"),
                title_gp = gpar(col = "black", fontface = "bold", fontsize = 12),
                #at = c(-2, 2),
                #labels = c("low", "high"),
                labels_gp = gpar(col = "black", fontface = "bold", fontsize = 12),
                legend_height = unit(5, "cm"))
)

HT

# Other visualizations ----
p1 <- draw(HT, annotation_legend_side = "bottom",)

# Save figure ----
tiff('heatmap_AV_01_03_2023.tiff', units="in", width=12, height=8, res=800, compression = 'lzw')
print(HT)
dev.off()

# Session info ----
sessionInfo()
