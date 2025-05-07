suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(gridExtra)
  library(viridis)
  library(devtools)
  library(monocle3)
  library(data.table)
  library(gridExtra)
  library(gghighlight)
  library(RColorBrewer)
  library(cowplot)
  library(ggridges)
})
#####READ in RDSs and META DATA#########################################################
#read in gast_cds
gast_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_gastrula.RDS")
head(rowData(gast_cds))
#read in CPbla_cds
CPbla_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_cleavage_blastula.RDS")

#read in bla_cds
Bla_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_blastula_only.RDS")

#read in annotation file
LOCannot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies_genetable.txt", header = T) #script for joining LOC ID to homologies in folder with this file

###top marker for GASTRULA##########################################

###back to LOC with row data - so you don't get duplicates
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
rowData(gast_cds)$gene_short_name<-rowDat_df$id

gast_marker_test_res <- monocle3::top_markers(gast_cds, group_cells_by = "cluster", 
                                              reference_cells = 1000, cores=8)
#pull top25 according to marker_score
gast_top_specific_markers <- gast_marker_test_res %>%
  filter(fraction_expressing >= 0.2) %>%
  group_by(cell_group) %>%
  top_n(25, marker_score)

#add annotations information
head(gast_top_specific_markers)
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
gast_top_specific_markers_annot<-left_join(gast_top_specific_markers, LOCannot, by = c("gene_short_name" = "LOC_ID"))
#quick look
View(gast_top_specific_markers_annot)

#save the file
fwrite(gast_top_specific_markers_annot,"~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis/gast_top_specific_markers_top25byMarkerScore.txt",sep = "\t")


###top marker for CLEAVAGE BLASTULA##########################################
head(rowData(CPbla_cds))
###back to LOC with row data - so you don't get duplicates
rowDat_df<-rowData(CPbla_cds) %>% as.data.frame()
rowData(CPbla_cds)$gene_short_name<-rowDat_df$id

CPbla_marker_test_res <- monocle3::top_markers(CPbla_cds, group_cells_by = "cluster", 
                                               reference_cells = 1000, cores=8)
#pull top25 according to marker_score
CPbla_top_specific_markers <- CPbla_marker_test_res %>%
  filter(fraction_expressing >= 0.2) %>%
  group_by(cell_group) %>%
  top_n(25, marker_score)

#add annotations information
head(CPbla_top_specific_markers)
rowDat_df<-rowData(CPbla_cds) %>% as.data.frame()
CPbla_top_specific_markers_annot<-left_join(CPbla_top_specific_markers, LOCannot, by = c("gene_short_name" = "LOC_ID"))
#quick look
View(CPbla_top_specific_markers_annot)

#save the file
fwrite(CPbla_top_specific_markers_annot,"~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis/CPbla_top_specific_markers_top25byMarkerScore.txt",sep = "\t")

###top marker for BLASTULA ONLY##########################################
head(rowData(Bla_cds))
###back to LOC with row data - so you don't get duplicates
rowDat_df<-rowData(Bla_cds) %>% as.data.frame()
rowData(Bla_cds)$gene_short_name<-rowDat_df$id

Bla_marker_test_res <- monocle3::top_markers(Bla_cds, group_cells_by = "cluster", 
                                             reference_cells = 1000, cores=8)
#pull top25 according to marker_score
Bla_top_specific_markers <- Bla_marker_test_res %>%
  filter(fraction_expressing >= 0.2) %>%
  group_by(cell_group) %>%
  top_n(25, marker_score)

#add annotations information
head(Bla_top_specific_markers)
rowDat_df<-rowData(Bla_cds) %>% as.data.frame()
Bla_top_specific_markers_annot<-left_join(Bla_top_specific_markers, LOCannot, by = c("gene_short_name" = "LOC_ID"))
#quick look
View(Bla_top_specific_markers_annot)

#save the file
fwrite(Bla_top_specific_markers_annot,"~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis/Bla_top_specific_markers_top25byMarkerScore.txt",sep = "\t")
