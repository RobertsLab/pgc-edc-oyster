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

#read in gast_cds
gast_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_gastrula.RDS")
head(colData(gast_cds))
gast_cds_df<-as.data.frame(colData(gast_cds))
head(gast_cds_df)
nrow(gast_cds_df)
median(gast_cds_df$num_genes_expressed)
mean(gast_cds_df$n.umi)

#read in annotation file
LOCannot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies_genetable.txt", header = T) #script for joining LOC ID to homologies in folder with this file
head(LOCannot)

plot_cells(gast_cds, 
                   color_cells_by = "num_gen_express_byUMI", 
                   group_label_size = 4, 
                   show_trajectory_graph = F,
                   label_groups_by_cluster = F, 
                   label_cell_groups = T,
                   cell_size = .1,
                   cell_stroke=0.2)

ggsave("numgeneexpressedUMAP_Fig5c.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#subcluster the muscle and then try to look for pseudotime-I feel comfortable anchoring this with num genes expressed
####taking the verticle selection and doing pseudotime analysis####
tmp_cds<-choose_cells(gast_cds) #take the vert selection and root the bottom right
tmp_cds = preprocess_cds(tmp_cds)
tmp_cds = reduce_dimension(tmp_cds)
tmp_cds = cluster_cells(tmp_cds)
tmp_cds = learn_graph(tmp_cds)
tmp_cds = order_cells(tmp_cds)
head(colData(tmp_cds))

#plot_cells(tmp_cds,
#           color_cells_by = "cluster",
#           group_label_size = 6,
#           label_cell_groups=T,
#           show_trajectory_graph=F,
#           scale_to_range = F,
#           min_expr = 0.2,
#           cell_stroke=0.8)

#plot_cells(tmp_cds, 
#           color_cells_by = "num_gen_express_byUMI", 
#           group_label_size = 4, 
#           show_trajectory_graph = F,
#           label_groups_by_cluster = F, 
#           label_cell_groups = T,
#           cell_size = .4,
#           graph_label_size=1.5)

#add back tmp cluster IDs to the full UMAP
#make a new column with the subcluster info
colData(tmp_cds)$assigned_cell_type <- as.character(clusters(tmp_cds)[colnames(tmp_cds)])
#make a new empty column in the original cds
colData(gast_cds)$assigned_cell_type<-'.'
#add the tmp cds data to the full cds
colData(gast_cds)[colnames(tmp_cds),]$assigned_cell_type <- colData(tmp_cds)$assigned_cell_type
#plot it
plot_cells(gast_cds, group_cells_by="cluster", 
           color_cells_by="assigned_cell_type", 
           labels_per_group=5,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           cell_size = .1,
           cell_stroke=0.2) +
  scale_color_hue(direction=-1)

ggsave("muscle_subclusterUMAP_Fig5c.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#add a new meta column to the tmp cds with the subcluster ID
#cluster 1 is base of tail
#cluster 2 is the left side of cat ear
#cluster 3 is the base of right ear
#cluster 4 is the top of right ear


###back to LOC with row data
rowDat_df<-rowData(tmp_cds) %>% as.data.frame()
rowData(tmp_cds)$gene_short_name<-rowDat_df$id

###top marker for all
tmp_marker_test_res <- monocle3::top_markers(tmp_cds, group_cells_by = "cluster", 
                                   reference_cells = 1000, cores=8)
#by markerscore 
tmp_top_specific_markers <- tmp_marker_test_res %>%
  filter(fraction_expressing >= 0.20) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)

tmp_top_specific_marker_ids <- unique(tmp_top_specific_markers %>% pull(gene_id))
tmp_topmarker_df<-as.data.frame(tmp_top_specific_marker_ids)
base <- left_join(tmp_top_specific_markers, LOCannot, by = c("gene_id" = "LOC_ID")) 
tmp_topmarkers_annot_cand<-base %>% pull(gene_short_name)
fwrite(base,"~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/gast_MuscleSubCluster_top5byMarkerScore.txt",sep = "\t",row.names = T)



plot_genes_by_group(tmp_cds,
                             tmp_topmarkers_annot_cand,
                             group_cells_by="cluster",
                             ordering_type="maximal_on_diag",
                             max.size=5)

ggsave("subclusUMAP_Fig5c.png",
       dpi = 750,
       height = 5,
       width = 8,
       bg = "transparent")
