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

#read in CPbla_cds
CPbla_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_cleavage_blastula.RDS")

#read in bla_cds
bla_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_blastula_only.RDS")

#read in CP_cds
CP_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_CP_only.RDS")

#read in annotation file
LOCannot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies_genetable.txt", header = T) #script for joining LOC ID to homologies in folder with this file


#read in gastrula cluster meta data
clusterAnnot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/gastrula_cluster_layer_annot.txt")

#join the cluster metadata to the gastula cds
colDat_df<-colData(gast_cds) %>% as.data.frame()
clusterAnnot$cluster<-as.factor(clusterAnnot$cluster)
colDat_df<-left_join(colDat_df, clusterAnnot, by = c("cluster" = "cluster"))
head(colDat_df)
colData(gast_cds)$germ_layer<-colDat_df$germ_layer
colData(gast_cds)$cell_type<-colDat_df$cell_type
colData(gast_cds)$ClusterID<-colDat_df$ClusterID
colData(gast_cds)$GEM_short<-colDat_df$GEM_short
colData(gast_cds)$CelTyp_short<-colDat_df$CelTyp_short
colData(gast_cds)$CelTyp_dotAnnot<-colDat_df$CelTyp_dotAnnot
colData(gast_cds)$shell_subcluster<-colDat_df$shell_subcluster
colData(gast_cds)$shell_annot<-colDat_df$shell_annot
colData(gast_cds)$cell_type_broadNeural<-colDat_df$cell_type_broadNeural
colData(gast_cds)$neural_annot<-colDat_df$neural_annot
colData(gast_cds)$dev_dotplot<-colDat_df$dev_dotplot
colData(gast_cds)$CelTyp_General<-colDat_df$CelTyp_General
head(colData(gast_cds))


#####MAIN FIGURES#########################################################
########FIGURE 1. UMAPS Cleav_bla and gastrula
plot_cells(gast_cds, 
           color_cells_by = "cluster", 
           group_label_size = 5, 
           show_trajectory_graph = F,
           label_groups_by_cluster = T, 
           label_cell_groups = F,
           cell_size = .1,
           cell_stroke=.5) 
          

ggsave("UMAPgastrulae.png",
       dpi = 750,
       height = 8,
       width = 8,
       bg = "transparent")

plot_cells(CPbla_cds,
           color_cells_by = "cluster", 
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("UMAPcleav_bla.png",
       dpi = 750,
       height = 8,
       width = 8,
       bg = "transparent")

########FIGURE 2. PGC 
#vasa plot B
vasa<- c("LOC105335166") # vasa

plot_cells(gast_cds,
           genes=vasa,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("vasa_UMAPgastrulae.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#spPHI plot B
spPHI<-c("LOC105327445") #spPHI
plot_cells(gast_cds,
           genes=spPHI,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("spPHI_UMAPgastrulae.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#unchr8839 plot B
unchr8839<-c("LOC105328839") #unchr8839
plot_cells(gast_cds,
           genes=unchr8839,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("unchr8839_UMAPgastrulae.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")


########cleavage_blastulae B
plot_cells(CPbla_cds,
           genes=vasa,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("vasa_UMAPcleav_bla.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#spPHI plot
plot_cells(CPbla_cds,
           genes=spPHI,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("spPHI_UMAPcleav_bla.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#unchr8839 plot
plot_cells(CPbla_cds,
           genes=unchr8839,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("unchr8839_UMAPcleav_bla.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#semiq is in a different script - barplot_semiquant.R

# dot plot

#pull the top 10 marker genes for the PGC cluster (cell group 19 = PGC)
PGC_markers<-c("LOC105328839","LOC105347430","LOC105327445","LOC105318512","LOC105335618","LOC105321539","LOC105345191","LOC105321610","LOC105346676","LOC109619157","LOC105335166")

#plot it 
plot_genes_by_group(gast_cds,
                    PGC_markers,
                    group_cells_by="CelTyp_short",
                    ordering_type = "maximal_on_diag",
                    max.size = 6, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right") +
  xlab("")

ggsave("PGCdotplot_Fig2.png",
       dpi = 750,
       height = 5,
       width = 10,
       bg = "transparent")


#Fig3A. large clustered marker gene diagram - THIS IS SAVED AS SEPARATE SCRIPT scRNA_dotplot_clustered.R

#c. UMAPs of genes that characterize global patterns
twist<-c("LOC105325904")
BMP3<-c("LOC105319691")
brachyury<-c("LOC105339662")
root_CROCC<-c("LOC105331920")
snail2b<-c("LOC105320154") 
Tkt12_3<-c("LOC105328534")
caudal<-c("LOC105323015")  
goosecoid<-c("LOC105327239") 

#update gene and file name 1 at a time Fig 3b
plot_cells(gast_cds,
           genes=Tkt12_3,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.3)

ggsave("Tkt12_3UMAP_Fig3.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")



#Figure 4. Biomineralization (paneled umap plots)
#big picture
plot_cells(gast_cds, 
           color_cells_by="shell_annot",
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 1,
           cell_size = .2,
           cell_stroke=0.3) +
           scale_color_hue(direction=-1)

ggsave("shell_clustersUMAP_Fig4.png",
       dpi = 750,
       height = 5,
       width = 7,
       bg = "transparent")

#genes for plotting in umaps - plot individual and combine later to scale each gene separately
chitinsynthase2<-c("LOC105324613")
prisilkin<-c("LOC117689669")
mantle<-c("LOC105326721")

#subset UMAPS
shell_tmp_cds<-choose_cells(gast_cds)

#plot and save individual genes separately
plot_cells(shell_tmp_cds,
           genes=prisilkin,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.6)

ggsave("prisilkin_Fig4.png",
       dpi = 750,
       height = 4,
       width = 6,
       bg = "transparent")

#shell dotplot
#as dot plot?
#back to LOC_ID for gene short name
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
rowData(gast_cds)$gene_short_name<-rowDat_df$id

shellgenes<-c("LOC105342064","LOC117689669","LOC105318856","LOC105326721","LOC105324157","LOC105324613","LOC105317751")

cell_order <- c("SF1","SF4","SF3", "SF2")

# subsetting to shell clusters for dot plot of select genes
plot_shell_cds <- gast_cds[,(colData(gast_cds)$shell_subcluster == "shell field")] 
colData(plot_shell_cds)$CelTyp_short <- factor(colData(plot_shell_cds)$CelTyp_short, levels = cell_order)
colData(gast_cds)$CelTyp_short <- factor(colData(gast_cds)$CelTyp_short, levels = cell_order)

# plot
monocle3::plot_genes_by_group(gast_cds,
                              shellgenes,
                              ordering_type = 'none',
                              group_cells_by = "CelTyp_short",
                              max.size = 6, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  xlab("")

ggsave("shell_dotplot_Fig4.png",
       dpi = 750,
       height = 4,
       width = 7,
       bg = "transparent")

#Figure 5. Muscle clusters (target gene umaps and subcluster analysis)
#umaps
#forpanelA
MELHst<-c("LOC105317061")
MELHsmooth<-c("LOC105320997")
#forpanelB
otp<-c("LOC105334613")
snail<-c("LOC105328595")
MHC<-c("LOC105338907")

#plot each one at a time
plot_cells(gast_cds,
           genes=otp,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.3)

ggsave("otp_Fig5b.png",
       dpi = 750,
       height = 6,
       width = 7,
       bg = "transparent")

#traj analysis
plot_cells(gast_cds, 
           color_cells_by = "num_gen_express_byUMI", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)


#1. see monocle_subcluster_analysis.R for now

#Figure 6. Neural clusters###################################################
#big picture
#plot_cells(gast_cds, 
#           color_cells_by="neural_annot",
#           group_label_size = 4,
#           label_cell_groups=F,
#           show_trajectory_graph=F,
#           scale_to_range = F,
#           min_expr = 1,
#           cell_size = .2,
#           cell_stroke=0.4)


#main figure save one at a time
dachshund<-c("LOC105340880")
sodium_serotonin<-c("LOC105340406")
six3_6<-c("LOC105327496")
orthopedia<-c("LOC105334613")
prospero<-c("LOC105345931")
HTa<-c("LOC105348009")

#plot each - updated gene name and save file
monocle3::plot_cells(gast_cds,
           genes=dachshund,
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.5,
           cell_size = .2,
           cell_stroke=.3)

ggsave("dachs_Fig6.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")


#as dot plot?
#back to LOC_ID for gene short name
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
rowData(gast_cds)$gene_short_name<-rowDat_df$id

neural_dot<-c("LOC105340880","LOC105327496","LOC105340406","LOC105348058","LOC105334171","LOC105348737","LOC105332393","LOC105344310","LOC105348692","LOC105345931","LOC105334613","LOC105348009","LOC105345175")

cell_order <- c("Mus1","Mus2","Circ1","Circ2","Circ3",
                "Neu8","Neu7","Neu1","Neu5","Neu6","Neu3","Neu2","Neu4",
                "Cil4","Cil2","Cil1","Cil3",
                "SF3","SF2","SF1","SF4",
                "PGC")
#or
# subsetting to neural clusters for dot plot of select genes
plot_neural_cds <- gast_cds[,(colData(gast_cds)$cell_type_broadNeural == "neural")] 
colData(plot_neural_cds)$CelTyp_short <- factor(colData(plot_neural_cds)$CelTyp_short, levels = cell_order)
colData(gast_cds)$CelTyp_short <- factor(colData(gast_cds)$CelTyp_short, levels = cell_order)


# plot
monocle3::plot_genes_by_group(gast_cds,
                              neural_dot,
                              ordering_type = 'none',
                              group_cells_by = "CelTyp_short",
                              max.size = 8, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  xlab("")


ggsave("dotplot_neural_Fig6.png",
       dpi = 750,
       height = 5,
       width = 10,
       bg = "transparent")


#Figure 7 cleavage blastula

#top 3 marker genes for clusters 9 and 10 - Figure 7c
#back to LOC_ID for gene short name
rowDat_df<-rowData(CPbla_cds) %>% as.data.frame()
rowData(CPbla_cds)$gene_short_name<-rowDat_df$id

clus9and10<-c("LOC105333366","LOC109617465","LOC105335420","LOC105318387","LOC105346085","LOC105339662")

plot_genes_by_group(CPbla_cds,
                    clus9and10,
                    group_cells_by="cluster",
                    ordering_type = "none",
                    max.size = 6, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right") +
  xlab("")

ggsave("dotplot_cleavbla_TFsdotplot_Fig7.png",
       dpi = 750,
       height = 5,
       width = 5,
       bg = "transparent")

#Fig7b - umap for cilliated cells
efHand<- c("LOC105337196")

monocle3::plot_cells(CPbla_cds,
                     genes=efHand,
                     group_label_size = 4,
                     label_cell_groups=F,
                     show_trajectory_graph=F,
                     scale_to_range = F,
                     min_expr = 0.5,
                     cell_size = .2,
                     cell_stroke=.4)

ggsave("spef1UMAP_Fig7.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")



#Fig 7A norm gene expressed 
plot_cells(CPbla_cds, 
          color_cells_by = "norm_log10geneexpress", 
          group_label_size = 12, 
          show_trajectory_graph = F,
          label_groups_by_cluster = T, 
          label_cell_groups = F)

ggsave("cleav_blaUMAP_numgeneexpressed_Fig7.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")


#####SUPP FIGURES#########################################################

########SUPP FIGURE 1. PDFs from cell ranger

######Supp 2 dotplot of PGC in cleavage blastula

#dot plot
#back to LOC_ID for gene short name
rowDat_df<-rowData(CPbla_cds) %>% as.data.frame()
rowData(CPbla_cds)$gene_short_name<-rowDat_df$id
#pull the top 10 marker genes for the PGC cluster (cell group 19 = PGC)
PGC_markers<-c("LOC105327445", "LOC105347430", "LOC105328839", "LOC105321610", "LOC105346676", "LOC105321539", "LOC105345191", "LOC105318512", "LOC109619157", "LOC105335618", "LOC105335166")

#plot it 
plot_genes_by_group(CPbla_cds,
                    PGC_markers,
                    group_cells_by="cluster",
                    ordering_type = "none",
                    max.size = 6, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right") +
  xlab("")

ggsave("dotplot_PGC_FigS2.png",
       dpi = 750,
       height = 8,
       width = 6,
       bg = "transparent")

#SUPP FIGURE 3. colored by cleav and bla ----add new column that is either CP or Bla and plot
#make a new column with either CP or Bla to color cells by bla sample
coldata_df = colData(CPbla_cds) %>% as.data.frame()
head(coldata_df)
coldata_df = coldata_df %>% 
  mutate(sample_group = case_when(ID %in% c("CP1", "CP2", "CP3") ~ "CP",
                                  ID == "early_Bla" ~ "Bla",
                                  TRUE ~ "Other"))

colData(CPbla_cds)$sample_group = coldata_df$sample_group

plot_cells(CPbla_cds, 
              color_cells_by = "sample_group", 
              group_label_size = 4, 
              show_trajectory_graph = F,
              label_groups_by_cluster = F, 
              label_cell_groups = F)

ggsave("cleav_or_bla_FigS3.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")


#SUPP FIGURE 4. negative expression for neural markers
#supp Figure 4 save 1 at a time
FMRF<-c("LOC105348058")
hex<-c("LOC105339244")
      
      
#plot each - updated gene name and save file
monocle3::plot_cells(gast_cds,
                     genes=hex,
                     group_label_size = 4,
                     label_cell_groups=F,
                     show_trajectory_graph=F,
                     scale_to_range = F,
                      min_expr = 0.5,
                     cell_size = .2,
                     cell_stroke=.4)
      
      ggsave("hex_FigS4.png",
             dpi = 750,
             height = 5,
             width = 6,
             bg = "transparent")
      
#SUPP FIGURE 5. blastula only and cleavage-blastula colored by blastula cluster
      
plot_cells(bla_cds,
            color_cells_by = "cluster", 
            group_label_size = 4,
            label_cell_groups=F,
            show_trajectory_graph=F,
            scale_to_range = F,
            min_expr = 0.2,
            cell_size = .2,
            cell_stroke=0.4)
      
ggsave("blastula_UMAP_FigS5.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")


#color cells by cluster of the blastula only subset
bla_cds$barcode <-paste0(bla_cds$barcode,"_4")
colData(CPbla_cds)$blastula_cluster <- NA
colData(CPbla_cds)[bla_cds$barcode,]$blastula_cluster <- colData(bla_cds)$cluster
colData(CPbla_cds)$blastula_cluster <- as.factor(colData(CPbla_cds)$blastula_cluster)


plot_cells(CPbla_cds,
           color_cells_by = "blastula_cluster", 
           group_label_size = 4,
           label_cell_groups=F,
           show_trajectory_graph=F,
           scale_to_range = F,
           min_expr = 0.2,
           cell_size = .2,
           cell_stroke=0.4)

ggsave("clev_blastula_UMAP_byblaOnlyCluster_FigS5.png",
       dpi = 750,
       height = 5,
       width = 6,
       bg = "transparent")

#make a new column with either CP or Bla to color cells by bla sample
      coldata_df = colData(CPbla_cds) %>% as.data.frame()
      head(coldata_df)
      coldata_df = coldata_df %>% 
        mutate(sample_group = case_when(ID %in% c("CP1", "CP2", "CP3") ~ "CP",
                                        ID == "early_Bla" ~ "Bla",
                                        TRUE ~ "Other"))
      
      colData(CPbla_cds)$sample_group = coldata_df$sample_group
      
      plot_cells(CPbla_cds, 
                 color_cells_by = "sample_group", 
                 group_label_size = 4, 
                 show_trajectory_graph = F,
                 label_groups_by_cluster = F, 
                 label_cell_groups = F)
      
      ggsave("bla_FigS3.png",
             dpi = 750,
             height = 5,
             width = 6,
             bg = "transparent")      
      

#############################################NOT IN PAPER BUT SAVE
###save this because I like it
plot_cells(gast_cds, 
           color_cells_by = "CelTyp_General", 
           group_label_size = 5, 
           show_trajectory_graph = F,
           label_groups_by_cluster = T, 
           label_cell_groups = F,
           cell_size = .1,
           cell_stroke=0.3)+
  scale_color_manual(values=c("#73B000","#FE6D8D","#998EFF","#E28900","#00A8FF","#EE67EC","#808080"))

#temp 3d umap ~~~to see neural cell types better (picking colors based on cell type annotations<need to adjust some of these)
#3d
cds_3d <- reduce_dimension(gast_cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)

gg_color_hue_list_ordered<-c("#EE67EC","#00B9E3","#998EFF","#BF80FF","#E28900","#73B000","#ACA300","#FE6D8D","#00A8FF","#808080")
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="cell_type",color_palette = gg_color_hue_list_ordered)
cds_3d_plot_obj




