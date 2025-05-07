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
  library(Polychrome)
  library(circlize)
  library(presto)
  library(ComplexHeatmap)
})

#read in gast_cds
#*start with the gast_cds from the R script: monocle_figure_May2024

###back to LOC with row data - so you don't get duplicates
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
rowData(gast_cds)$gene_short_name<-rowDat_df$id

gast_marker_test_res <- monocle3::top_markers(gast_cds, group_cells_by = "cluster", 
                                    reference_cells = 1000, cores=8)
#by marker score metric
gast_top_specific_markers <- gast_marker_test_res %>%
  filter(fraction_expressing >= 0.2) %>%
  group_by(cell_group) %>%
  top_n(3, marker_score)

gast_top_specific_markers_LOC<-gast_top_specific_markers %>% pull(gene_short_name)
head(gast_top_specific_markers_LOC)


############
#clustered dot plot for gastrula - Lauren's code line #444 https://github.com/lsaund11/zfish-skin/blob/main/Aman_Saunders_2023_Fig1.R
###top marker for all
#default plot
top_gast<-plot_genes_by_group(gast_cds,
                              gast_top_specific_markers_LOC,
                              group_cells_by="CelTyp_dotAnnot",
                              ordering_type="cluster_row_col",
                              max.size=5)

# save genes and cell order <- I don't have a gene order preference yet to blocking that out
overall_dotplot<-c("LOC105322588","LOC105319984",#"LOC105345029","LOC105341434","LOC105334836",
                   "LOC105322697","LOC105322336","LOC105324613","LOC105325948","LOC105340661","LOC105334469",
                   "LOC105345408","LOC105339321","LOC105345833","LOC105333508","LOC105327243","LOC105336945",
                   "LOC105338907","LOC105330417","LOC105318379","LOC105342146","LOC105339152","LOC105324245",
                   "LOC105332393","LOC105344310","LOC105344183","LOC105340413","LOC105327018","LOC105348692",
                   "LOC105339411","LOC105348477","LOC105325311","LOC117689456","LOC105318857","LOC105331844",
                   "LOC105334171","LOC105336795","LOC105348737","LOC105333253","LOC105323015","LOC105347171",
                   "LOC105328839","LOC105347430","LOC105327445","LOC117683259","LOC105349206","LOC105338896",
                   "LOC105332654","LOC105340469","LOC105339354","LOC105343336","LOC105333580","LOC105339545",
                   "LOC105334990","LOC105328000","LOC105321897","LOC105337048","LOC105344501","LOC105336057",
                   "LOC105322816","LOC105326722","LOC105344238","LOC105338820","LOC105318663","LOC105333710",
                   "LOC105332606","LOC105328005","LOC105318652")

cell_order <- c("Circ3","Circ2","Circ1",
                "Neu8","Neu7","Neu1","Neu5","Neu6","Neu3","Neu2","Neu4",
                "Cil4","Cil2","Cil1","Cil3",
                "SF3","SF2","SF1","SF4",
                "PGC",
                "Mus1","Mus2")

# clean and apply factors to order cells
plot_cds <- gast_cds[,colData(gast_cds)$CelTyp_dotAnnot != "Un"]
colData(plot_cds)$CelTyp_short <- factor(colData(plot_cds)$CelTyp_short, levels = cell_order)

# plot
monocle3::plot_genes_by_group(plot_cds,
                              overall_dotplot,
                              ordering_type = 'none',
                              group_cells_by = "CelTyp_short",
                              max.size = 4, scale_max = 2, scale_min = -2) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6)) +
  xlab("")

ggsave("gastrula_Fig3dotplot.png",
       dpi = 750,
       height = 9,
       width =8,
       bg = "transparent")

