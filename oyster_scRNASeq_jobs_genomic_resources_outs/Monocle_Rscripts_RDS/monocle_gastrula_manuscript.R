###################RDS generation and top marker analysis for oyster gastrula single cell data

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
})

####READ IN ANNOTATION FILE####
LOCannot<-fread("~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies.txt", header = T) #script for joining LOC ID to homologies in folder with this file

#read in cell ranger data
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

E1 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_E1_redo2_roslin-mito/") 
E2 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_E2_redo2_roslin-mito/") 
E3 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_E3_redo2_roslin-mito/") 
E4 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_E4_redo2_roslin-mito/") 

colData(E1)$ID = "E1"
colData(E2)$ID = "E2"
colData(E3)$ID = "E3"
colData(E4)$ID = "E4"


gast_cds = combine_cds(list(E1,E2,E3,E4), cell_names_unique = F) # use when there are multiple samples
colData(gast_cds)$n.umi = Matrix::colSums(assay(gast_cds))

#filtercells
gast_cds <- gast_cds[,colData(gast_cds)$n.umi >= 2000]

# don't need to re-run if loading processed data 
gast_cds = estimate_size_factors(gast_cds)
gast_cds = detect_genes(gast_cds)

gast_cds = preprocess_cds(gast_cds, num_dim = 20 ) 
gast_cds = align_cds(gast_cds, alignment_group = "sample") # use when there are multiple samples
gast_cds = reduce_dimension(gast_cds, umap.min_dist = 0.2, umap.n_neighbors = 15L)
gast_cds = cluster_cells(gast_cds)

#####################Adding metadata###########################################

#1. add annotations as gene_short_name?
rowDat_df<-rowData(gast_cds) %>% as.data.frame()
rowDat_df<-left_join(rowDat_df, LOCannot, by = c("gene_short_name" = "LOC_ID"))
head(rowDat_df)
rowData(gast_cds)$gene_short_name<-rowDat_df$gene_ID
head(rowData(gast_cds))

#2. add samples,clusters and partitions as metadata
colData(gast_cds)$sample = factor(colData(gast_cds)$sample) 
colData(gast_cds)$cluster = clusters(gast_cds)
colData(gast_cds)$partitions = partitions(gast_cds)

#3. normalized # genes per cell by umi
colData(gast_cds)$num_gen_express_byUMI = colData(gast_cds)$num_genes_expressed/colData(gast_cds)$n.umi
colData(gast_cds)$norm_log10geneexpress = log10(colData(gast_cds)$num_genes_expressed)/colData(gast_cds)$n.umi

#4. log10n.umi
colData(gast_cds)$log10n.umi = log10(colData(gast_cds)$n.umi)

#5. mito fraction
## Annotate fraction of reads that are mitochondrial ##
mito<-c("CYTB", "COX2", "ATP6", "ND2", "ND4", "ND5", "ND6", "ND3", "ND1", "ND4L", "COX1", "COX3")
cds1 = gast_cds[rowData(gast_cds)$gene_short_name %in% mito,]
head(rowData(cds1))
#mt.genes = subset(rowData(gast_cds), grepl("^mt-", gene_short_name_mito_annot))$id
mt.genes = rowData(cds1)$id
head(mt.genes)
norm.mat = counts(gast_cds)
mito.mat = norm.mat[mt.genes,]
colData(gast_cds)$mito.sum = Matrix::colSums(mito.mat)
colData(gast_cds)$mito.frac = with(colData(gast_cds), mito.sum / n.umi)

#take a quick look
hist(colData(gast_cds)$mito.frac)


##########taking a look at RDS########
#how many cells?
nrow(colData(gast_cds)) #18511
#how many genes?
nrow(rowData(gast_cds)) #38,222

#how many cells per cluster and paritition?
cells_per_clus<-as.data.frame(colData(gast_cds)) %>% group_by(cluster) %>% dplyr::summarize(count = n())
cells_per_clus

cells_per_part<-as.data.frame(colData(gast_cds)) %>% group_by(partitions) %>% dplyr::summarize(count = n())
cells_per_part

####PLOT CELLS####
plot_cells(gast_cds, 
           color_cells_by = "cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)


#####################Saving RDS###########################################
saveRDS(gast_cds, "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_gastrula.RDS")
