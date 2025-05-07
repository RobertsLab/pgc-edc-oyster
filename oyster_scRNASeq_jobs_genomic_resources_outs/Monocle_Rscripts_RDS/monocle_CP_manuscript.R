###################RDS generation and top marker analysis for oyster CP_bla single cell data

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
LOCannot<-fread("~/Documents/gigas_annotations/RoslinGenome/LOCannot_homologies.txt", header = T) #script for joining LOC ID to homologies in folder with this file

#read in cell ranger data
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

CP1 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/CellRanger_outs_roslin-mito-CRv3/oyster_r1and2_CP1_roslin-mito-CRv3/") # cleavage pool 1
CP2 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/CellRanger_outs_roslin-mito-CRv3/oyster_r1and2_CP2_roslin-mito-CRv3/") # cleavage pool 2
CP3 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/CellRanger_outs_roslin-mito-CRv3/oyster_r1and2_CP3_roslin-mito-CRv3/") # cleavage pool 3

colData(CP1)$ID = "CP1"
colData(CP2)$ID = "CP2"
colData(CP3)$ID = "CP3"


CP_cds = combine_cds(list(CP1,CP2,CP3), cell_names_unique = F) # use when there are multiple samples
colData(CP_cds)$n.umi = Matrix::colSums(assay(CP_cds))

#filtercells
CP_cds <- CP_cds[,colData(CP_cds)$n.umi >= 2000]

# don't need to re-run if loading processed data 
CP_cds = estimate_size_factors(CP_cds)
CP_cds = detect_genes(CP_cds)

CP_cds = preprocess_cds(CP_cds, num_dim = 20 ) 
CP_cds = align_cds(CP_cds, alignment_group = "sample") # use when there are multiple samples
CP_cds = reduce_dimension(CP_cds, umap.min_dist = 0.2, umap.n_neighbors = 15L)
CP_cds = cluster_cells(CP_cds)

#####################Adding metadata###########################################

#1. add annotations as gene_short_name?
rowDat_df<-rowData(CP_cds) %>% as.data.frame()
rowDat_df<-left_join(rowDat_df, LOCannot, by = c("gene_short_name" = "LOC_ID"))
head(rowDat_df)
rowData(CP_cds)$gene_short_name<-rowDat_df$gene_ID
head(rowData(CP_cds))

#2. add samples,clusters and partitions as metadata
colData(CP_cds)$sample = factor(colData(CP_cds)$sample) 
colData(CP_cds)$cluster = clusters(CP_cds)
colData(CP_cds)$partitions = partitions(CP_cds)

#3. normalized # genes per cell by umi
colData(CP_cds)$num_gen_express_byUMI = colData(CP_cds)$num_genes_expressed/colData(CP_cds)$n.umi
colData(CP_cds)$norm_log10geneexpress = log10(colData(CP_cds)$num_genes_expressed)/colData(CP_cds)$n.umi

#4. log10n.umi
colData(CP_cds)$log10n.umi = log10(colData(CP_cds)$n.umi)

#5. mito fraction
## Annotate fraction of reads that are mitochondrial ##
mito<-c("CYTB", "COX2", "ATP6", "ND2", "ND4", "ND5", "ND6", "ND3", "ND1", "ND4L", "COX1", "COX3")
cds1 = CP_cds[rowData(CP_cds)$gene_short_name %in% mito,]
head(rowData(cds1))
#mt.genes = subset(rowData(CPbla_cds), grepl("^mt-", gene_short_name_mito_annot))$id
mt.genes = rowData(cds1)$id
head(mt.genes)
norm.mat = counts(CP_cds)
mito.mat = norm.mat[mt.genes,]
colData(CP_cds)$mito.sum = Matrix::colSums(mito.mat)
colData(CP_cds)$mito.frac = with(colData(CP_cds), mito.sum / n.umi)

##########taking a look at RDS########
#how many cells?
nrow(colData(CP_cds)) #7,658
#how many genes?
nrow(rowData(CP_cds)) #38,222

#how many cells per cluster and paritition?
cells_per_clus<-as.data.frame(colData(CP_cds)) %>% group_by(cluster) %>% dplyr::summarize(count = n())
cells_per_clus

cells_per_part<-as.data.frame(colData(CP_cds)) %>% group_by(partitions) %>% dplyr::summarize(count = n())
cells_per_part

####PLOT CELLS####
plot_cells(CP_cds, 
           color_cells_by = "partition", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)


#####################Saving RDS###########################################
saveRDS(CP_cds, "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_CP_only.RDS")
