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
LOCannot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies.txt", header = T) #script for joining LOC ID to homologies in folder with this file

#read in cell ranger data
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

Bla_cds = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_r1and2_Bla_roslin-mito-CRv3/") # bla

colData(Bla_cds)$ID = "Bla"

colData(Bla_cds)$n.umi = Matrix::colSums(assay(Bla_cds))

#filtercells
Bla_cds <- Bla_cds[,colData(Bla_cds)$n.umi >= 2000]

#how many cells
nrow(colData(Bla_cds)) #2836

# don't need to re-run if loading processed data 
Bla_cds = estimate_size_factors(Bla_cds)
Bla_cds = detect_genes(Bla_cds)

Bla_cds = preprocess_cds(Bla_cds, num_dim = 20 ) 
Bla_cds = reduce_dimension(Bla_cds, umap.min_dist = 0.2, umap.n_neighbors = 15L)
Bla_cds = cluster_cells(Bla_cds)

#####################Adding metadata###########################################

#1. add annotations as gene_short_name?
rowDat_df<-rowData(Bla_cds) %>% as.data.frame()
rowDat_df<-left_join(rowDat_df, LOCannot, by = c("gene_short_name" = "LOC_ID"))
head(rowDat_df)
rowData(Bla_cds)$gene_short_name<-rowDat_df$gene_ID
head(rowData(Bla_cds))

#2. add clusters and partitions as metadata
colData(Bla_cds)$cluster = clusters(Bla_cds)
colData(Bla_cds)$partitions = partitions(Bla_cds)

#3. normalized # genes per cell by umi
colData(Bla_cds)$num_gen_express_byUMI = colData(Bla_cds)$num_genes_expressed/colData(Bla_cds)$n.umi
colData(Bla_cds)$norm_log10geneexpress = log10(colData(Bla_cds)$num_genes_expressed)/colData(Bla_cds)$n.umi

#4. log10n.umi
colData(Bla_cds)$log10n.umi = log10(colData(Bla_cds)$n.umi)

#5. mito fraction
## Annotate fraction of reads that are mitochondrial ##
mito<-c("CYTB", "COX2", "ATP6", "ND2", "ND4", "ND5", "ND6", "ND3", "ND1", "ND4L", "COX1", "COX3")
cds1 = Bla_cds[rowData(Bla_cds)$gene_short_name %in% mito,]
head(rowData(cds1))
#mt.genes = subset(rowData(Bla_cds), grepl("^mt-", gene_short_name_mito_annot))$id
mt.genes = rowData(cds1)$id
head(mt.genes)
norm.mat = counts(Bla_cds)
mito.mat = norm.mat[mt.genes,]
colData(Bla_cds)$mito.sum = Matrix::colSums(mito.mat)
colData(Bla_cds)$mito.frac = with(colData(Bla_cds), mito.sum / n.umi)

saveRDS(bla_cds, "~/Documents/scRNA_Seq_gigas/jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_blastula_only.RDS")


plot_cells(Bla_cds, 
           color_cells_by = "cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)



##########taking a look at RDS########
#how many cells?
nrow(colData(Bla_cds)) #10,494
#how many genes?
nrow(rowData(Bla_cds)) #38,222

#how many cells per cluster and paritition?
cells_per_clus<-as.data.frame(colData(Bla_cds)) %>% group_by(cluster) %>% dplyr::summarize(count = n())
cells_per_clus

cells_per_part<-as.data.frame(colData(Bla_cds)) %>% group_by(partitions) %>% dplyr::summarize(count = n())
cells_per_part

####PLOT CELLS####
plot_cells(Bla_cds, 
           color_cells_by = "cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)


#####################Saving RDS###########################################
saveRDS(Bla_cds, "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_blastula_only.RDS")
