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

CP1 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_r1and2_CP1_roslin-mito-CRv3/") # cleavage pool 1
CP2 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_r1and2_CP2_roslin-mito-CRv3/") # cleavage pool 2
CP3 = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_r1and2_CP3_roslin-mito-CRv3/") # cleavage pool 3
Bla = load_cellranger_data(pipestance_path = "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam/oyster_r1and2_Bla_roslin-mito-CRv3/") # bla

colData(Bla)$ID = "early_Bla"
colData(CP1)$ID = "CP1"
colData(CP2)$ID = "CP2"
colData(CP3)$ID = "CP3"


CPbla_cds = combine_cds(list(CP1,CP2,CP3,Bla), cell_names_unique = F) # use when there are multiple samples
colData(CPbla_cds)$n.umi = Matrix::colSums(assay(CPbla_cds))

#filtercells
CPbla_cds <- CPbla_cds[,colData(CPbla_cds)$n.umi >= 2000]

# don't need to re-run if loading processed data 
CPbla_cds = estimate_size_factors(CPbla_cds)
CPbla_cds = detect_genes(CPbla_cds)

CPbla_cds = preprocess_cds(CPbla_cds, num_dim = 20 ) 
CPbla_cds = align_cds(CPbla_cds, alignment_group = "sample") # use when there are multiple samples
CPbla_cds = reduce_dimension(CPbla_cds, umap.min_dist = 0.2, umap.n_neighbors = 15L)
CPbla_cds = cluster_cells(CPbla_cds)

#####################Adding metadata###########################################

#1. add annotations as gene_short_name?
rowDat_df<-rowData(CPbla_cds) %>% as.data.frame()
rowDat_df<-left_join(rowDat_df, LOCannot, by = c("gene_short_name" = "LOC_ID"))
head(rowDat_df)
rowData(CPbla_cds)$gene_short_name<-rowDat_df$gene_ID
head(rowData(CPbla_cds))

#2. add samples,clusters and partitions as metadata
colData(CPbla_cds)$sample = factor(colData(CPbla_cds)$sample) 
colData(CPbla_cds)$cluster = clusters(CPbla_cds)
colData(CPbla_cds)$partitions = partitions(CPbla_cds)

#3. normalized # genes per cell by umi
colData(CPbla_cds)$num_gen_express_byUMI = colData(CPbla_cds)$num_genes_expressed/colData(CPbla_cds)$n.umi
colData(CPbla_cds)$norm_log10geneexpress = log10(colData(CPbla_cds)$num_genes_expressed)/colData(CPbla_cds)$n.umi

#4. log10n.umi
colData(CPbla_cds)$log10n.umi = log10(colData(CPbla_cds)$n.umi)

#5. mito fraction
## Annotate fraction of reads that are mitochondrial ##
mito<-c("CYTB", "COX2", "ATP6", "ND2", "ND4", "ND5", "ND6", "ND3", "ND1", "ND4L", "COX1", "COX3")
cds1 = CPbla_cds[rowData(CPbla_cds)$gene_short_name %in% mito,]
head(rowData(cds1))
#mt.genes = subset(rowData(CPbla_cds), grepl("^mt-", gene_short_name_mito_annot))$id
mt.genes = rowData(cds1)$id
head(mt.genes)
norm.mat = counts(CPbla_cds)
mito.mat = norm.mat[mt.genes,]
colData(CPbla_cds)$mito.sum = Matrix::colSums(mito.mat)
colData(CPbla_cds)$mito.frac = with(colData(CPbla_cds), mito.sum / n.umi)

#6. add blastula only cluster info

#load bla only cds and add cluster IDs as meta
bla_cds<-readRDS("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_blastula_only.RDS")


plot_cells(bla_cds, 
           color_cells_by = "cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)

#add relevant joining info
bla_cds$barcode <-paste0(bla_cds$barcode,"_4")
head(colData(bla_cds))
#add column of bla cluster
colData(CPbla_cds)$blastula_cluster <- NA
colData(CPbla_cds)[bla_cds$barcode,]$blastula_cluster <- colData(bla_cds)$cluster
colData(CPbla_cds)$blastula_cluster <- as.factor(colData(CPbla_cds)$blastula_cluster)

head(rowData(CPbla_cds))

plot_cells(CPbla_cds, 
           color_cells_by = "blastula_cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)

##########taking a look at RDS########
#how many cells?
nrow(colData(CPbla_cds)) #10,494
#how many genes?
nrow(rowData(CPbla_cds)) #38,222

#how many cells per cluster and paritition?
cells_per_clus<-as.data.frame(colData(CPbla_cds)) %>% group_by(cluster) %>% dplyr::summarize(count = n())
cells_per_clus

cells_per_part<-as.data.frame(colData(CPbla_cds)) %>% group_by(partitions) %>% dplyr::summarize(count = n())
cells_per_part

####PLOT CELLS####
plot_cells(CPbla_cds, 
           color_cells_by = "cluster", 
           group_label_size = 4, 
           show_trajectory_graph = F,
           label_groups_by_cluster = F, 
           label_cell_groups = T)


#####################Saving RDS###########################################
saveRDS(CPbla_cds, "~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/oyster_cleavage_blastula.RDS")
