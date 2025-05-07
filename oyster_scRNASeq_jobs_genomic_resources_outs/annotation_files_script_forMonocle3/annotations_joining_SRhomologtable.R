LOCannot<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/GCF_902806645.1_cgigas_uk_roslin_v1_rna.fna_carrotfiltered.reordered.nodupes.txt", header = T)
head(LOCannot)

LOC_homologies<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/Cg_hits-red_SpurSupp.tab", header = T)
head(LOC_homologies)

trim_homologies<-LOC_homologies[,c(1,3:10)]
head(trim_homologies)

LOCannot_homologies <- left_join(LOCannot, trim_homologies, by = c("LOC_ID" = "LOC")) 
head(LOCannot_homologies)

#01/19/24 join to the new "gene table" available on the genome websit to link LOC to a protein ID
LOC_proteintable<-fread("~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/GCF_902806645.1_genetable_downloaded_011924.tsv", header = T)
head(LOC_proteintable)

trim_genetable<-LOC_proteintable[,c(7,9,11)]
head(trim_genetable)

#join back to the homologies table
LOCannot_homologies_genetable<- left_join(LOCannot_homologies, trim_genetable, by = c("LOC_ID" = "Symbol")) 
head(LOCannot_homologies_genetable)

fwrite(LOCannot_homologies_genetable,"~/Documents/scRNA_Seq_gigas/oyster_scRNASeq_jobs_genomic_resources_outs/annotation_files_script_forMonocle3/LOCannot_homologies_genetable.txt")
