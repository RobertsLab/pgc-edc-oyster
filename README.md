### pgc-edc-oyster

_Repository associated with manuscript:_

# Primordial Germ Cell Specification and Early Developmental Cell States in Pacific Oyster 

Authors: Mackenzie R. Gavery, Lauren Vandepas, Lauren M. Saunders, Brent Vadopalas, J. Adam Luckenbach, Cole Trapnell, Steven Roberts




### Directory Contents

1. annotation_files_script_forMonocle3 <-includes C. gigas LOC# annotation files that are read into Monocle3

2. CellRnager_outputs_nobam <-outputs of CellRanger (Monocle3 scripts (see 3.) point here), .bam files not included to save space

3. Monocle_Rscript_RDS <- Monocle3 scripts used for 'gastrula', 'CPbla' and 'figures', folder also includes RDS files that are generated and used in 'figures' script
							
4. sc_gigas_jobs_genomicresources <-CellRanger jobs for CP,bla and gastrula (gastrula labeled as E1-E4), also includes the genome info used in the CellRanger script

5. SRA_submission_docs <- SRA meta data information. Does not need to be included for reproducibility