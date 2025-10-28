# pgc-edc-oyster

## GitHub Pages: Oyster scRNA‑seq Explorer

An interactive, static explorer is available in `docs/` and can be published via GitHub Pages. It visualizes UMAP/t‑SNE embeddings from Cell Ranger analysis and colors cells by clusters.


### Data wiring

The explorer automatically references these paths for each dataset:

- `outs/analysis/umap/2_components/projection.csv`
- `outs/analysis/tsne/2_components/projection.csv`
- `outs/analysis/clustering/graphclust/clusters.csv`

Datasets currently included:

- `oyster_r1and2_CP3_roslin-mito-CRv3`
- `oyster_r1and2_CP2_roslin-mito-CRv3`
- `oyster_r1and2_CP1_roslin-mito-CRv3`
- `oyster_r1and2_Bla_roslin-mito-CRv3`
- `oyster_E4_redo2_roslin-mito`
- `oyster_E3_redo2_roslin-mito`
- `oyster_E2_redo2_roslin-mito`
- `oyster_E1_redo2_roslin-mito`


### Optional: Gene expression visualization

Expression layers are not bundled for page weight reasons. To add expression:

1. Export per‑cell expression for a set of genes into CSV files with headers `Barcode,<GENE>`.
2. Place them under `docs/data/<dataset>/<gene>.csv` and extend `loadCoordinates`/`loadClusters` with a `loadExpression` function.
3. Wire the `Show` button in `docs/index.html` to overlay expression values as a color scale.


--- 

Repository associated with manuscript:

**Primordial Germ Cell Specification and Early Developmental Cell States in Pacific Oyster**

>Mackenzie R. Gavery, Lauren Vandepas, Lauren M. Saunders, Brent Vadopalas, J. Adam Luckenbach, Cole Trapnell, Steven Roberts. Primordial germ cell specification and early developmental cell states in Pacific oyster. BMC Genomics 26, 951 (2025). https://doi.org/10.1186/s12864-025-12122-7
---

## Repository Structure

This repository contains single-cell RNA sequencing (scRNA-seq) analysis data and scripts for studying primordial germ cell specification and early developmental cell states in Pacific oyster (*Crassostrea gigas*).

### Main Directory: `oyster_scRNASeq_jobs_genomic_resources_outs/`

#### 1. `annotation_files_script_forMonocle3/`
Contains *C. gigas* gene annotation files used for Monocle3 analysis:
- Gene annotation tables with LOC identifiers
- Homology tables linking to other species (e.g., *Strongylocentrotus purpuratus*)
- R script for joining annotation files (`annotations_joining_SRhomologtable.R`)
- NCBI gene table downloads

**Key files:**
- `LOCannot_homologies.txt` - Gene homology annotations
- `GCF_902806645.1_genetable_downloaded_011924.tsv` - NCBI gene table
- `Cg_hits-red_SpurSupp.tab` - Cross-species gene comparisons

#### 2. `CellRanger_outputs_nobam/`
Output files from 10x Genomics CellRanger analysis (BAM files excluded to save space):
- **Gastrula stage samples:** `oyster_E1_redo2_roslin-mito` through `oyster_E4_redo2_roslin-mito`
- **Cleavage and blastula samples:**
  - `oyster_r1and2_Bla_roslin-mito-CRv3` (Blastula)
  - `oyster_r1and2_CP1_roslin-mito-CRv3` through `oyster_r1and2_CP3_roslin-mito-CRv3` (Cleavage period)

These directories contain filtered/raw feature-barcode matrices and analysis outputs used as input for Monocle3 analysis.

#### 3. `Monocle_Rscripts_RDS/`
Monocle3 analysis scripts and RDS (R data serialization) files:

**Main analysis scripts:**
- `monocle_gastrula_manuscript.R` - Gastrula stage analysis
- `monocle_CP_manuscript.R` - Cleavage period analysis  
- `monocle_bla_manuscript.R` - Blastula stage analysis
- `monocle_CPbla_manuscript.R` - Combined cleavage and blastula analysis
- `monocle_figure_April2025.R` - Figure generation script
- `monocle_subcluster_analysis.R` - Subcluster analysis
- `MarkerGenes_markerscore.R` - Marker gene scoring analysis
- `scRNA_dotplot_clustered.R` - Dot plot visualization
- `barplot_semiquant.R` - Semi-quantitative bar plots

**RDS data objects:**
- `oyster_CP_only.RDS` - Cleavage period cell dataset
- `oyster_blastula_only.RDS` - Blastula stage cell dataset
- `oyster_cleavage_blastula.RDS` - Combined dataset

**Subdirectories:**
- `RDSgeneration_scripts_andRDSfiles/` - Scripts for generating RDS objects
- `MarkerScoreAnalysis/` - Marker gene analysis results (top markers by score for different stages)

#### 4. `sc_gigas_jobs_genomicresources/`
CellRanger job scripts and genome reference files:

**Job scripts:**
- `count_CellRanger_E1_redo2.job` through `count_CellRanger_E4_redo2.job` - Gastrula samples
- `count_CellRanger_Bla_Roslin-mito.job` - Blastula sample
- `count_CellRanger_CP1_Roslin-mito.job` through `count_CellRanger_CP3_Roslin-mito.job` - Cleavage period samples
- `CellRanger_mkref_Roslin.job` - Reference genome creation

**Subdirectory:**
- `Cgigas_Roslin_mkref_genome-CRv3/` - CellRanger genome reference (STAR index and GTF files)

Note: Large genome/reference files (.fa, .gtf, STAR index files) are excluded via `.gitignore`

#### 5. `SRA_submission_docs/`
Sequence Read Archive (SRA) submission metadata files:
- Sample metadata spreadsheets
- FASTQ file metadata and quality checks
- Multiple versions tracking submission revisions

**Note:** These files document data submission but are not required for reproducing the analysis.

#### 6. `manuscript_docs/`
Supplementary materials for the manuscript:
- `supplementary_figures_singlecell_oyster.docx` - Supplementary figures
- `supplementary_files_singlecell_oyster.xlsx` - Supplementary data tables

---

## Additional Files

- `pgc-edc-oyster.Rproj` - RStudio project file
- `.gitignore` - Excludes large genome files, CellRanger outputs, and R workspace files

---

## Data Analysis Workflow

1. **CellRanger processing:** Raw sequencing data processed using scripts in `sc_gigas_jobs_genomicresources/`
2. **Gene annotation:** Annotation files prepared in `annotation_files_script_forMonocle3/`
3. **Monocle3 analysis:** Cell clustering, trajectory analysis, and marker identification using scripts in `Monocle_Rscripts_RDS/`
4. **Visualization:** Figures generated using various R scripts for manuscript preparation
