Expression CSV placement
========================

Place per‑gene expression CSVs here to enable the Gene (feature) search overlay.

Directory layout:

docs/data/<dataset>/expr/<GENE>.csv

CSV format (header row required):

Barcode,<GENE>
AAACCCAAGACGGTTG-1,0.0
AAACCCAAGAGGATCC-1,1.23
...

Notes:
- <dataset> must match one of the dataset IDs used in docs/app.js (e.g., oyster_E1_redo2_roslin-mito).
- <GENE> is case-insensitive; the app tries the exact text, UPPER, and lower case.
- When deployed on GitHub Pages, these files are read from https://raw.githubusercontent.com/RobertsLab/pgc-edc-oyster/main/docs/data/…

