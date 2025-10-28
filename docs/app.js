/* global Papa, Plotly */
(function () {
  const datasetSelect = document.getElementById('datasetSelect');
  const layoutSelect = document.getElementById('layoutSelect');
  const clusterSelect = document.getElementById('clusterSelect');
  const geneInput = document.getElementById('geneInput');
  const geneShowBtn = document.getElementById('geneShowBtn');
  const geneInfoBtn = document.getElementById('geneInfoBtn');
  const info = document.getElementById('info');

  const state = {
    datasets: [],
    currentDataset: null,
    currentLayout: 'umap',
    coordinates: null, // { x:[], y:[], barcode:[] }
    clusters: null, // Map barcode -> clusterId
    clusterList: [],
    expression: null, // { gene: string, valuesByBarcode: Map, min: number, max: number }
    annotations: { loaded: false, index: null, records: null },
    markers: { loadedFor: null, byCluster: new Map() }, // Map clusterId -> [{ geneId, shortName, score }]
  };

  init();

  async function init() {
    info.textContent = 'Loading datasets…';
    const datasets = await buildDatasets();
    state.datasets = datasets;
    populateDatasetSelect(datasets);
    if (datasets.length > 0) {
      datasetSelect.value = datasets[0].id;
      await loadDataset(datasets[0].id);
    }
    wireEvents();
  }

  function wireEvents() {
    datasetSelect.addEventListener('change', async () => {
      await loadDataset(datasetSelect.value);
    });
    layoutSelect.addEventListener('change', async () => {
      state.currentLayout = layoutSelect.value;
      await renderEmbedding();
    });
    clusterSelect.addEventListener('change', () => {
      renderEmbedding();
    });
    geneShowBtn.addEventListener('click', async () => {
      const gene = geneInput.value.trim();
      if (!gene) {
        // Clear expression overlay
        state.expression = null;
        await renderEmbedding();
        info.textContent = 'Cleared expression overlay.';
        return;
      }
      try {
        await loadExpression(state.currentDataset.id, gene);
        await renderEmbedding();
        // Try to show annotations too
        try {
          await ensureAnnotationsLoaded();
          const rec = findAnnotation(gene);
          if (rec) {
            renderAnnotationRow(rec);
          } else {
            info.textContent = `Expression: ${state.expression.gene} loaded. No annotation found.`;
          }
        } catch (e2) {
          info.textContent = `Expression: ${state.expression.gene} loaded.`;
        }
      } catch (err) {
        console.error(err);
        info.textContent = `Expression not found for gene "${gene}" in ${state.currentDataset.label}. Place CSV at docs/data/${state.currentDataset.id}/expr/${gene}.csv`;
      }
    });

    geneInfoBtn?.addEventListener('click', async () => {
      const gene = geneInput.value.trim();
      if (!gene) {
        info.textContent = 'Enter a gene symbol or LOC ID (e.g., LOC105326593 or Tyr).';
        return;
      }
      try {
        await ensureAnnotationsLoaded();
        const rec = findAnnotation(gene);
        renderAnnotationRow(rec);
      } catch (err) {
        console.error(err);
        info.textContent = 'Failed to load annotations index.';
      }
    });
  }

  function populateDatasetSelect(datasets) {
    datasetSelect.innerHTML = '';
    for (const ds of datasets) {
      const opt = document.createElement('option');
      opt.value = ds.id;
      opt.textContent = ds.label;
      datasetSelect.appendChild(opt);
    }
  }

  async function loadDataset(datasetId) {
    const ds = state.datasets.find(d => d.id === datasetId);
    state.currentDataset = ds;
    state.currentLayout = layoutSelect.value;
    info.textContent = `Loading ${ds.label}…`;

    try {
      const [coords, clusters] = await Promise.all([
        loadCoordinates(ds.paths[state.currentLayout]),
        loadClusters(ds.paths.clusters)
      ]);
      state.coordinates = coords;
      state.clusters = clusters.map;
      state.clusterList = clusters.unique;

      populateClusterSelect(state.clusterList);
      await renderEmbedding();
      info.textContent = `${coords.barcode.length.toLocaleString()} cells loaded.`;
      // Load markers for this dataset (non-blocking for plot)
      try {
        await loadMarkersForDataset(ds.id);
        // Refresh legend to include markers
        renderLegend(buildPalette());
      } catch (e) {
        // Silently ignore if marker files are unavailable
      }
    } catch (err) {
      console.error(err);
      Plotly.purge('scatter');
      document.getElementById('legend').innerHTML = '';
      info.textContent = `Failed to load data: ${err && err.message ? err.message : err}`;
    }
  }

  function populateClusterSelect(uniqueClusters) {
    clusterSelect.innerHTML = '';
    const allOpt = document.createElement('option');
    allOpt.value = 'all';
    allOpt.textContent = 'All clusters';
    clusterSelect.appendChild(allOpt);

    for (const cl of uniqueClusters) {
      const opt = document.createElement('option');
      opt.value = String(cl);
      opt.textContent = `Cluster ${cl}`;
      clusterSelect.appendChild(opt);
    }
    clusterSelect.value = 'all';
  }

  async function renderEmbedding() {
    const coords = state.coordinates;
    const clusters = state.clusters;
    if (!coords || !clusters) return;

    const filterCluster = clusterSelect.value;

    const x = [];
    const y = [];
    const text = [];
    const color = [];
    const z = []; // expression values when present

    // Build color palette
    const palette = buildPalette();

    for (let i = 0; i < coords.barcode.length; i++) {
      const bc = coords.barcode[i];
      const cl = clusters.get(bc);
      if (filterCluster !== 'all' && String(cl) !== filterCluster) continue;
      x.push(coords.x[i]);
      y.push(coords.y[i]);
      text.push(`${bc} | cluster ${cl}`);
      if (state.expression) {
        const v = state.expression.valuesByBarcode.get(bc);
        z.push(typeof v === 'number' ? v : NaN);
      } else {
        color.push(palette.get(String(cl)) || '#cccccc');
      }
    }

    const trace = state.expression ? {
      x,
      y,
      text,
      type: 'scattergl',
      mode: 'markers',
      marker: {
        size: 4,
        color: z,
        colorscale: 'Viridis',
        opacity: 0.9,
        cmin: state.expression.min,
        cmax: state.expression.max,
        colorbar: { title: state.expression.gene }
      },
      hovertemplate: '%{text}<extra></extra>'
    } : {
      x,
      y,
      text,
      type: 'scattergl',
      mode: 'markers',
      marker: { size: 4, color, opacity: 0.85 },
      hovertemplate: '%{text}<extra></extra>'
    };

    const layout = {
      title: state.expression ? `${state.currentDataset.label} – ${state.currentLayout.toUpperCase()} · ${state.expression.gene}` : `${state.currentDataset.label} – ${state.currentLayout.toUpperCase()}`,
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: '#e6ecff' },
      xaxis: { zeroline: false, showgrid: false },
      yaxis: { zeroline: false, showgrid: false },
      margin: { l: 40, r: 10, t: 40, b: 40 }
    };

    Plotly.react('scatter', [trace], layout, { displayModeBar: true, responsive: true });

    renderLegend(palette);
  }

  function renderLegend(palette) {
    const legendEl = document.getElementById('legend');
    legendEl.innerHTML = '';
    if (state.expression) {
      const title = document.createElement('div');
      title.textContent = `Expression: ${state.expression.gene}`;
      legendEl.appendChild(title);
      const hint = document.createElement('div');
      hint.textContent = 'Colorbar shown on plot.';
      hint.style.color = '#a5b4d6';
      hint.style.marginTop = '6px';
      legendEl.appendChild(hint);
      return;
    }
    const title = document.createElement('div');
    title.textContent = 'Clusters';
    title.style.marginBottom = '8px';
    legendEl.appendChild(title);

    const list = document.createElement('div');
    for (const cl of state.clusterList) {
      const row = document.createElement('div');
      row.style.display = 'flex';
      row.style.alignItems = 'center';
      row.style.gap = '8px';
      row.style.marginBottom = '6px';

      const sw = document.createElement('span');
      sw.style.display = 'inline-block';
      sw.style.width = '14px';
      sw.style.height = '14px';
      sw.style.borderRadius = '3px';
      sw.style.background = palette.get(String(cl)) || '#cccccc';

      const label = document.createElement('span');
      label.textContent = `Cluster ${cl}`;

      row.appendChild(sw);
      row.appendChild(label);
      list.appendChild(row);
    }
    legendEl.appendChild(list);

    // Markers section
    const hr = document.createElement('div');
    hr.style.borderTop = '1px solid #233055';
    hr.style.margin = '10px 0';
    legendEl.appendChild(hr);

    const mtitle = document.createElement('div');
    mtitle.textContent = 'Top markers';
    mtitle.style.marginBottom = '8px';
    legendEl.appendChild(mtitle);

    const selectedCluster = clusterSelect.value;
    if (!state.markers.byCluster || state.markers.byCluster.size === 0) {
      const note = document.createElement('div');
      note.textContent = 'No marker table available for this dataset.';
      note.style.color = '#a5b4d6';
      legendEl.appendChild(note);
    } else if (selectedCluster === 'all') {
      const note = document.createElement('div');
      note.textContent = 'Select a cluster to view its top markers.';
      note.style.color = '#a5b4d6';
      legendEl.appendChild(note);
    } else {
      const genes = state.markers.byCluster.get(Number(selectedCluster)) || [];
      const ul = document.createElement('div');
      for (const g of genes.slice(0, 20)) { // limit for brevity
        const row = document.createElement('div');
        row.style.display = 'flex';
        row.style.alignItems = 'center';
        row.style.justifyContent = 'space-between';
        row.style.gap = '8px';
        row.style.marginBottom = '6px';

        const name = document.createElement('button');
        name.textContent = g.shortName || g.geneId;
        name.style.background = 'transparent';
        name.style.border = '1px solid #233055';
        name.style.color = '#e6ecff';
        name.style.borderRadius = '6px';
        name.style.padding = '4px 8px';
        name.style.cursor = 'pointer';
        name.title = `Overlay expression and show annotation for ${g.geneId}`;
        name.addEventListener('click', async () => {
          try {
            await loadExpression(state.currentDataset.id, g.geneId);
            await renderEmbedding();
            try {
              await ensureAnnotationsLoaded();
              const rec = findAnnotation(g.geneId);
              renderAnnotationRow(rec);
            } catch {}
          } catch (err) {
            info.textContent = `Expression not found for gene "${g.geneId}" in ${state.currentDataset.label}.`;
          }
        });

        const score = document.createElement('span');
        score.textContent = g.score != null ? g.score.toFixed(3) : '';
        score.style.color = '#a5b4d6';
        score.style.fontSize = '12px';

        row.appendChild(name);
        row.appendChild(score);
        ul.appendChild(row);
      }
      legendEl.appendChild(ul);
    }
  }

  function buildPalette() {
    const colors = [
      '#4e79a7','#f28e2b','#e15759','#76b7b2','#59a14f','#edc949','#af7aa1','#ff9da7','#9c755f','#bab0ab',
      '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
    ];
    const map = new Map();
    let idx = 0;
    for (const cl of state.clusterList) {
      const key = String(cl);
      if (!map.has(key)) {
        map.set(key, colors[idx % colors.length]);
        idx++;
      }
    }
    return map;
  }

  function getDataBase() {
    const isPages = /github\.io$/i.test(window.location.hostname);
    const repoRawBase = 'https://raw.githubusercontent.com/RobertsLab/pgc-edc-oyster/main';
    return isPages ? `${repoRawBase}/docs/data` : 'data';
  }

  async function ensureAnnotationsLoaded() {
    if (state.annotations.loaded) return;
    const url = `${getDataBase()}/annotations/index.json`;
    const res = await fetch(url);
    if (!res.ok) throw new Error(`GET ${url} → ${res.status}`);
    const payload = await res.json();
    state.annotations.index = payload.index || {};
    state.annotations.records = payload.records || [];
    state.annotations.loaded = true;
  }

  function findAnnotation(geneRaw) {
    if (!geneRaw) return null;
    const key = String(geneRaw).toUpperCase();
    const idx = state.annotations.index && state.annotations.index[key];
    if (typeof idx !== 'number') return null;
    return state.annotations.records[idx] || null;
  }

  function renderAnnotationRow(r) {
    if (!r) {
      info.textContent = 'No annotation found for this gene.';
      return;
    }
    const parts = [];
    parts.push(`Gene: ${r.loc_id || ''}${r.gene_name ? ` — ${r.gene_name}` : ''}`.trim());
    if (r.gene_type) parts.push(`Type: ${r.gene_type}`);
    if (r.protein_accession) parts.push(`NCBI Protein: ${r.protein_accession}`);
    if (r.uniprot) parts.push(`UniProt: ${r.uniprot}${r.uniprot_e ? ` (e=${r.uniprot_e})` : ''}`);
    if (r.dmel) parts.push(`Dmel: ${r.dmel}${r.dmel_e ? ` (e=${r.dmel_e})` : ''}`);
    if (r.cel) parts.push(`Cel: ${r.cel}${r.cel_e ? ` (e=${r.cel_e})` : ''}`);
    if (r.spur) parts.push(`Spur: ${r.spur}${r.spur_e ? ` (e=${r.spur_e})` : ''}`);

    const links = [];
    if (r.uniprot) links.push(`<a href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(r.uniprot)}/entry" target="_blank" rel="noopener">UniProt</a>`);
    if (r.protein_accession) links.push(`<a href="https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(r.protein_accession)}" target="_blank" rel="noopener">NCBI Protein</a>`);

    info.innerHTML = [
      parts.join(' · '),
      links.length ? `Links: ${links.join(' | ')}` : ''
    ].filter(Boolean).join('<br>');
  }

  function getMarkersBase() {
    const isPages = /github\.io$/i.test(window.location.hostname);
    const repoRawBase = 'https://raw.githubusercontent.com/RobertsLab/pgc-edc-oyster/main';
    return isPages
      ? `${repoRawBase}/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis`
      : '../oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis';
  }

  function getMarkerFileForDataset(datasetId) {
    if (/^oyster_E[1-4]_redo2_roslin-mito$/.test(datasetId)) return 'gast_top_specific_markers_top25byMarkerScore.txt';
    if (datasetId === 'oyster_r1and2_Bla_roslin-mito-CRv3') return 'Bla_top_specific_markers_top25byMarkerScore.txt';
    if (/^oyster_r1and2_CP[1-3]_roslin-mito-CRv3$/.test(datasetId)) return 'CPbla_top_specific_markers_top25byMarkerScore.txt';
    return null;
  }

  async function loadMarkersForDataset(datasetId) {
    if (state.markers.loadedFor === datasetId) return;
    const fname = getMarkerFileForDataset(datasetId);
    state.markers.byCluster = new Map();
    state.markers.loadedFor = null;
    if (!fname) return;
    const url = `${getMarkersBase()}/${fname}`;
    const res = await fetch(url);
    if (!res.ok) return; // silently ignore
    const txt = await res.text();
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true, delimiter: '\t' });
    const fields = parsed.meta && parsed.meta.fields ? parsed.meta.fields : [];
    const clusterField = fields.find(f => /cell[ _]?group/i.test(f)) || 'cell_group';
    const geneField = fields.find(f => /^gene_id$/i.test(f)) || fields.find(f => /^gene_ID$/i.test(f)) || 'gene_id';
    const shortField = fields.find(f => /short.*name/i.test(f)) || fields.find(f => /gene_short_name/i.test(f));
    const scoreField = fields.find(f => /marker_score/i.test(f));
    for (const row of parsed.data) {
      if (!row || row[clusterField] == null) continue;
      const cl = Number(row[clusterField]);
      const geneId = String(row[geneField] || '').trim();
      if (!geneId) continue;
      const shortName = shortField ? String(row[shortField] || '').trim() : '';
      const score = scoreField != null ? Number(row[scoreField]) : null;
      if (!state.markers.byCluster.has(cl)) state.markers.byCluster.set(cl, []);
      state.markers.byCluster.get(cl).push({ geneId, shortName, score });
    }
    // Sort each cluster by score desc when available
    for (const [cl, arr] of state.markers.byCluster.entries()) {
      arr.sort((a, b) => (Number(b.score || 0) - Number(a.score || 0)));
    }
    state.markers.loadedFor = datasetId;
  }

  async function loadCoordinates(csvPath) {
    const res = await fetch(csvPath);
    if (!res.ok) {
      throw new Error(`GET ${csvPath} → ${res.status}`);
    }
    const txt = await res.text();
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const x = [];
    const y = [];
    const barcode = [];
    // Expect headers: Barcode,UMAP-1,UMAP-2 or Barcode,tSNE-1,tSNE-2
    const xKey = parsed.meta.fields.find(f => /umap-1|tsne-1/i.test(f));
    const yKey = parsed.meta.fields.find(f => /umap-2|tsne-2/i.test(f));
    for (const row of parsed.data) {
      if (!row || row.Barcode == null) continue;
      barcode.push(String(row.Barcode));
      x.push(Number(row[xKey]));
      y.push(Number(row[yKey]));
    }
    return { x, y, barcode };
  }

  async function loadClusters(csvPath) {
    const res = await fetch(csvPath);
    if (!res.ok) {
      throw new Error(`GET ${csvPath} → ${res.status}`);
    }
    const txt = await res.text();
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const map = new Map();
    const unique = new Set();
    for (const row of parsed.data) {
      if (!row || row.Barcode == null) continue;
      map.set(String(row.Barcode), row.Cluster);
      unique.add(row.Cluster);
    }
    return { map, unique: Array.from(unique).sort((a, b) => Number(a) - Number(b)) };
  }

  async function buildDatasets() {
    // Build a base path that works both locally and on GitHub Pages
    const isPages = /github\.io$/i.test(window.location.hostname);
    const repoRawBase = 'https://raw.githubusercontent.com/RobertsLab/pgc-edc-oyster/main';
    const base = isPages
      ? `${repoRawBase}/oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam`
      : '../oyster_scRNASeq_jobs_genomic_resources_outs/CellRanger_outputs_nobam';
    const dataBase = isPages
      ? `${repoRawBase}/docs/data`
      : 'data';
    const ids = [
      'oyster_r1and2_CP3_roslin-mito-CRv3',
      'oyster_r1and2_CP2_roslin-mito-CRv3',
      'oyster_r1and2_CP1_roslin-mito-CRv3',
      'oyster_r1and2_Bla_roslin-mito-CRv3',
      'oyster_E4_redo2_roslin-mito',
      'oyster_E3_redo2_roslin-mito',
      'oyster_E2_redo2_roslin-mito',
      'oyster_E1_redo2_roslin-mito'
    ];
    const friendly = new Map([
      ['oyster_r1and2_CP3_roslin-mito-CRv3','CP3 roslin-mito (CRv3)'],
      ['oyster_r1and2_CP2_roslin-mito-CRv3','CP2 roslin-mito (CRv3)'],
      ['oyster_r1and2_CP1_roslin-mito-CRv3','CP1 roslin-mito (CRv3)'],
      ['oyster_r1and2_Bla_roslin-mito-CRv3','Bla roslin-mito (CRv3)'],
      ['oyster_E4_redo2_roslin-mito','E4 redo2 roslin-mito'],
      ['oyster_E3_redo2_roslin-mito','E3 redo2 roslin-mito'],
      ['oyster_E2_redo2_roslin-mito','E2 redo2 roslin-mito'],
      ['oyster_E1_redo2_roslin-mito','E1 redo2 roslin-mito']
    ]);
    return ids.map(id => ({
      id,
      label: friendly.get(id) || id,
      paths: {
        umap: `${base}/${id}/outs/analysis/umap/2_components/projection.csv`,
        tsne: `${base}/${id}/outs/analysis/tsne/2_components/projection.csv`,
        clusters: `${base}/${id}/outs/analysis/clustering/graphclust/clusters.csv`,
        exprDir: `${dataBase}/${id}/expr`
      }
    }));
  }

  async function loadExpression(datasetId, geneRaw) {
    const ds = state.datasets.find(d => d.id === datasetId);
    if (!ds) throw new Error('Dataset not found');
    const candidates = [geneRaw, geneRaw.toUpperCase(), geneRaw.toLowerCase()].filter(Boolean);
    let csvText = null;
    let geneUsed = null;
    for (const g of candidates) {
      const url = `${ds.paths.exprDir}/${encodeURIComponent(g)}.csv`;
      const res = await fetch(url);
      if (res.ok) {
        csvText = await res.text();
        geneUsed = g;
        break;
      }
    }
    if (!csvText) throw new Error('Expression CSV not found');

    const parsed = Papa.parse(csvText, { header: true, dynamicTyping: true });
    // Determine expression column: prefer a column literally named geneUsed; otherwise, use second column
    let exprCol = parsed.meta.fields.find(f => f.toLowerCase() === String(geneUsed).toLowerCase());
    if (!exprCol && parsed.meta.fields.length >= 2) exprCol = parsed.meta.fields[1];
    const valuesByBarcode = new Map();
    let min = Infinity, max = -Infinity;
    for (const row of parsed.data) {
      if (!row || row.Barcode == null) continue;
      const v = Number(row[exprCol]);
      if (Number.isFinite(v)) {
        valuesByBarcode.set(String(row.Barcode), v);
        if (v < min) min = v;
        if (v > max) max = v;
      }
    }
    if (!Number.isFinite(min) || !Number.isFinite(max)) { min = 0; max = 1; }
    state.expression = { gene: geneUsed, valuesByBarcode, min, max };
  }
})();


