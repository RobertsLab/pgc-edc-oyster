/* global Papa, Plotly */
(function () {
  const datasetSelect = document.getElementById('datasetSelect');
  const layoutSelect = document.getElementById('layoutSelect');
  const clusterSelect = document.getElementById('clusterSelect');
  const geneInput = document.getElementById('geneInput');
  const geneShowBtn = document.getElementById('geneShowBtn');
  const info = document.getElementById('info');

  const state = {
    datasets: [],
    currentDataset: null,
    currentLayout: 'umap',
    coordinates: null, // { x:[], y:[], barcode:[] }
    clusters: null, // Map barcode -> clusterId
    clusterList: [],
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
    geneShowBtn.addEventListener('click', () => {
      // Placeholder: expression not wired (Matrix is large). Could wire via API or precomputed features.
      const gene = geneInput.value.trim();
      if (!gene) return;
      alert('Expression visualization is not bundled in this static demo. Provide per-cell expression CSVs to enable.');
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

    // Build color palette
    const palette = buildPalette();

    for (let i = 0; i < coords.barcode.length; i++) {
      const bc = coords.barcode[i];
      const cl = clusters.get(bc);
      if (filterCluster !== 'all' && String(cl) !== filterCluster) continue;
      x.push(coords.x[i]);
      y.push(coords.y[i]);
      text.push(`${bc} | cluster ${cl}`);
      color.push(palette.get(String(cl)) || '#cccccc');
    }

    const trace = {
      x,
      y,
      text,
      type: 'scattergl',
      mode: 'markers',
      marker: { size: 4, color, opacity: 0.85 },
      hovertemplate: '%{text}<extra></extra>'
    };

    const layout = {
      title: `${state.currentDataset.label} – ${state.currentLayout.toUpperCase()}`,
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
        clusters: `${base}/${id}/outs/analysis/clustering/graphclust/clusters.csv`
      }
    }));
  }
})();


