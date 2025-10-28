(function () {
  /* global Papa, Plotly */

  const dsSelect = document.getElementById('datasetSelect2');
  const qcPlot = document.getElementById('qc_plot');
  const clusterSizesPlot = document.getElementById('cluster_sizes');
  const pcaScreePlot = document.getElementById('pca_scree');
  const pcaScatterPlot = document.getElementById('pca_scatter');
  const hvfTable = document.getElementById('hvf_table');
  const markerStageSel = document.getElementById('markerStage');
  const markerTopNSel = document.getElementById('markerTopN');
  const markersTable = document.getElementById('markers_table');
  const deClusterSel = document.getElementById('deCluster');
  const deLoadBtn = document.getElementById('deLoadBtn');
  const deTable = document.getElementById('de_table');
  const deVolcano = document.getElementById('de_volcano');
  const smallMultiplesBtn = document.getElementById('renderSmallMultiples');
  const smallMultiplesGrid = document.getElementById('small_multiples');

  const state = {
    datasets: [],
    currentDataset: null,
    clusterList: [],
    // cache
    qcMetrics: new Map(), // id -> metrics object
  };

  init();

  async function init() {
    state.datasets = buildDatasets();
    populateDatasetSelect(state.datasets);
    wireEvents();
    // Render QC overview comparing datasets
    await renderQcOverview(state.datasets);
    // Select first dataset and render its figures
    if (state.datasets.length > 0) {
      dsSelect.value = state.datasets[0].id;
      await onDatasetChange();
    }
  }

  function wireEvents() {
    dsSelect.addEventListener('change', async () => {
      await onDatasetChange();
    });
    markerStageSel.addEventListener('change', () => {
      renderMarkerTable();
    });
    markerTopNSel.addEventListener('change', () => {
      renderMarkerTable();
    });
    deLoadBtn.addEventListener('click', () => {
      renderDifferentialExpression();
    });
    smallMultiplesBtn.addEventListener('click', () => {
      renderSmallMultiples();
    });
  }

  function populateDatasetSelect(datasets) {
    dsSelect.innerHTML = '';
    for (const ds of datasets) {
      const opt = document.createElement('option');
      opt.value = ds.id;
      opt.textContent = ds.label;
      dsSelect.appendChild(opt);
    }
  }

  async function onDatasetChange() {
    const ds = state.datasets.find(d => d.id === dsSelect.value);
    state.currentDataset = ds;
    await renderClusterSizes(ds);
    await renderPcaScree(ds);
    await renderPcaScatter(ds);
    await renderHvfTable(ds);
    await populateDeClusters(ds);
    renderMarkerTable();
  }

  function isPagesHost() {
    return /github\.io$/i.test(window.location.hostname);
  }

  function getRepoRawBase() {
    return 'https://raw.githubusercontent.com/RobertsLab/pgc-edc-oyster/main';
  }

  function buildDatasets() {
    const isPages = isPagesHost();
    const repoRawBase = getRepoRawBase();
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
      base,
      paths: {
        metrics: `${base}/${id}/outs/metrics_summary.csv`,
        clusters: `${base}/${id}/outs/analysis/clustering/graphclust/clusters.csv`,
        umap: `${base}/${id}/outs/analysis/umap/2_components/projection.csv`,
        tsne: `${base}/${id}/outs/analysis/tsne/2_components/projection.csv`,
        pcaVariance: `${base}/${id}/outs/analysis/pca/10_components/variance.csv`,
        pcaProjection: `${base}/${id}/outs/analysis/pca/10_components/projection.csv`,
        hvf: `${base}/${id}/outs/analysis/pca/10_components/features_selected.csv`,
        diffexp: `${base}/${id}/outs/analysis/diffexp/graphclust/differential_expression.csv`,
        webSummary: `${base}/${id}/outs/web_summary_${id.includes('E') ? id.split('_')[1] : 'summary'}.html`
      }
    }));
  }

  // web summary link removed per request

  async function renderQcOverview(datasets) {
    const rows = await Promise.all(datasets.map(async ds => {
      try {
        const metrics = await loadMetrics(ds.paths.metrics);
        state.qcMetrics.set(ds.id, metrics);
        return { id: ds.id, label: ds.label, metrics };
      } catch (e) {
        return { id: ds.id, label: ds.label, metrics: null };
      }
    }));

    const labels = [];
    const estCells = [];
    const medGenes = [];
    const meanReads = [];
    const fracInCells = [];
    const validBarcodes = [];
    const seqSat = [];

    for (const r of rows) {
      labels.push(r.label);
      if (r.metrics) {
        estCells.push(r.metrics.EstimatedNumberOfCells || 0);
        medGenes.push(r.metrics.MedianGenesPerCell || 0);
        meanReads.push(r.metrics.MeanReadsPerCell || 0);
        fracInCells.push(r.metrics.FractionReadsInCells || 0);
        validBarcodes.push(r.metrics.ValidBarcodes || 0);
        seqSat.push(r.metrics.SequencingSaturation || 0);
      } else {
        estCells.push(0); medGenes.push(0); meanReads.push(0); fracInCells.push(0); validBarcodes.push(0); seqSat.push(0);
      }
    }

    const traces = [
      { x: labels, y: estCells, type: 'bar', name: 'Estimated cells' },
      { x: labels, y: medGenes, type: 'bar', name: 'Median genes/cell' },
      { x: labels, y: meanReads, type: 'bar', name: 'Mean reads/cell' },
      { x: labels, y: fracInCells, type: 'bar', name: 'Fraction reads in cells' },
      { x: labels, y: validBarcodes, type: 'bar', name: 'Valid barcodes' },
      { x: labels, y: seqSat, type: 'bar', name: 'Sequencing saturation' }
    ];
    Plotly.react(qcPlot, traces, {
      barmode: 'group',
      title: 'QC summary (selected metrics)',
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: { color: '#e6ecff' },
      margin: { l: 60, r: 10, t: 40, b: 100 },
      xaxis: { tickangle: -30 }
    }, { displayModeBar: true, responsive: true });
  }

  async function renderClusterSizes(ds) {
    const clusters = await loadClusters(ds.paths.clusters);
    state.clusterList = clusters.unique;
    const counts = new Map();
    for (const cl of clusters.list) {
      counts.set(cl, (counts.get(cl) || 0) + 1);
    }
    const x = Array.from(counts.keys()).map(c => `Cluster ${c}`);
    const y = Array.from(counts.values());
    Plotly.react(clusterSizesPlot, [{ x, y, type: 'bar', marker: { color: '#4e79a7' } }], {
      title: `${ds.label} – Cluster sizes`,
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: '#e6ecff' }
    }, { displayModeBar: true, responsive: true });
  }

  async function renderPcaScree(ds) {
    const v = await loadPcaVariance(ds.paths.pcaVariance);
    const cum = [];
    let acc = 0;
    for (const p of v.proportion) { acc += p; cum.push(acc); }
    const traces = [
      { x: v.pc, y: v.proportion, type: 'bar', name: 'Proportion' },
      { x: v.pc, y: cum, type: 'scatter', mode: 'lines+markers', name: 'Cumulative' }
    ];
    Plotly.react(pcaScreePlot, traces, {
      title: `${ds.label} – PCA scree`,
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: '#e6ecff' }
    }, { displayModeBar: true, responsive: true });
  }

  async function renderPcaScatter(ds) {
    const [proj, cldata] = await Promise.all([
      loadPcaProjection(ds.paths.pcaProjection),
      loadClusters(ds.paths.clusters)
    ]);
    const palette = buildPalette(cldata.unique);
    const x = [], y = [], text = [], color = [];
    for (let i = 0; i < proj.barcode.length; i++) {
      const bc = proj.barcode[i];
      const cl = cldata.map.get(bc);
      x.push(proj.x[i]);
      y.push(proj.y[i]);
      text.push(`${bc} | cluster ${cl}`);
      color.push(palette.get(String(cl)) || '#cccccc');
    }
    Plotly.react(pcaScatterPlot, [{
      x, y, text, type: 'scattergl', mode: 'markers', marker: { size: 4, color, opacity: 0.85 }, hovertemplate: '%{text}<extra></extra>'
    }], {
      title: `${ds.label} – PCA PC1 vs PC2`, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: '#e6ecff' }, xaxis: { zeroline: false }, yaxis: { zeroline: false }
    }, { displayModeBar: true, responsive: true });
  }

  async function renderHvfTable(ds) {
    const features = await loadHvf(ds.paths.hvf);
    const top = features.slice(0, 100);
    hvfTable.innerHTML = '';
    const table = document.createElement('table');
    table.style.width = '100%';
    const thead = document.createElement('thead');
    thead.innerHTML = '<tr><th>#</th><th>Feature</th></tr>';
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    for (const row of top) {
      const tr = document.createElement('tr');
      const tdIdx = document.createElement('td'); tdIdx.textContent = String(row.idx);
      const tdFeat = document.createElement('td'); tdFeat.textContent = row.feature;
      tr.appendChild(tdIdx); tr.appendChild(tdFeat); tbody.appendChild(tr);
    }
    table.appendChild(tbody);
    hvfTable.appendChild(table);
  }

  async function populateDeClusters(ds) {
    const clusters = await loadClusters(ds.paths.clusters);
    deClusterSel.innerHTML = '';
    for (const cl of clusters.unique) {
      const opt = document.createElement('option');
      opt.value = String(cl);
      opt.textContent = `Cluster ${cl}`;
      deClusterSel.appendChild(opt);
    }
  }

  async function renderDifferentialExpression() {
    const ds = state.currentDataset;
    if (!ds) return;
    const clusterId = Number(deClusterSel.value);
    const de = await loadDiffExp(ds.paths.diffexp);
    // Build table for selected cluster
    const col = de.columnsForCluster(clusterId);
    const rows = de.rows
      .map(r => ({
        feature: r['Feature ID'] || r['Feature Name'] || r.feature || '',
        mean: toNum(r[col.mean]),
        log2fc: toNum(r[col.log2fc]),
        padj: toNum(r[col.padj])
      }))
      .filter(r => Number.isFinite(r.log2fc) && Number.isFinite(r.padj))
      .sort((a, b) => a.padj - b.padj)
      .slice(0, 200);

    // Table
    deTable.innerHTML = '';
    const table = document.createElement('table');
    table.style.width = '100%';
    const thead = document.createElement('thead');
    thead.innerHTML = '<tr><th>Feature</th><th>Mean</th><th>log2FC</th><th>adj p</th></tr>';
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    for (const r of rows) {
      const tr = document.createElement('tr');
      tr.innerHTML = `<td>${escapeHtml(r.feature)}</td><td>${fmt(r.mean)}</td><td>${fmt(r.log2fc)}</td><td>${fmt(r.padj)}</td>`;
      tbody.appendChild(tr);
    }
    table.appendChild(tbody);
    deTable.appendChild(table);

    // Volcano
    const x = rows.map(r => r.log2fc);
    const y = rows.map(r => -Math.log10(Math.max(r.padj, 1e-300)));
    Plotly.react(deVolcano, [{ x, y, mode: 'markers', type: 'scattergl', marker: { size: 5, color: '#e15759', opacity: 0.7 } }], {
      title: `${ds.label} – Cluster ${clusterId} volcano (log2FC vs -log10(q))`,
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: '#e6ecff' },
      xaxis: { title: 'log2 fold change' }, yaxis: { title: '-log10(adj p)' }
    }, { displayModeBar: true, responsive: true });
  }

  function renderMarkerTable() {
    const stage = markerStageSel.value; // 'gast' | 'Bla' | 'CPbla'
    const topN = markerTopNSel.value; // '5' | '25'
    const filename = `${stage}_top_specific_markers_top${topN}byMarkerScore.txt`;
    loadMarkerTable(filename).then(rows => {
      markersTable.innerHTML = '';
      const table = document.createElement('table');
      table.style.width = '100%';
      const thead = document.createElement('thead');
      thead.innerHTML = '<tr><th>gene_id</th><th>gene_short_name</th><th>cell_group</th><th>marker_score</th><th>mean_expression</th><th>fraction_expressing</th><th>specificity</th><th>marker_test_q_value</th></tr>';
      table.appendChild(thead);
      const tbody = document.createElement('tbody');
      for (const r of rows.slice(0, 200)) {
        const tr = document.createElement('tr');
        tr.innerHTML = `<td>${escapeHtml(r.gene_id || '')}</td><td>${escapeHtml(r.gene_short_name || '')}</td><td>${escapeHtml(r.cell_group || '')}</td><td>${fmt(r.marker_score)}</td><td>${fmt(r.mean_expression)}</td><td>${fmt(r.fraction_expressing)}</td><td>${fmt(r.specificity)}</td><td>${fmt(r.marker_test_q_value)}</td>`;
        tbody.appendChild(tr);
      }
      table.appendChild(tbody);
      markersTable.appendChild(table);
    }).catch(() => {
      markersTable.textContent = 'Marker table not available.';
    });
  }

  async function renderSmallMultiples() {
    smallMultiplesGrid.innerHTML = '';
    const datasets = state.datasets;
    // Limit to 8 for layout
    for (const ds of datasets) {
      const cell = document.createElement('div');
      cell.style.background = 'rgba(255,255,255,0.02)';
      cell.style.border = '1px solid rgba(255,255,255,0.08)';
      cell.style.borderRadius = '8px';
      cell.style.padding = '6px';
      const title = document.createElement('div');
      title.textContent = ds.label;
      title.style.margin = '4px 0 6px 6px';
      title.style.fontSize = '12px';
      const div = document.createElement('div');
      div.style.width = '100%';
      div.style.height = '240px';
      cell.appendChild(title);
      cell.appendChild(div);
      smallMultiplesGrid.appendChild(cell);
      try {
        const [coords, clusters] = await Promise.all([
          loadCoordinates(ds.paths.umap),
          loadClusters(ds.paths.clusters)
        ]);
        const palette = buildPalette(clusters.unique);
        const x = [], y = [], color = [];
        for (let i = 0; i < coords.barcode.length; i++) {
          const bc = coords.barcode[i];
          const cl = clusters.map.get(bc);
          x.push(coords.x[i]); y.push(coords.y[i]);
          color.push(palette.get(String(cl)) || '#cccccc');
        }
        Plotly.react(div, [{ x, y, type: 'scattergl', mode: 'markers', marker: { size: 3, color, opacity: 0.85 } }], {
          paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)', font: { color: '#e6ecff' }, margin: { l: 30, r: 5, t: 5, b: 30 }, xaxis: { showgrid: false, zeroline: false }, yaxis: { showgrid: false, zeroline: false }
        }, { displayModeBar: false, responsive: true });
      } catch (e) {
        div.textContent = 'Failed to load.';
      }
    }
  }

  // Data loaders
  async function fetchText(url) {
    const res = await fetch(url);
    if (!res.ok) throw new Error(`GET ${url} → ${res.status}`);
    return await res.text();
  }

  async function loadMetrics(url) {
    const txt = await fetchText(url);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: false });
    if (!parsed.data || parsed.data.length === 0) throw new Error('Empty metrics');
    const row = parsed.data[0];
    return {
      EstimatedNumberOfCells: toNum(row['Estimated Number of Cells']),
      MeanReadsPerCell: toNum(row['Mean Reads per Cell']),
      MedianGenesPerCell: toNum(row['Median Genes per Cell']),
      NumberOfReads: toNum(row['Number of Reads']),
      ValidBarcodes: toNum(row['Valid Barcodes']),
      SequencingSaturation: toNum(row['Sequencing Saturation']),
      FractionReadsInCells: toNum(row['Fraction Reads in Cells'])
    };
  }

  async function loadClusters(url) {
    const txt = await fetchText(url);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const map = new Map();
    const list = [];
    const unique = new Set();
    for (const row of parsed.data) {
      if (!row || row.Barcode == null) continue;
      map.set(String(row.Barcode), row.Cluster);
      list.push(row.Cluster);
      unique.add(row.Cluster);
    }
    return { map, list, unique: Array.from(unique).sort((a, b) => Number(a) - Number(b)) };
  }

  async function loadCoordinates(csvPath) {
    const txt = await fetchText(csvPath);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const x = [], y = [], barcode = [];
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

  async function loadPcaProjection(csvPath) {
    const txt = await fetchText(csvPath);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const x = [], y = [], barcode = [];
    for (const row of parsed.data) {
      if (!row || row.Barcode == null) continue;
      barcode.push(String(row.Barcode));
      x.push(Number(row['PC-1']));
      y.push(Number(row['PC-2']));
    }
    return { x, y, barcode };
  }

  async function loadPcaVariance(csvPath) {
    const txt = await fetchText(csvPath);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const pc = [], proportion = [];
    for (const row of parsed.data) {
      if (!row || row.PC == null) continue;
      pc.push(Number(row.PC));
      proportion.push(Number(row['Proportion.Variance.Explained']));
    }
    return { pc, proportion };
  }

  async function loadHvf(csvPath) {
    const txt = await fetchText(csvPath);
    const parsed = Papa.parse(txt, { header: false, dynamicTyping: false });
    // Expect header: Feature in first cell; subsequent lines like "1,LOC..."
    const out = [];
    for (let i = 1; i < parsed.data.length; i++) {
      const row = parsed.data[i];
      if (!row || row.length < 2) continue;
      out.push({ idx: toNum(row[0]), feature: String(row[1]) });
    }
    return out;
  }

  async function loadMarkerTable(filename) {
    const isPages = isPagesHost();
    const repoRawBase = getRepoRawBase();
    const url = isPages
      ? `${repoRawBase}/oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis/${filename}`
      : `../oyster_scRNASeq_jobs_genomic_resources_outs/Monocle_Rscripts_RDS/MarkerScoreAnalysis/${filename}`;
    const txt = await fetchText(url);
    const parsed = Papa.parse(txt, { header: true, delimiter: '\t', dynamicTyping: true });
    return parsed.data || [];
  }

  async function loadDiffExp(csvPath) {
    const txt = await fetchText(csvPath);
    const parsed = Papa.parse(txt, { header: true, dynamicTyping: true });
    const headers = parsed.meta.fields;
    function columnsForCluster(clusterId) {
      // Match like: 'Cluster 3 Mean Counts', 'Cluster 3 Log2 fold change', 'Cluster 3 Adjusted p value'
      const mean = headers.find(h => new RegExp(`^Cluster ${clusterId}\\s+Mean Counts$`).test(h));
      const log2fc = headers.find(h => new RegExp(`^Cluster ${clusterId}\\s+Log2 fold change$`).test(h));
      const padj = headers.find(h => new RegExp(`^Cluster ${clusterId}\\s+Adjusted p value$`).test(h));
      return { mean, log2fc, padj };
    }
    return { rows: parsed.data || [], columnsForCluster };
  }

  function buildPalette(uniqueClusters) {
    const colors = [
      '#4e79a7','#f28e2b','#e15759','#76b7b2','#59a14f','#edc949','#af7aa1','#ff9da7','#9c755f','#bab0ab',
      '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'
    ];
    const map = new Map();
    let idx = 0;
    for (const cl of uniqueClusters) {
      const key = String(cl);
      if (!map.has(key)) {
        map.set(key, colors[idx % colors.length]);
        idx++;
      }
    }
    return map;
  }

  // utils
  function toNum(v) {
    if (v == null || v === '') return NaN;
    if (typeof v === 'number') return v;
    const s = String(v).trim().replace(/%/g, '').replace(/,/g, '');
    const n = Number(s);
    if (!Number.isFinite(n)) return NaN;
    if (/[%]/.test(String(v))) return n / 100;
    return n;
  }

  function fmt(v) {
    if (!Number.isFinite(v)) return '';
    if (Math.abs(v) >= 1000) return v.toFixed(0);
    if (Math.abs(v) >= 1) return v.toFixed(3);
    return v.toExponential(2);
  }

  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;')
      .replace(/'/g, '&#039;');
  }
})();


