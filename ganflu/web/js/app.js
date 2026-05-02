import { createMiniprotManager } from './app/miniprot.js';
import { createPyodideManager } from './app/pyodide.js';

const $ = (selector) => document.querySelector(selector);

const state = {
  outputs: {},
  selectedOutput: '',
  resultSummary: null,
  running: false
};

const AUTO_TARGETS = ['IAV', 'IBV', 'ICV', 'IDV'];
const AUTO_PREFIXES = {
  IAV: 'MPIA',
  IBV: 'MPIB',
  ICV: 'MPIC',
  IDV: 'MPID'
};

const elements = {
  form: $('#run-form'),
  fastaFile: $('#fasta-file'),
  fastaText: $('#fasta-text'),
  target: $('#target'),
  isolate: $('#isolate'),
  outputStem: $('#output-stem'),
  preserveOriginalId: $('#preserve-original-id'),
  runButton: $('#run-button'),
  status: $('#status'),
  summary: $('#summary'),
  outputTabs: $('#output-tabs'),
  autoDashboard: $('#auto-dashboard'),
  outputText: $('#output-text'),
  downloadButton: $('#download-button'),
  clearButton: $('#clear-button'),
  logText: $('#log-text')
};

const getStatusTone = (message) => {
  const normalized = String(message || '').toLowerCase();
  if (!normalized || normalized === 'idle') return 'idle';
  if (normalized.includes('error')) return 'error';
  if (normalized.includes('completed')) return 'success';
  return 'working';
};

const setStatus = (message) => {
  const value = message || 'Idle';
  elements.status.textContent = value;
  elements.status.dataset.status = getStatusTone(value);
};

const setRunning = (running) => {
  state.running = running;
  elements.runButton.disabled = running;
  elements.runButton.textContent = running ? 'Running...' : 'Run annotation';
  elements.form.classList.toggle('is-running', running);
};

const makeSafeStem = (value, fallback = 'ganflu') => {
  const cleaned = String(value || fallback)
    .replace(/\.[^.]+$/, '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};

const escapeHtml = (value) =>
  String(value ?? '').replace(/[&<>"']/g, (char) => ({
    '&': '&amp;',
    '<': '&lt;',
    '>': '&gt;',
    '"': '&quot;',
    "'": '&#39;'
  })[char]);

const readInputFasta = async () => {
  const file = elements.fastaFile.files?.[0];
  if (file) return { text: await file.text(), name: file.name };
  return { text: elements.fastaText.value, name: elements.outputStem.value || 'ganflu.fasta' };
};

const downloadTextFile = (filename, text) => {
  const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
};

const renderTabs = () => {
  elements.outputTabs.innerHTML = '';
  Object.keys(state.outputs).forEach((name) => {
    const button = document.createElement('button');
    button.type = 'button';
    button.className = name === state.selectedOutput ? 'tab is-active' : 'tab';
    button.textContent = name;
    button.addEventListener('click', () => {
      state.selectedOutput = name;
      renderOutputs();
    });
    elements.outputTabs.appendChild(button);
  });
};

const isSummaryOutput = (name) => /(?:\.auto)?\.summary\.json$/i.test(String(name || ''));

const parseJsonOutput = (text) => {
  try {
    return JSON.parse(String(text || ''));
  } catch {
    return null;
  }
};

const parseTsv = (text) => {
  const lines = String(text || '').trim().split(/\r?\n/).filter(Boolean);
  if (lines.length < 2) return [];
  const headers = lines[0].split('\t');
  return lines.slice(1).map((line) => {
    const values = line.split('\t');
    return headers.reduce((row, header, index) => {
      row[header] = values[index] ?? '';
      return row;
    }, {});
  });
};

const formatNumber = (value, digits = 3) => {
  const numberValue = Number(value);
  if (!Number.isFinite(numberValue)) return value ?? '-';
  return numberValue.toFixed(digits).replace(/\.?0+$/, '');
};

const getTsvRows = (pattern) => {
  const tsvName = Object.keys(state.outputs).find((name) => pattern.test(name));
  return tsvName ? parseTsv(state.outputs[tsvName]) : [];
};

const getAutoTsvRows = () => getTsvRows(/\.auto\.tsv$/i);

const getTargetHitRows = () => getTsvRows(/\.hits\.tsv$/i);

const renderTable = ({ columns, rows, emptyText }) => {
  if (!rows.length) {
    return `<div class="dashboard-empty">${escapeHtml(emptyText)}</div>`;
  }
  const headerHtml = columns.map((column) => `<th>${escapeHtml(column.label)}</th>`).join('');
  const rowsHtml = rows.map((row) => {
    const cells = columns.map((column) => {
      const value = typeof column.value === 'function' ? column.value(row) : row[column.value];
      const displayValue = value === '' || value == null ? '-' : value;
      return `<td>${escapeHtml(displayValue)}</td>`;
    }).join('');
    return `<tr>${cells}</tr>`;
  }).join('');
  return `<div class="dashboard-table-wrap"><table class="dashboard-table"><thead><tr>${headerHtml}</tr></thead><tbody>${rowsHtml}</tbody></table></div>`;
};

const renderAutoDashboard = (summary) => {
  const counts = summary?.counts || {};
  const byTarget = summary?.by_target || {};
  const byStatus = summary?.by_status || {};
  const thresholds = summary?.thresholds || {};
  const calls = getAutoTsvRows();
  const targetRows = Object.entries(byTarget).map(([target, values]) => ({
    target,
    accepted: values.accepted ?? 0,
    rejected: values.rejected ?? 0,
    complete: values.complete ?? 0,
    partial: values.partial ?? 0,
    segments: Object.entries(values.segments || {})
      .map(([segment, count]) => `${segment}: ${count}`)
      .join(', ') || '-'
  }));
  const statusRows = Object.entries(byStatus).map(([status, count]) => ({ status, count }));
  const thresholdRows = Object.entries(thresholds).map(([name, value]) => ({
    name: name.replaceAll('_', ' '),
    value: formatNumber(value)
  }));
  const cards = [
    ['Input contigs', counts.input_contigs ?? calls.length ?? '-'],
    ['Accepted', counts.accepted ?? '-'],
    ['Rejected', counts.rejected ?? '-']
  ].map(([label, value]) => (
    `<div class="dashboard-card"><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></div>`
  )).join('');

  elements.autoDashboard.innerHTML = `
    <div class="dashboard-grid">${cards}</div>
    <section class="dashboard-section">
      <h3>Targets</h3>
      ${renderTable({
        columns: [
          { label: 'Target', value: 'target' },
          { label: 'Accepted', value: 'accepted' },
          { label: 'Rejected', value: 'rejected' },
          { label: 'Complete', value: 'complete' },
          { label: 'Partial', value: 'partial' },
          { label: 'Segments', value: 'segments' }
        ],
        rows: targetRows,
        emptyText: 'No target calls.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Contig Calls</h3>
      ${renderTable({
        columns: [
          { label: 'Contig', value: 'contig_id' },
          { label: 'Call', value: 'call' },
          { label: 'Target', value: 'target' },
          { label: 'Segment', value: 'segment' },
          { label: 'Status', value: 'status' },
          { label: 'Score', value: (row) => formatNumber(row.best_score) },
          { label: 'Identity', value: (row) => formatNumber(row.best_identity) },
          { label: 'Coverage', value: (row) => formatNumber(row.best_aa_coverage) },
          { label: 'Flags', value: 'flags' }
        ],
        rows: calls,
        emptyText: 'No contig calls.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Statuses</h3>
      ${renderTable({
        columns: [
          { label: 'Status', value: 'status' },
          { label: 'Count', value: 'count' }
        ],
        rows: statusRows,
        emptyText: 'No statuses.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Thresholds</h3>
      ${renderTable({
        columns: [
          { label: 'Threshold', value: 'name' },
          { label: 'Value', value: 'value' }
        ],
        rows: thresholdRows,
        emptyText: 'No thresholds.'
      })}
    </section>
  `;
};

const renderTargetDashboard = (summary) => {
  const recordIds = Array.isArray(summary?.record_ids) ? summary.record_ids : [];
  const outputs = summary?.outputs || {};
  const bySegment = summary?.by_segment || {};
  const byStatus = summary?.by_status || {};
  const hits = getTargetHitRows();
  const outputRows = Object.entries(outputs).map(([name]) => ({ name }));
  const segmentRows = Object.entries(bySegment).map(([segment, count]) => ({ segment, count }));
  const statusRows = Object.entries(byStatus).map(([status, count]) => ({ status, count }));
  const cards = [
    ['Target', summary?.target || '-'],
    ['Input contigs', summary?.input_contigs ?? summary?.record_count ?? '-'],
    ['Annotated', summary?.annotated_contigs ?? summary?.record_count ?? '-'],
    ['Skipped', summary?.skipped_contigs ?? 0]
  ].map(([label, value]) => (
    `<div class="dashboard-card"><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></div>`
  )).join('');

  elements.autoDashboard.innerHTML = `
    <div class="dashboard-grid">${cards}</div>
    <section class="dashboard-section">
      <h3>Segments</h3>
      ${renderTable({
        columns: [
          { label: 'Segment', value: 'segment' },
          { label: 'Annotated contigs', value: 'count' }
        ],
        rows: segmentRows,
        emptyText: 'No segment-compatible hits.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Hit Calls</h3>
      ${renderTable({
        columns: [
          { label: 'Contig', value: 'contig_id' },
          { label: 'Call', value: 'call' },
          { label: 'Segment', value: 'segment' },
          { label: 'Status', value: 'status' },
          { label: 'Product', value: 'best_product' },
          { label: 'Identity', value: (row) => formatNumber(row.best_identity) },
          { label: 'Coverage', value: (row) => formatNumber(row.best_aa_coverage) },
          { label: 'Score', value: (row) => formatNumber(row.best_score) },
          { label: 'Query range', value: 'best_query_range' },
          { label: 'Flags', value: 'flags' }
        ],
        rows: hits,
        emptyText: 'No hit calls.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Statuses</h3>
      ${renderTable({
        columns: [
          { label: 'Status', value: 'status' },
          { label: 'Count', value: 'count' }
        ],
        rows: statusRows,
        emptyText: 'No statuses.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Record IDs</h3>
      ${renderTable({
        columns: [{ label: 'Record ID', value: 'recordId' }],
        rows: recordIds.map((recordId) => ({ recordId })),
        emptyText: 'No record IDs.'
      })}
    </section>
    <section class="dashboard-section">
      <h3>Outputs</h3>
      ${renderTable({
        columns: [{ label: 'File', value: 'name' }],
        rows: outputRows,
        emptyText: 'No output files.'
      })}
    </section>
  `;
};

const renderSummaryDashboard = (summary) => {
  if (summary?.auto) {
    renderAutoDashboard(summary);
  } else {
    renderTargetDashboard(summary);
  }
};

const renderOutputs = () => {
  renderTabs();
  const selected = state.selectedOutput;
  const isDashboard = isSummaryOutput(selected);
  const summary = isDashboard
    ? (parseJsonOutput(state.outputs[selected]) || state.resultSummary || {})
    : null;
  if (isDashboard && summary) {
    renderSummaryDashboard(summary);
  } else {
    elements.autoDashboard.innerHTML = '';
  }
  elements.autoDashboard.hidden = !isDashboard;
  elements.outputText.hidden = isDashboard;
  elements.outputText.value = !isDashboard && selected ? state.outputs[selected] || '' : '';
  elements.downloadButton.disabled = !selected;
};

const getPreferredOutputName = (outputs) => {
  const names = Object.keys(outputs || {});
  const summaryName = names.find(isSummaryOutput);
  return summaryName || names[0] || '';
};

const renderSummary = (result = {}, gff3Text = '') => {
  const summary = result?.summary || {};
  let rows;
  if (summary.auto) {
    const counts = summary.counts || {};
    const byTarget = summary.by_target || {};
    const targetText = Object.entries(byTarget)
      .map(([targetName, targetSummary]) => `${targetName}: ${targetSummary.accepted ?? 0}`)
      .join(', ') || '-';
    const statusText = Object.entries(summary.by_status || {})
      .map(([status, count]) => `${status}: ${count}`)
      .join(', ') || '-';
    rows = [
      ['Accepted', counts.accepted ?? '-'],
      ['Rejected', counts.rejected ?? '-'],
      ['Targets', targetText],
      ['Statuses', statusText]
    ];
  } else {
    rows = [
      ['Records', summary.record_count ?? '-'],
      ['CDS', summary.cds_count ?? '-'],
      ['GFF3 lines', summary.gff3_lines ?? String(gff3Text || '').split(/\r?\n/).filter(Boolean).length],
      ['Record IDs', Array.isArray(summary.record_ids) ? summary.record_ids.join(', ') : '-']
    ];
  }
  elements.summary.innerHTML = rows
    .map(([label, value]) => `<div><span>${escapeHtml(label)}</span><strong>${escapeHtml(value)}</strong></div>`)
    .join('');
};

const pyodideManager = createPyodideManager({ onStatus: setStatus });
const miniprotManager = createMiniprotManager({ onStatus: setStatus });

const getReferencePayload = async (target) => {
  const pyodide = await pyodideManager.init();
  pyodide.globals.set('GANFLU_TARGET', target);
  const payload = JSON.parse(pyodide.runPython('get_ganflu_reference_json(GANFLU_TARGET)'));
  if (payload.error) throw new Error(payload.error);
  return payload;
};

const runGff3ToOutputs = async ({
  fastaText,
  gff3Text,
  target,
  isolate,
  outputStem,
  preserveOriginalId
}) => {
  const pyodide = await pyodideManager.init();
  pyodide.globals.set('GANFLU_INPUT_FASTA', fastaText);
  pyodide.globals.set('GANFLU_GFF3_TEXT', gff3Text);
  pyodide.globals.set('GANFLU_TARGET', target);
  pyodide.globals.set('GANFLU_ISOLATE', isolate);
  pyodide.globals.set('GANFLU_OUTPUT_STEM', outputStem);
  pyodide.globals.set('GANFLU_PRESERVE_ORIGINAL_ID', Boolean(preserveOriginalId));
  setStatus('Converting GFF3 to GenBank');
  const result = JSON.parse(
    pyodide.runPython(
      'run_ganflu_web(GANFLU_INPUT_FASTA, GANFLU_GFF3_TEXT, GANFLU_TARGET, GANFLU_ISOLATE, GANFLU_OUTPUT_STEM, GANFLU_PRESERVE_ORIGINAL_ID)'
    )
  );
  if (result.error) {
    const detail = result.error.traceback || result.error.message || 'ganflu failed';
    throw new Error(detail);
  }
  return result;
};

const runAutoGff3ToOutputs = async ({
  fastaText,
  gff3ByTarget,
  isolate,
  outputStem,
  preserveOriginalId
}) => {
  const pyodide = await pyodideManager.init();
  pyodide.globals.set('GANFLU_INPUT_FASTA', fastaText);
  pyodide.globals.set('GANFLU_AUTO_GFF3_JSON', JSON.stringify(gff3ByTarget));
  pyodide.globals.set('GANFLU_ISOLATE', isolate);
  pyodide.globals.set('GANFLU_OUTPUT_STEM', outputStem);
  pyodide.globals.set('GANFLU_PRESERVE_ORIGINAL_ID', Boolean(preserveOriginalId));
  setStatus('Classifying contigs');
  const result = JSON.parse(
    pyodide.runPython(
      'run_ganflu_auto_web(GANFLU_INPUT_FASTA, GANFLU_AUTO_GFF3_JSON, GANFLU_ISOLATE, GANFLU_OUTPUT_STEM, GANFLU_PRESERVE_ORIGINAL_ID)'
    )
  );
  if (result.error) {
    const detail = result.error.traceback || result.error.message || 'ganflu auto mode failed';
    throw new Error(detail);
  }
  return result;
};

elements.form.addEventListener('submit', async (event) => {
  event.preventDefault();
  if (state.running) return;

  setRunning(true);
  elements.logText.textContent = '';
  state.outputs = {};
  state.selectedOutput = '';
  state.resultSummary = null;
  renderSummary();
  renderOutputs();

  try {
    const target = elements.target.value;

    const { text: fastaText, name } = await readInputFasta();
    if (!fastaText.trim()) throw new Error('FASTA input is empty.');

    const outputStem = makeSafeStem(elements.outputStem.value, makeSafeStem(name));
    const isolate = elements.isolate.value.trim() || outputStem;
    let result;
    let gff3Text = '';
    if (target === 'auto') {
      const gff3ByTarget = {};
      for (const autoTarget of AUTO_TARGETS) {
        const reference = await getReferencePayload(autoTarget);
        gff3ByTarget[autoTarget] = await miniprotManager.run({
          genomeFasta: fastaText,
          proteinFasta: reference.protein_fasta,
          prefix: AUTO_PREFIXES[autoTarget] || 'MP',
          intronOpenPenalty: 15,
          statusLabel: `Running Miniprot (${autoTarget})`
        });
      }
      result = await runAutoGff3ToOutputs({
        fastaText,
        gff3ByTarget,
        isolate,
        outputStem,
        preserveOriginalId: elements.preserveOriginalId.checked
      });
      state.outputs = result.outputs || {};
    } else {
      const reference = await getReferencePayload(target);
      gff3Text = await miniprotManager.run({
        genomeFasta: fastaText,
        proteinFasta: reference.protein_fasta,
        prefix: 'MP',
        intronOpenPenalty: 15
      });
      result = await runGff3ToOutputs({
        fastaText,
        gff3Text,
        target,
        isolate,
        outputStem,
        preserveOriginalId: elements.preserveOriginalId.checked
      });
      state.outputs = {
        [`${outputStem}.gff3`]: gff3Text,
        ...(result.outputs || {})
      };
    }
    state.resultSummary = result.summary || null;
    state.selectedOutput = getPreferredOutputName(state.outputs);
    renderSummary(result, gff3Text);
    renderOutputs();
    elements.logText.textContent = result.log || 'Completed.';
    setStatus('Completed');
  } catch (error) {
    const message = error?.message ? String(error.message) : String(error || 'Unknown error');
    elements.logText.textContent = message;
    setStatus('Error');
  } finally {
    setRunning(false);
  }
});

elements.downloadButton.addEventListener('click', () => {
  if (!state.selectedOutput) return;
  downloadTextFile(state.selectedOutput, state.outputs[state.selectedOutput] || '');
});

elements.clearButton.addEventListener('click', () => {
  elements.fastaFile.value = '';
  elements.fastaText.value = '';
  elements.target.value = 'auto';
  elements.outputStem.value = '';
  elements.isolate.value = '';
  state.outputs = {};
  state.selectedOutput = '';
  state.resultSummary = null;
  elements.logText.textContent = '';
  renderSummary();
  renderOutputs();
  setStatus('Idle');
});

elements.fastaFile.addEventListener('change', () => {
  const file = elements.fastaFile.files?.[0];
  if (file && !elements.outputStem.value.trim()) {
    elements.outputStem.value = makeSafeStem(file.name);
  }
});

elements.target.value = 'auto';
setStatus('Idle');
renderSummary();
renderOutputs();
