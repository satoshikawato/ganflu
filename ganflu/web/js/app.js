import { createMiniprotManager } from './app/miniprot.js';
import { createPyodideManager } from './app/pyodide.js';

const $ = (selector) => document.querySelector(selector);

const state = {
  outputs: {},
  selectedOutput: '',
  running: false
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

const renderOutputs = () => {
  renderTabs();
  const selected = state.selectedOutput;
  elements.outputText.value = selected ? state.outputs[selected] || '' : '';
  elements.downloadButton.disabled = !selected;
};

const renderSummary = (result = {}, gff3Text = '') => {
  const summary = result?.summary || {};
  const rows = [
    ['Records', summary.record_count ?? '-'],
    ['CDS', summary.cds_count ?? '-'],
    ['GFF3 lines', String(gff3Text || '').split(/\r?\n/).filter(Boolean).length],
    ['Record IDs', Array.isArray(summary.record_ids) ? summary.record_ids.join(', ') : '-']
  ];
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

elements.form.addEventListener('submit', async (event) => {
  event.preventDefault();
  if (state.running) return;

  setRunning(true);
  elements.logText.textContent = '';
  state.outputs = {};
  state.selectedOutput = '';
  renderSummary();
  renderOutputs();

  try {
    const target = elements.target.value;
    const isolate = elements.isolate.value.trim();
    if (!isolate) throw new Error('Isolate is required.');

    const { text: fastaText, name } = await readInputFasta();
    if (!fastaText.trim()) throw new Error('FASTA input is empty.');

    const outputStem = makeSafeStem(elements.outputStem.value, makeSafeStem(name));
    const reference = await getReferencePayload(target);
    const gff3Text = await miniprotManager.run({
      genomeFasta: fastaText,
      proteinFasta: reference.protein_fasta,
      prefix: 'MP',
      intronOpenPenalty: 15
    });
    const result = await runGff3ToOutputs({
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
    state.selectedOutput = Object.keys(state.outputs)[0] || '';
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
  elements.outputStem.value = '';
  state.outputs = {};
  state.selectedOutput = '';
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

setStatus('Idle');
renderSummary();
renderOutputs();
