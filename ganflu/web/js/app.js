import { createMiniprotManager } from './app/miniprot.js';
import { createPyodideManager } from './app/pyodide.js';

const $ = (selector) => document.querySelector(selector);

const AUTO_TARGETS = ['IAV', 'IBV', 'ICV', 'IDV'];
const DOWNLOAD_ALL_VALUE = '__all__';
const UNCLASSIFIED_TARGET_KEY = '__unclassified__';
const UNCLASSIFIED_TARGET_LABEL = 'Unclassified';

const SAMPLE_PROFILES = {
  IAV: {
    target: 'IAV',
    outputStem: 'PR8',
    isolate: 'A/Puerto Rico/8/1934',
    files: ['./samples/IAV_PR8.fasta']
  },
  IBV: {
    target: 'IBV',
    outputStem: 'B_Victoria_2_1987',
    isolate: 'B/Victoria/2/1987',
    files: ['./samples/IBV_B_Victoria_2_1987.fa']
  },
  ICV: {
    target: 'ICV',
    outputStem: 'Ann_Arbor_1_1950',
    isolate: 'C/Ann Arbor/1/50',
    files: ['./samples/ICV_Ann_Arbor_1_1950.fna']
  },
  IDV: {
    target: 'IDV',
    outputStem: 'swine_Oklahoma_1334_2011',
    isolate: 'D/swine/Oklahoma/1334/2011',
    files: ['./samples/IDV_swine_Oklahoma_1334_2011.fna']
  },
  auto: {
    target: 'auto',
    outputStem: 'ganflu_sample_auto',
    isolate: 'ganflu_sample_auto',
    files: [
      './samples/IAV_PR8.fasta',
      './samples/IBV_B_Victoria_2_1987.fa',
      './samples/ICV_Ann_Arbor_1_1950.fna',
      './samples/IDV_swine_Oklahoma_1334_2011.fna'
    ]
  }
};

const state = {
  outputs: {},
  resultSummary: null,
  running: false,
  downloadSelection: DOWNLOAD_ALL_VALUE,
  selectedTargetTab: '',
  statusMessage: 'Idle',
  loadedSampleKey: '',
  sampleInputDirty: false
};

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
  advancedHitFilters: $('#advanced-hit-filters'),
  minIdentity: $('#min-identity'),
  minAaCoverage: $('#min-aa-coverage'),
  minScore: $('#min-score'),
  minMargin: $('#min-margin'),
  completeAaCoverage: $('#complete-aa-coverage'),
  maxSecondaryAlignments: $('#max-secondary-alignments'),
  outputScoreRatio: $('#output-score-ratio'),
  secondaryToPrimaryRatio: $('#secondary-to-primary-ratio'),
  preserveOriginalId: $('#preserve-original-id'),
  runButton: $('#run-button'),
  status: $('#status'),
  runSummary: $('#run-summary'),
  downloadSelect: $('#download-select'),
  downloadButton: $('#download-button'),
  workOverlay: $('#work-overlay'),
  workOverlayTitle: $('#work-overlay-title'),
  workOverlayMessage: $('#work-overlay-message'),
  advancedRunDetails: $('#advanced-run-details'),
  advancedRunBody: $('#advanced-run-body'),
  sampleTabs: document.querySelectorAll('[data-sample-key]'),
  clearButton: $('#clear-button'),
  logText: $('#log-text')
};

const HIT_SETTING_DEFINITIONS = [
  { key: 'minIdentity', element: elements.minIdentity, label: 'Minimum identity', min: 0, max: 1 },
  { key: 'minAaCoverage', element: elements.minAaCoverage, label: 'Minimum AA coverage', min: 0, max: 1 },
  { key: 'minScore', element: elements.minScore, label: 'Minimum score', min: 0, max: 1 },
  { key: 'minMargin', element: elements.minMargin, label: 'Minimum margin', min: 0, max: 1 },
  { key: 'completeAaCoverage', element: elements.completeAaCoverage, label: 'Complete coverage', min: 0, max: 1 },
  { key: 'maxSecondaryAlignments', element: elements.maxSecondaryAlignments, label: 'Max secondary alignments', min: 1, integer: true },
  { key: 'outputScoreRatio', element: elements.outputScoreRatio, label: 'Output score ratio', min: 0, max: 1 },
  { key: 'secondaryToPrimaryRatio', element: elements.secondaryToPrimaryRatio, label: 'Secondary/primary ratio', min: 0, max: 1 }
];

const getStatusTone = (message) => {
  const normalized = String(message || '').toLowerCase();
  if (!normalized || normalized === 'idle') return 'idle';
  if (normalized.includes('error')) return 'error';
  if (normalized.includes('completed') || normalized.includes('loaded')) return 'success';
  return 'working';
};

const setStatus = (message) => {
  const value = message || 'Idle';
  state.statusMessage = value;
  elements.status.textContent = value;
  elements.status.dataset.status = getStatusTone(value);
  updateWorkOverlay();
};

function getWorkOverlayTitle(message) {
  const normalized = String(message || '').toLowerCase();
  if (
    normalized.includes('loading') ||
    normalized.includes('installing') ||
    normalized.includes('runtime') ||
    normalized.includes('ready')
  ) {
    return 'Initializing ganflu...';
  }
  return 'Running annotation...';
}

function updateWorkOverlay() {
  if (!elements.workOverlay) return;
  const visible = Boolean(state.running);
  elements.workOverlay.hidden = !visible;
  elements.workOverlay.setAttribute('aria-busy', visible ? 'true' : 'false');
  if (!visible) return;
  const message = state.statusMessage && state.statusMessage !== 'Idle'
    ? state.statusMessage
    : 'Preparing annotation';
  if (elements.workOverlayTitle) {
    elements.workOverlayTitle.textContent = getWorkOverlayTitle(message);
  }
  if (elements.workOverlayMessage) {
    elements.workOverlayMessage.textContent = message;
  }
}

const setRunning = (running) => {
  state.running = running;
  elements.runButton.disabled = running;
  elements.runButton.textContent = running ? 'Running...' : 'Run annotation';
  elements.form.classList.toggle('is-running', running);
  updateWorkOverlay();
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

const downloadBlobFile = (filename, blob) => {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
};

const outputMimeType = (name) => {
  const lower = String(name || '').toLowerCase();
  if (lower.endsWith('.json')) return 'application/json;charset=utf-8';
  if (lower.endsWith('.tsv')) return 'text/tab-separated-values;charset=utf-8';
  if (lower.endsWith('.gff3') || lower.endsWith('.fna') || lower.endsWith('.faa') || lower.endsWith('.log')) {
    return 'text/plain;charset=utf-8';
  }
  if (lower.endsWith('.gbk')) return 'chemical/x-genbank;charset=utf-8';
  return 'text/plain;charset=utf-8';
};

const downloadTextFile = (filename, text) => {
  downloadBlobFile(filename, new Blob([String(text ?? '')], { type: outputMimeType(filename) }));
};

const zipEncoder = new TextEncoder();

const makeCrc32Table = () => {
  const table = new Uint32Array(256);
  for (let index = 0; index < 256; index += 1) {
    let value = index;
    for (let bit = 0; bit < 8; bit += 1) {
      value = (value & 1) ? (0xedb88320 ^ (value >>> 1)) : (value >>> 1);
    }
    table[index] = value >>> 0;
  }
  return table;
};

const crc32Table = makeCrc32Table();

const crc32 = (bytes) => {
  let crc = 0xffffffff;
  for (const byte of bytes) {
    crc = crc32Table[(crc ^ byte) & 0xff] ^ (crc >>> 8);
  }
  return (crc ^ 0xffffffff) >>> 0;
};

const getDosDateTime = (date = new Date()) => ({
  time: (
    (date.getHours() << 11) |
    (date.getMinutes() << 5) |
    Math.floor(date.getSeconds() / 2)
  ) & 0xffff,
  date: (
    ((date.getFullYear() - 1980) << 9) |
    ((date.getMonth() + 1) << 5) |
    date.getDate()
  ) & 0xffff
});

const sanitizeZipPath = (name) => {
  const path = String(name || 'output.txt')
    .replace(/\\/g, '/')
    .split('/')
    .filter((part) => part && part !== '.' && part !== '..')
    .join('/');
  return path || 'output.txt';
};

const uniquifyZipPath = (name, usedNames) => {
  const sanitized = sanitizeZipPath(name);
  if (!usedNames.has(sanitized)) {
    usedNames.add(sanitized);
    return sanitized;
  }
  const dotIndex = sanitized.lastIndexOf('.');
  const base = dotIndex > 0 ? sanitized.slice(0, dotIndex) : sanitized;
  const extension = dotIndex > 0 ? sanitized.slice(dotIndex) : '';
  let suffix = 2;
  let candidate = `${base}-${suffix}${extension}`;
  while (usedNames.has(candidate)) {
    suffix += 1;
    candidate = `${base}-${suffix}${extension}`;
  }
  usedNames.add(candidate);
  return candidate;
};

const makeZipHeader = (size) => {
  const bytes = new Uint8Array(size);
  return { bytes, view: new DataView(bytes.buffer) };
};

const createZipBlob = (files) => {
  const localChunks = [];
  const centralChunks = [];
  const usedNames = new Set();
  const { time, date } = getDosDateTime();
  let offset = 0;
  let centralSize = 0;

  files.forEach(({ name, text }) => {
    const filename = uniquifyZipPath(name, usedNames);
    const nameBytes = zipEncoder.encode(filename);
    const dataBytes = zipEncoder.encode(String(text ?? ''));
    const checksum = crc32(dataBytes);
    const local = makeZipHeader(30 + nameBytes.length);
    local.view.setUint32(0, 0x04034b50, true);
    local.view.setUint16(4, 20, true);
    local.view.setUint16(6, 0x0800, true);
    local.view.setUint16(8, 0, true);
    local.view.setUint16(10, time, true);
    local.view.setUint16(12, date, true);
    local.view.setUint32(14, checksum, true);
    local.view.setUint32(18, dataBytes.length, true);
    local.view.setUint32(22, dataBytes.length, true);
    local.view.setUint16(26, nameBytes.length, true);
    local.view.setUint16(28, 0, true);
    local.bytes.set(nameBytes, 30);

    const central = makeZipHeader(46 + nameBytes.length);
    central.view.setUint32(0, 0x02014b50, true);
    central.view.setUint16(4, 20, true);
    central.view.setUint16(6, 20, true);
    central.view.setUint16(8, 0x0800, true);
    central.view.setUint16(10, 0, true);
    central.view.setUint16(12, time, true);
    central.view.setUint16(14, date, true);
    central.view.setUint32(16, checksum, true);
    central.view.setUint32(20, dataBytes.length, true);
    central.view.setUint32(24, dataBytes.length, true);
    central.view.setUint16(28, nameBytes.length, true);
    central.view.setUint16(30, 0, true);
    central.view.setUint16(32, 0, true);
    central.view.setUint16(34, 0, true);
    central.view.setUint16(36, 0, true);
    central.view.setUint32(38, 0, true);
    central.view.setUint32(42, offset, true);
    central.bytes.set(nameBytes, 46);

    localChunks.push(local.bytes, dataBytes);
    centralChunks.push(central.bytes);
    offset += local.bytes.length + dataBytes.length;
    centralSize += central.bytes.length;
  });

  const end = makeZipHeader(22);
  end.view.setUint32(0, 0x06054b50, true);
  end.view.setUint16(4, 0, true);
  end.view.setUint16(6, 0, true);
  end.view.setUint16(8, files.length, true);
  end.view.setUint16(10, files.length, true);
  end.view.setUint32(12, centralSize, true);
  end.view.setUint32(16, offset, true);
  end.view.setUint16(20, 0, true);

  return new Blob([...localChunks, ...centralChunks, end.bytes], { type: 'application/zip' });
};

const getResultArchiveName = () => {
  const stem = state.resultSummary?.run?.output_stem ||
    state.resultSummary?.output_stem ||
    elements.outputStem.value ||
    'ganflu-results';
  return `${makeSafeStem(stem, 'ganflu-results')}.results.zip`;
};

const downloadResultsZip = () => {
  const downloads = getDownloads();
  if (!downloads.length) return;
  const files = downloads.map((download) => ({
    name: download.name,
    text: state.outputs[download.name] ?? ''
  }));
  downloadBlobFile(getResultArchiveName(), createZipBlob(files));
};

const readNumericSetting = ({ element, label, min = -Infinity, max = Infinity, integer = false }) => {
  const fallback = Number(element?.dataset?.default);
  const rawValue = String(element?.value ?? '').trim();
  const value = rawValue === '' ? fallback : Number(rawValue);
  if (!Number.isFinite(value)) {
    throw new Error(`${label} must be a number.`);
  }
  if (value < min && max === Infinity) {
    throw new Error(`${label} must be at least ${min}.`);
  }
  if (value < min || value > max) {
    throw new Error(`${label} must be between ${min} and ${max}.`);
  }
  if (integer && !Number.isInteger(value)) {
    throw new Error(`${label} must be a whole number.`);
  }
  return integer ? Math.trunc(value) : value;
};

const getHitSettings = () => HIT_SETTING_DEFINITIONS.reduce((settings, definition) => {
  settings[definition.key] = readNumericSetting(definition);
  return settings;
}, {});

const resetHitSettings = () => {
  HIT_SETTING_DEFINITIONS.forEach(({ element }) => {
    if (element) element.value = element.dataset.default ?? '';
  });
  if (elements.advancedHitFilters) elements.advancedHitFilters.open = false;
};

const sampleTextCache = new Map();

const renderSampleTabs = () => {
  elements.sampleTabs.forEach((button) => {
    const active = Boolean(state.loadedSampleKey && !state.sampleInputDirty && button.dataset.sampleKey === state.loadedSampleKey);
    button.classList.toggle('is-active', active);
    button.setAttribute('aria-selected', active ? 'true' : 'false');
  });
};

const resetResultState = () => {
  state.outputs = {};
  state.resultSummary = null;
  state.downloadSelection = DOWNLOAD_ALL_VALUE;
  state.selectedTargetTab = '';
  elements.logText.textContent = '';
  renderResults();
};

const fetchSampleText = async (path) => {
  if (!sampleTextCache.has(path)) {
    const response = await fetch(path, { cache: 'no-cache' });
    if (!response.ok) {
      throw new Error(`Unable to load sample FASTA: ${path}`);
    }
    sampleTextCache.set(path, await response.text());
  }
  return sampleTextCache.get(path);
};

const getSampleProfile = (key) => SAMPLE_PROFILES[key] || SAMPLE_PROFILES.auto;

const loadSampleProfile = async (key = elements.target.value, { announce = false } = {}) => {
  const profile = getSampleProfile(key);
  setStatus('Loading sample data');
  elements.logText.textContent = '';
  elements.sampleTabs.forEach((button) => button.setAttribute('aria-busy', 'true'));
  try {
    const fastaParts = await Promise.all(profile.files.map(fetchSampleText));
    elements.fastaFile.value = '';
    elements.fastaText.value = `${fastaParts.map((text) => text.trim()).filter(Boolean).join('\n\n')}\n`;
    elements.target.value = profile.target;
    elements.outputStem.value = profile.outputStem;
    elements.isolate.value = profile.isolate;
    elements.preserveOriginalId.checked = false;
    resetHitSettings();
    state.loadedSampleKey = profile.target;
    state.sampleInputDirty = false;
    renderSampleTabs();
    resetResultState();
    setStatus(announce ? 'Sample loaded' : 'Idle');
  } finally {
    elements.sampleTabs.forEach((button) => button.removeAttribute('aria-busy'));
  }
};

const handleSampleLoadError = (error) => {
  const message = error?.message ? String(error.message) : String(error || 'Sample load failed');
  state.loadedSampleKey = '';
  renderSampleTabs();
  elements.runSummary.innerHTML = `<div class="run-summary-empty">Sample load failed: ${escapeHtml(message)}</div>`;
  elements.logText.textContent = message;
  setStatus('Sample load error');
};

const markSampleInputDirty = () => {
  state.sampleInputDirty = true;
  state.loadedSampleKey = '';
  renderSampleTabs();
};

const formatNumber = (value, digits = 3) => {
  const numberValue = Number(value);
  if (!Number.isFinite(numberValue)) return value ?? '-';
  if (digits === 0) return String(Math.round(numberValue));
  return numberValue.toFixed(digits).replace(/\.?0+$/, '');
};

const renderMetric = (value, digits = 3) => {
  if (value === '' || value == null) return '-';
  return formatNumber(value, digits);
};

const renderPercentMetric = (value, digits = 1) => {
  if (value === '' || value == null) return '-';
  const numberValue = Number(value);
  if (!Number.isFinite(numberValue)) return value ?? '-';
  return `${formatNumber(numberValue * 100, digits)}%`;
};

const asList = (value) => {
  if (Array.isArray(value)) return value.filter((item) => item != null && item !== '');
  if (value == null || value === '') return [];
  return [value];
};

const inferDownloadKind = (name) => {
  const lower = String(name || '').toLowerCase();
  if (lower.endsWith('.summary.json')) return 'summary_json';
  if (lower.endsWith('.gff3')) return 'gff3';
  if (lower.endsWith('.gbk')) return 'genbank';
  if (lower.endsWith('.cds.fna')) return 'cds_fasta';
  if (lower.endsWith('.faa')) return 'protein_fasta';
  if (lower.endsWith('.hits.tsv') || lower.endsWith('.auto.tsv')) return 'contig_report';
  if (lower.endsWith('.log')) return 'log';
  return 'output';
};

const DOWNLOAD_KIND_LABELS = {
  summary_json: 'Summary JSON',
  gff3: 'GFF3',
  genbank: 'GenBank',
  cds_fasta: 'CDS FASTA',
  protein_fasta: 'Protein FASTA',
  contig_report: 'Contig report',
  log: 'Run log',
  output: 'Output'
};

const getDownloads = (summary = state.resultSummary) => {
  const downloads = Array.isArray(summary?.downloads) ? summary.downloads : [];
  const seen = new Set();
  if (downloads.length) {
    const manifestDownloads = downloads.filter((download) => {
      const available = Object.prototype.hasOwnProperty.call(state.outputs, download.name);
      if (available) seen.add(download.name);
      return available;
    });
    const extraDownloads = Object.keys(state.outputs)
      .filter((name) => !seen.has(name))
      .map((name) => ({
        name,
        kind: inferDownloadKind(name),
        target: ''
      }));
    return [...manifestDownloads, ...extraDownloads];
  }
  return Object.keys(state.outputs).map((name) => ({
    name,
    kind: inferDownloadKind(name),
    target: ''
  }));
};

const formatDownloadLabel = (download) => {
  const kind = DOWNLOAD_KIND_LABELS[download.kind] || DOWNLOAD_KIND_LABELS.output;
  const target = String(download.target || '').trim();
  const prefix = target ? `${target} / ` : '';
  return `${prefix}${kind} - ${download.name}`;
};

const appendDownloadOption = (value, label) => {
  const option = document.createElement('option');
  option.value = value;
  option.textContent = label;
  elements.downloadSelect.appendChild(option);
};

const renderDownloadOptions = (summary = state.resultSummary) => {
  const downloads = getDownloads(summary);
  const hasDownloads = downloads.length > 0;
  elements.downloadSelect.innerHTML = '';
  if (hasDownloads) {
    appendDownloadOption(DOWNLOAD_ALL_VALUE, `All results (.zip) - ${downloads.length} files`);
    downloads.forEach((download) => appendDownloadOption(download.name, formatDownloadLabel(download)));
  } else {
    appendDownloadOption('', 'No outputs');
  }
  const validSelections = new Set([DOWNLOAD_ALL_VALUE, ...downloads.map((download) => download.name)]);
  if (!validSelections.has(state.downloadSelection)) {
    state.downloadSelection = DOWNLOAD_ALL_VALUE;
  }
  elements.downloadSelect.disabled = !hasDownloads;
  elements.downloadSelect.value = hasDownloads ? state.downloadSelection : '';
  elements.downloadButton.disabled = !hasDownloads;
  const selectedDownload = downloads.find((download) => download.name === state.downloadSelection);
  elements.downloadButton.title = !hasDownloads
    ? 'No output files to download.'
    : selectedDownload
      ? `Download ${selectedDownload.name}.`
      : `Download ${downloads.length} output files as a ZIP archive.`;
};

const downloadSelectedResult = () => {
  const downloads = getDownloads();
  if (!downloads.length) return;
  const selected = elements.downloadSelect.value || DOWNLOAD_ALL_VALUE;
  if (selected === DOWNLOAD_ALL_VALUE) {
    downloadResultsZip();
    return;
  }
  const download = downloads.find((item) => item.name === selected);
  if (!download) {
    throw new Error('Selected output file is unavailable.');
  }
  downloadTextFile(download.name, state.outputs[download.name]);
};

const renderBadge = (value, tone = '') => {
  if (!value) return '';
  const className = tone ? `summary-badge ${tone}` : 'summary-badge';
  return `<span class="${className}">${escapeHtml(value)}</span>`;
};

const metricValue = (feature, contig, key) => {
  const value = feature?.metrics?.[key];
  if (value !== '' && value != null) return value;
  return contig?.best_hit?.[key];
};

const formatStrandLabel = (strand) => {
  if (strand === '+') return 'Forward (+)';
  if (strand === '-') return 'Reverse (-)';
  return strand || '-';
};

const sequenceLengthLabel = (sequence, unit) => {
  const length = String(sequence || '').length;
  return length ? `${formatNumber(length, 0)} ${unit}` : `No ${unit}`;
};

const normalizeSequence = (sequence) => String(sequence || '').replace(/\s+/g, '').toUpperCase();

const wrapSequence = (sequence, width = 60) => {
  const normalized = normalizeSequence(sequence);
  const lines = [];
  for (let index = 0; index < normalized.length; index += width) {
    lines.push(normalized.slice(index, index + width));
  }
  return lines.join('\n');
};

const fastaId = (feature, contig, sequenceType) => {
  const fallback = [
    contig?.record_id || contig?.input_id,
    feature?.gene || feature?.product || feature?.type,
    sequenceType
  ].filter(Boolean).join('_');
  return String(feature?.id || fallback || sequenceType || 'sequence')
    .trim()
    .replace(/\s+/g, '_')
    .replace(/[^A-Za-z0-9_.:-]+/g, '_') || 'sequence';
};

const fastaHeaderValue = (value) => String(value || '').replace(/\s+/g, ' ').trim();

const formatSequenceFasta = (feature, contig, sequenceType) => {
  const sequence = normalizeSequence(feature?.sequences?.[sequenceType]);
  if (!sequence) return '';
  const label = sequenceType === 'aa' ? 'AA' : 'CDS_nt';
  const attributes = [
    label,
    feature?.gene ? `gene=${fastaHeaderValue(feature.gene)}` : '',
    feature?.product ? `product=${fastaHeaderValue(feature.product)}` : '',
    feature?.location ? `location=${fastaHeaderValue(feature.location)}` : ''
  ].filter(Boolean).join(' ');
  return `>${fastaId(feature, contig, sequenceType)} ${attributes}\n${wrapSequence(sequence)}\n`;
};

const renderFeatureQualifiers = (notes) => {
  if (!notes.length) return '';
  const rows = notes.map((note) => {
    const isSlippage = note === 'ribosomal slippage';
    const qualifier = isSlippage ? 'GenBank /ribosomal_slippage' : 'GenBank /note';
    const value = isSlippage ? 'present' : note;
    return `<div><dt>${escapeHtml(qualifier)}</dt><dd>${escapeHtml(value)}</dd></div>`;
  }).join('');
  return `<dl class="feature-qualifiers">${rows}</dl>`;
};

const renderSequenceBox = ({
  title,
  copyLabel,
  sequenceType,
  fasta,
  lengthLabel,
  contigIndex,
  featureIndex
}) => `
  <section class="sequence-box${fasta ? '' : ' is-empty'}">
    <div class="sequence-box-head">
      <div>
        <h5>${escapeHtml(title)}</h5>
        <span>${escapeHtml(lengthLabel)}</span>
      </div>
      ${fasta ? `<button class="sequence-copy" type="button" data-copy-label="${escapeHtml(copyLabel)}" data-contig-index="${contigIndex}" data-feature-index="${featureIndex}" data-sequence="${sequenceType}">${escapeHtml(copyLabel)}</button>` : ''}
    </div>
    <pre>${escapeHtml(fasta || 'No sequence')}</pre>
  </section>
`;

const renderFeatureSummary = (feature, contig, contigIndex, featureIndex) => {
  const reference = feature.reference || {};
  const sequences = feature.sequences || {};
  const cdsNt = sequences.cds_nt || '';
  const aa = sequences.aa || '';
  const notes = asList(feature.notes);
  const refProduct = reference.product || contig.best_hit?.product || '-';
  const targetRange = reference.target_range || feature.target_range || contig.best_hit?.target_range || '-';
  const queryRange = feature.query_range || contig.best_hit?.query_range || '-';
  const featureTitle = [feature.gene, feature.product].filter(Boolean).join(' / ') || feature.id || 'Feature';
  const identity = renderPercentMetric(metricValue(feature, contig, 'identity'));
  const aaCoverage = renderPercentMetric(metricValue(feature, contig, 'aa_coverage'));
  const cdsFasta = formatSequenceFasta(feature, contig, 'cds_nt');
  const aaFasta = formatSequenceFasta(feature, contig, 'aa');

  return `
    <details class="feature-summary">
      <summary class="feature-title-row">
        <div class="feature-title-main">
          <h4>${escapeHtml(featureTitle)}</h4>
          <p>Location ${escapeHtml(feature.location || '-')} | Identity ${escapeHtml(identity)} | AA coverage ${escapeHtml(aaCoverage)}</p>
        </div>
        <div class="feature-badges">
          ${renderBadge(feature.type || 'feature')}
          <span class="feature-toggle-icon" aria-hidden="true"></span>
        </div>
      </summary>
      <div class="feature-body">
        <dl class="feature-grid">
          <div><dt>Location</dt><dd>${escapeHtml(feature.location || '-')}</dd></div>
          <div><dt>Reference</dt><dd>${escapeHtml(refProduct)}</dd></div>
          <div><dt>Strand</dt><dd>${escapeHtml(formatStrandLabel(feature.strand))}</dd></div>
          <div><dt>Identity</dt><dd>${escapeHtml(identity)}</dd></div>
          <div><dt>AA coverage</dt><dd>${escapeHtml(aaCoverage)}</dd></div>
          <div><dt>Score</dt><dd>${escapeHtml(renderMetric(metricValue(feature, contig, 'score')))}</dd></div>
          <div><dt>Query range</dt><dd>${escapeHtml(queryRange)}</dd></div>
          <div><dt>Target range</dt><dd>${escapeHtml(targetRange)}</dd></div>
        </dl>
        ${renderFeatureQualifiers(notes)}
        <div class="sequence-boxes">
          ${renderSequenceBox({
            title: 'CDS nucleotide FASTA',
            copyLabel: 'Copy CDS',
            sequenceType: 'cds_nt',
            fasta: cdsFasta,
            lengthLabel: sequenceLengthLabel(cdsNt, 'nt'),
            contigIndex,
            featureIndex
          })}
          ${renderSequenceBox({
            title: 'Amino acid FASTA',
            copyLabel: 'Copy AA',
            sequenceType: 'aa',
            fasta: aaFasta,
            lengthLabel: sequenceLengthLabel(aa, 'aa'),
            contigIndex,
            featureIndex
          })}
        </div>
      </div>
    </details>
  `;
};

const renderContigSummary = (contig, contigIndex) => {
  const features = Array.isArray(contig.features) ? contig.features : [];
  const notes = [...asList(contig.flags), ...asList(contig.notes)];
  const qcTone = contig.qc_result === 'pass' ? 'pass' : contig.qc_result === 'fail' ? 'fail' : '';
  const featureHtml = features.length
    ? features.map((feature, featureIndex) => renderFeatureSummary(feature, contig, contigIndex, featureIndex)).join('')
    : '<div class="run-summary-empty compact">No annotated features.</div>';

  return `
    <section class="contig-summary">
      <div class="contig-head">
        <div class="contig-head-main">
          <div class="contig-title-main">
            <h3>${escapeHtml(contig.input_id || contig.record_id || 'Contig')}</h3>
            <p>${escapeHtml(contig.record_id || 'No GenBank record')}</p>
          </div>
          <dl class="contig-head-metrics">
            <div><dt>Length</dt><dd>${escapeHtml(formatNumber(contig.length, 0))}</dd></div>
            <div><dt>Call</dt><dd>${escapeHtml(contig.call || '-')}</dd></div>
            <div><dt>Confidence</dt><dd>${escapeHtml(contig.confidence || '-')}</dd></div>
          </dl>
        </div>
        <div class="contig-badges">
          ${renderBadge(contig.target || 'unassigned')}
          ${renderBadge(contig.segment || 'no segment')}
          ${renderBadge(contig.status || 'unknown')}
          ${renderBadge(contig.qc_result || 'qc unknown', qcTone)}
        </div>
      </div>
      ${notes.length ? `<div class="contig-notes">${notes.map((note) => escapeHtml(note)).join('; ')}</div>` : ''}
      <div class="feature-list">${featureHtml}</div>
    </section>
  `;
};

const isAutoSummary = (summary) => {
  const mode = String(summary?.mode || '').toLowerCase();
  const target = String(summary?.run?.target || summary?.target || '').toLowerCase();
  return Boolean(summary?.auto) || mode === 'auto' || target === 'auto';
};

const displayTargetLabel = (target) => (
  target === UNCLASSIFIED_TARGET_KEY ? UNCLASSIFIED_TARGET_LABEL : (target || 'Unassigned')
);

const getAutoContigTabKey = (contig) => {
  const target = String(contig?.target || '').trim();
  const call = String(contig?.call || '').toLowerCase();
  return call === 'accept' && target ? target : UNCLASSIFIED_TARGET_KEY;
};

const groupAutoContigsByTab = (contigs) => {
  const groupsByTarget = new Map();
  contigs.forEach((contig, index) => {
    const target = getAutoContigTabKey(contig);
    if (!groupsByTarget.has(target)) {
      groupsByTarget.set(target, []);
    }
    groupsByTarget.get(target).push({ contig, index });
  });

  const orderedTargets = [
    ...AUTO_TARGETS.filter((target) => groupsByTarget.has(target)),
    ...Array.from(groupsByTarget.keys())
      .filter((target) => !AUTO_TARGETS.includes(target) && target !== UNCLASSIFIED_TARGET_KEY)
      .sort((left, right) => left.localeCompare(right)),
    ...(groupsByTarget.has(UNCLASSIFIED_TARGET_KEY) ? [UNCLASSIFIED_TARGET_KEY] : [])
  ];

  return orderedTargets.map((target) => ({
    target,
    items: groupsByTarget.get(target)
  }));
};

const contigCountLabel = (count) => `${formatNumber(count, 0)} contig${count === 1 ? '' : 's'}`;

const getActiveTargetTab = (groups) => {
  const keys = groups.map((group) => group.target);
  if (!keys.includes(state.selectedTargetTab)) {
    state.selectedTargetTab = keys[0] || '';
  }
  return state.selectedTargetTab;
};

const renderAutoTargetTabs = (groups) => {
  const activeTarget = getActiveTargetTab(groups);
  const activeGroup = groups.find((group) => group.target === activeTarget) || groups[0];
  const tabs = groups.map(({ target, items }) => {
    const active = target === activeTarget;
    return `
      <button class="target-tab${active ? ' is-active' : ''}" type="button" role="tab" aria-selected="${active ? 'true' : 'false'}" data-target-tab="${escapeHtml(target)}">
        ${escapeHtml(displayTargetLabel(target))}<span>${escapeHtml(String(items.length))}</span>
      </button>
    `;
  }).join('');

  return `
    <div class="target-tab-view">
      <div class="target-tabs" role="tablist" aria-label="Auto result targets">
        ${tabs}
      </div>
      <section class="target-tab-panel" role="tabpanel" aria-label="${escapeHtml(displayTargetLabel(activeGroup.target))}">
        <div class="target-contig-group-head">
          <h3>${escapeHtml(displayTargetLabel(activeGroup.target))}</h3>
          <span>${escapeHtml(contigCountLabel(activeGroup.items.length))}</span>
        </div>
        <div class="target-contig-items">
          ${activeGroup.items.map(({ contig, index }) => renderContigSummary(contig, index)).join('')}
        </div>
      </section>
    </div>
  `;
};

const renderContigList = (summary, contigs) => {
  if (!contigs.length) return '<div class="run-summary-empty">No contig results.</div>';
  if (!isAutoSummary(summary)) return contigs.map(renderContigSummary).join('');
  return renderAutoTargetTabs(groupAutoContigsByTab(contigs));
};

const renderRunSummary = (summary = state.resultSummary) => {
  if (!summary) {
    elements.runSummary.innerHTML = '<div class="run-summary-empty">No run summary.</div>';
    return;
  }
  const run = summary.run || {};
  const counts = run.counts || summary.counts || {};
  const contigs = Array.isArray(summary.contigs) ? summary.contigs : [];
  const targetText = run.target || summary.target || (summary.auto ? 'auto' : '-');
  const countItems = [
    ['Input contigs', counts.input_contigs],
    ['Annotated', counts.annotated_contigs ?? counts.accepted],
    ['Passed', counts.passed],
    ['Failed', counts.failed]
  ];
  const overview = countItems.map(([label, value]) => (
    `<div><dt>${escapeHtml(label)}</dt><dd>${escapeHtml(value ?? '-')}</dd></div>`
  )).join('');

  elements.runSummary.innerHTML = `
    <section class="run-overview">
      <dl>
        <div><dt>Mode</dt><dd>${escapeHtml(summary.mode || (summary.auto ? 'auto' : 'target'))}</dd></div>
        <div><dt>Target</dt><dd>${escapeHtml(targetText)}</dd></div>
        <div><dt>Isolate</dt><dd>${escapeHtml(run.isolate || summary.isolate || '-')}</dd></div>
        <div><dt>Output stem</dt><dd>${escapeHtml(run.output_stem || summary.output_stem || '-')}</dd></div>
        ${overview}
      </dl>
    </section>
    <div class="contig-list">
      ${renderContigList(summary, contigs)}
    </div>
  `;
};

const renderSettingsTable = (title, settings) => {
  const rows = Object.entries(settings || {});
  if (!rows.length) return '';
  return `
    <section class="advanced-block">
      <h3>${escapeHtml(title)}</h3>
      <table class="advanced-table">
        <tbody>
          ${rows.map(([name, value]) => (
            `<tr><th>${escapeHtml(name.replaceAll('_', ' '))}</th><td>${escapeHtml(formatNumber(value))}</td></tr>`
          )).join('')}
        </tbody>
      </table>
    </section>
  `;
};

const renderAdvancedDetails = (summary = state.resultSummary) => {
  const runLog = summary?.run_log || {};
  elements.advancedRunDetails.hidden = !summary;
  elements.advancedRunBody.innerHTML = summary
    ? [
      renderSettingsTable('Thresholds', runLog.thresholds || summary.thresholds),
      renderSettingsTable('Miniprot', runLog.miniprot || summary.miniprot)
    ].join('') || '<div class="run-summary-empty compact">No advanced settings.</div>'
    : '';
  elements.logText.textContent = runLog.messages || '';
};

const renderResults = () => {
  renderDownloadOptions();
  renderRunSummary();
  renderAdvancedDetails();
};

const writeClipboardText = async (text) => {
  if (navigator.clipboard?.writeText) {
    await navigator.clipboard.writeText(text);
    return;
  }
  const textarea = document.createElement('textarea');
  textarea.value = text;
  textarea.setAttribute('readonly', '');
  textarea.style.position = 'fixed';
  textarea.style.left = '-9999px';
  document.body.appendChild(textarea);
  textarea.select();
  document.execCommand('copy');
  textarea.remove();
};

const copyResetTimers = new WeakMap();

const copySequence = async (button) => {
  const contigIndex = Number(button.dataset.contigIndex);
  const featureIndex = Number(button.dataset.featureIndex);
  const sequenceType = button.dataset.sequence;
  const contig = state.resultSummary?.contigs?.[contigIndex];
  const feature = contig?.features?.[featureIndex];
  const text = formatSequenceFasta(feature, contig, sequenceType);
  if (!text) return;
  const originalText = button.dataset.copyLabel || button.textContent;
  await writeClipboardText(text);
  button.textContent = 'Copied!';
  window.clearTimeout(copyResetTimers.get(button));
  copyResetTimers.set(button, window.setTimeout(() => {
    button.textContent = originalText;
  }, 2000));
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
  preserveOriginalId,
  hitSettings
}) => {
  const pyodide = await pyodideManager.init();
  pyodide.globals.set('GANFLU_INPUT_FASTA', fastaText);
  pyodide.globals.set('GANFLU_GFF3_TEXT', gff3Text);
  pyodide.globals.set('GANFLU_TARGET', target);
  pyodide.globals.set('GANFLU_ISOLATE', isolate);
  pyodide.globals.set('GANFLU_OUTPUT_STEM', outputStem);
  pyodide.globals.set('GANFLU_PRESERVE_ORIGINAL_ID', Boolean(preserveOriginalId));
  pyodide.globals.set('GANFLU_HIT_SETTINGS_JSON', JSON.stringify(hitSettings || {}));
  setStatus('Converting GFF3 to GenBank');
  const result = JSON.parse(
    pyodide.runPython(
      'run_ganflu_web(GANFLU_INPUT_FASTA, GANFLU_GFF3_TEXT, GANFLU_TARGET, GANFLU_ISOLATE, GANFLU_OUTPUT_STEM, GANFLU_PRESERVE_ORIGINAL_ID, GANFLU_HIT_SETTINGS_JSON)'
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
  preserveOriginalId,
  hitSettings
}) => {
  const pyodide = await pyodideManager.init();
  pyodide.globals.set('GANFLU_INPUT_FASTA', fastaText);
  pyodide.globals.set('GANFLU_AUTO_GFF3_JSON', JSON.stringify(gff3ByTarget));
  pyodide.globals.set('GANFLU_ISOLATE', isolate);
  pyodide.globals.set('GANFLU_OUTPUT_STEM', outputStem);
  pyodide.globals.set('GANFLU_PRESERVE_ORIGINAL_ID', Boolean(preserveOriginalId));
  pyodide.globals.set('GANFLU_HIT_SETTINGS_JSON', JSON.stringify(hitSettings || {}));
  setStatus('Classifying contigs');
  const result = JSON.parse(
    pyodide.runPython(
      'run_ganflu_auto_web(GANFLU_INPUT_FASTA, GANFLU_AUTO_GFF3_JSON, GANFLU_ISOLATE, GANFLU_OUTPUT_STEM, GANFLU_PRESERVE_ORIGINAL_ID, GANFLU_HIT_SETTINGS_JSON)'
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
  setStatus('Preparing annotation');
  elements.logText.textContent = '';
  state.outputs = {};
  state.resultSummary = null;
  state.downloadSelection = DOWNLOAD_ALL_VALUE;
  state.selectedTargetTab = '';
  renderResults();

  try {
    const target = elements.target.value;

    const { text: fastaText, name } = await readInputFasta();
    if (!fastaText.trim()) throw new Error('FASTA input is empty.');

    const outputStem = makeSafeStem(elements.outputStem.value, makeSafeStem(name));
    const isolate = elements.isolate.value.trim() || outputStem;
    const hitSettings = getHitSettings();
    let result;
    if (target === 'auto') {
      const gff3ByTarget = {};
      for (const autoTarget of AUTO_TARGETS) {
        const reference = await getReferencePayload(autoTarget);
        gff3ByTarget[autoTarget] = await miniprotManager.run({
          genomeFasta: fastaText,
          proteinFasta: reference.protein_fasta,
          prefix: AUTO_PREFIXES[autoTarget] || 'MP',
          intronOpenPenalty: 15,
          bestN: hitSettings.maxSecondaryAlignments,
          outputScoreRatio: hitSettings.outputScoreRatio,
          secondaryToPrimaryRatio: hitSettings.secondaryToPrimaryRatio,
          statusLabel: `Running Miniprot (${autoTarget})`
        });
      }
      result = await runAutoGff3ToOutputs({
        fastaText,
        gff3ByTarget,
        isolate,
        outputStem,
        preserveOriginalId: elements.preserveOriginalId.checked,
        hitSettings
      });
      state.outputs = result.outputs || {};
    } else {
      const reference = await getReferencePayload(target);
      const gff3Text = await miniprotManager.run({
        genomeFasta: fastaText,
        proteinFasta: reference.protein_fasta,
        prefix: 'MP',
        intronOpenPenalty: 15,
        bestN: hitSettings.maxSecondaryAlignments,
        outputScoreRatio: hitSettings.outputScoreRatio,
        secondaryToPrimaryRatio: hitSettings.secondaryToPrimaryRatio
      });
      result = await runGff3ToOutputs({
        fastaText,
        gff3Text,
        target,
        isolate,
        outputStem,
        preserveOriginalId: elements.preserveOriginalId.checked,
        hitSettings
      });
      state.outputs = result.outputs || {};
    }
    state.resultSummary = result.summary || null;
    renderResults();
    if (!elements.logText.textContent) {
      elements.logText.textContent = result.log || 'Completed.';
    }
    setStatus('Completed');
  } catch (error) {
    const message = error?.message ? String(error.message) : String(error || 'Unknown error');
    elements.advancedRunDetails.hidden = false;
    elements.advancedRunBody.innerHTML = '';
    elements.logText.textContent = message;
    setStatus('Error');
  } finally {
    setRunning(false);
  }
});

elements.downloadButton.addEventListener('click', () => {
  try {
    downloadSelectedResult();
  } catch (error) {
    elements.advancedRunDetails.hidden = false;
    elements.logText.textContent = error?.message || String(error);
  }
});

elements.downloadSelect.addEventListener('change', () => {
  state.downloadSelection = elements.downloadSelect.value || DOWNLOAD_ALL_VALUE;
  renderDownloadOptions();
});

elements.runSummary.addEventListener('click', (event) => {
  const tab = event.target.closest('.target-tab');
  if (tab) {
    state.selectedTargetTab = tab.dataset.targetTab || '';
    renderRunSummary();
    return;
  }
  const button = event.target.closest('.sequence-copy');
  if (!button) return;
  copySequence(button).catch((error) => {
    elements.advancedRunDetails.hidden = false;
    elements.logText.textContent = error?.message || String(error);
  });
});

elements.clearButton.addEventListener('click', () => {
  elements.fastaFile.value = '';
  elements.fastaText.value = '';
  elements.target.value = 'IAV';
  elements.outputStem.value = '';
  elements.isolate.value = '';
  elements.preserveOriginalId.checked = false;
  resetHitSettings();
  state.loadedSampleKey = '';
  state.sampleInputDirty = false;
  renderSampleTabs();
  resetResultState();
  setStatus('Idle');
});

elements.sampleTabs.forEach((button) => {
  button.addEventListener('click', (event) => {
    event.preventDefault();
    loadSampleProfile(button.dataset.sampleKey, { announce: true }).catch(handleSampleLoadError);
  });
});

elements.fastaFile.addEventListener('change', () => {
  markSampleInputDirty();
  const file = elements.fastaFile.files?.[0];
  if (file) {
    elements.outputStem.value = makeSafeStem(file.name);
  }
});

elements.fastaText.addEventListener('input', markSampleInputDirty);
elements.target.addEventListener('change', markSampleInputDirty);
elements.outputStem.addEventListener('input', markSampleInputDirty);
elements.isolate.addEventListener('input', markSampleInputDirty);
elements.preserveOriginalId.addEventListener('change', markSampleInputDirty);
HIT_SETTING_DEFINITIONS.forEach(({ element }) => {
  element?.addEventListener('input', markSampleInputDirty);
});

elements.target.value = 'IAV';
setStatus('Idle');
renderResults();
renderSampleTabs();
