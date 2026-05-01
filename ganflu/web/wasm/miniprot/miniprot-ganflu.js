import createMiniprotModule from './dist/miniprot-ganflu.mjs';

const encoder = new TextEncoder();
const decoder = new TextDecoder();

const toBytes = (value, label) => {
  if (typeof value === 'string') return encoder.encode(value);
  if (value instanceof Uint8Array) return value;
  if (value instanceof ArrayBuffer) return new Uint8Array(value);
  throw new TypeError(`${label} must be a string, Uint8Array, or ArrayBuffer.`);
};

const copyToWasm = (module, bytes) => {
  const ptr = module._miniprot_web_alloc(bytes.length);
  if (!ptr) throw new Error('miniprot wasm allocation failed.');
  module.HEAPU8.set(bytes, ptr);
  return { ptr, len: bytes.length };
};

const readWasmText = (module, ptr, len) => {
  if (!ptr || !len) return '';
  const bytes = new Uint8Array(module.HEAPU8.buffer, ptr, len).slice();
  return decoder.decode(bytes);
};

export const createGanfluMiniprot = async (moduleOptions = {}) => {
  const module = await createMiniprotModule(moduleOptions);

  const run = ({
    genomeFasta,
    proteinFasta,
    prefix = 'MP',
    intronOpenPenalty = 15
  } = {}) => {
    const allocations = [
      copyToWasm(module, toBytes(genomeFasta, 'genomeFasta')),
      copyToWasm(module, toBytes(proteinFasta, 'proteinFasta')),
      copyToWasm(module, encoder.encode(`${prefix || 'MP'}\0`))
    ];

    try {
      const exitCode = module._miniprot_ganflu_run(
        allocations[0].ptr,
        allocations[0].len,
        allocations[1].ptr,
        allocations[1].len,
        allocations[2].ptr,
        Number.isFinite(intronOpenPenalty) ? Math.trunc(intronOpenPenalty) : 15
      );
      const ptr = exitCode === 0
        ? module._miniprot_ganflu_result_ptr()
        : module._miniprot_ganflu_error_ptr();
      const len = exitCode === 0
        ? module._miniprot_ganflu_result_len()
        : module._miniprot_ganflu_error_len();
      const text = readWasmText(module, ptr, len);
      if (exitCode !== 0) {
        throw new Error(text || `miniprot exited with code ${exitCode}`);
      }
      return text;
    } finally {
      allocations.forEach(({ ptr, len }) => module._miniprot_web_dealloc(ptr, len));
      module._miniprot_ganflu_clear();
    }
  };

  return { module, run };
};

export const runGanfluMiniprot = async (options, moduleOptions = {}) => {
  const runtime = await createGanfluMiniprot(moduleOptions);
  return runtime.run(options);
};
