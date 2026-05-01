import { createGanfluMiniprot } from '../../wasm/miniprot/miniprot-ganflu.js';

export const createMiniprotManager = ({ onStatus = () => {} } = {}) => {
  let runtimePromise = null;

  const init = async () => {
    if (runtimePromise) return runtimePromise;
    runtimePromise = (async () => {
      onStatus('Loading Miniprot WebAssembly');
      const runtime = await createGanfluMiniprot();
      onStatus('Miniprot ready');
      return runtime;
    })();
    return runtimePromise;
  };

  const run = async ({ genomeFasta, proteinFasta, prefix = 'MP', intronOpenPenalty = 15 }) => {
    const runtime = await init();
    onStatus('Running Miniprot');
    return runtime.run({ genomeFasta, proteinFasta, prefix, intronOpenPenalty });
  };

  return {
    init,
    run
  };
};
