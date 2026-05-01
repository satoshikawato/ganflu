import {
  GANFLU_WHEEL_NAME,
  PYODIDE_INDEX_URL,
  PYODIDE_LOCAL_WHEELS
} from '../config.js';
import { PYTHON_HELPERS } from './python-helpers.js';

const resolveAssetUrl = (path) => new URL(path, window.location.href).toString();

const ensureAsset = async (path, label) => {
  const url = resolveAssetUrl(path);
  const response = await fetch(url, { method: 'HEAD', cache: 'no-store' });
  if (!response.ok) {
    throw new Error(`Missing packaged asset: ${label} (${response.status})`);
  }
  return url;
};

export const createPyodideManager = ({ onStatus = () => {} } = {}) => {
  let pyodidePromise = null;
  let pyodide = null;

  const init = async () => {
    if (pyodidePromise) return pyodidePromise;
    pyodidePromise = (async () => {
      onStatus('Loading Python runtime');
      const pyodideIndexUrl = resolveAssetUrl(PYODIDE_INDEX_URL);
      pyodide = await loadPyodide({
        indexURL: pyodideIndexUrl,
        packageBaseUrl: pyodideIndexUrl
      });

      onStatus('Loading Python installer');
      await pyodide.loadPackage('micropip');
      const micropip = pyodide.pyimport('micropip');

      onStatus('Installing Python dependencies');
      const dependencyUrls = await Promise.all(
        PYODIDE_LOCAL_WHEELS.map((path) => ensureAsset(path, path))
      );
      await micropip.install(dependencyUrls);

      onStatus('Installing ganflu');
      const wheelUrl = await ensureAsset(GANFLU_WHEEL_NAME, GANFLU_WHEEL_NAME);
      await micropip.install(wheelUrl);
      await pyodide.runPythonAsync(PYTHON_HELPERS);
      onStatus('Python ready');
      return pyodide;
    })();
    return pyodidePromise;
  };

  const getPyodide = () => pyodide;

  return {
    init,
    getPyodide
  };
};
