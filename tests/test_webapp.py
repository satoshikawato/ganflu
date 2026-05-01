import sys
import subprocess
import zipfile
from pathlib import Path

import ganflu
from ganflu.scripts import gff3togbk


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "ganflu" / "web"


def test_gff3togbk_accepts_raw_args_without_sys_argv(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["gff3togbk"])
    args = gff3togbk.parse_arguments(
        [
            "-i",
            "input.fa",
            "-g",
            "input.gff3",
            "--toml",
            "IAV.toml",
            "-o",
            "output.gbk",
            "--isolate",
            "A/Test/1/2026",
        ]
    )
    assert args.input == "input.fa"
    assert args.gff == "input.gff3"
    assert args.isolate == "A/Test/1/2026"


def test_webapp_assets_are_packaged_for_static_serving():
    required = [
        WEB_ROOT / "index.html",
        WEB_ROOT / "open-source-notices.html",
        WEB_ROOT / "js" / "app.js",
        WEB_ROOT / "js" / "app" / "pyodide.js",
        WEB_ROOT / "js" / "app" / "python-helpers.js",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "python_stdlib.zip",
        WEB_ROOT / "vendor" / "pyodide-wheels" / "numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl",
        WEB_ROOT / "vendor" / "pyodide-wheels" / "biopython-1.85-cp313-cp313-pyodide_2025_0_wasm32.whl",
        WEB_ROOT / "wasm" / "miniprot" / "dist" / "miniprot-ganflu.mjs",
        WEB_ROOT / "wasm" / "miniprot" / "dist" / "miniprot-ganflu.wasm",
    ]
    missing = [str(path.relative_to(REPO_ROOT)) for path in required if not path.exists()]
    assert not missing


def test_browser_wheel_contains_ganflu_db_without_web_assets():
    wheel_path = WEB_ROOT / f"ganflu-{ganflu.__version__}-py3-none-any.whl"
    assert wheel_path.exists()
    with zipfile.ZipFile(wheel_path) as zf:
        names = set(zf.namelist())
    assert "ganflu/db/IAV/IAV.toml" in names
    assert "ganflu/db/IDV/prot/IDV_proteome_consensus.faa" in names
    assert not any(name.startswith("ganflu/web/") for name in names)


def test_normal_wheel_contains_webapp_assets(tmp_path):
    dist_dir = tmp_path / "dist"
    subprocess.run(
        [sys.executable, "setup.py", "bdist_wheel", "--dist-dir", str(dist_dir)],
        cwd=REPO_ROOT,
        check=True,
    )
    wheel_path = dist_dir / f"ganflu-{ganflu.__version__}-py3-none-any.whl"
    assert wheel_path.exists()
    with zipfile.ZipFile(wheel_path) as zf:
        names = set(zf.namelist())
    assert "ganflu/web/index.html" in names
    assert f"ganflu/web/ganflu-{ganflu.__version__}-py3-none-any.whl" in names
    assert "ganflu/web/wasm/miniprot/dist/miniprot-ganflu.wasm" in names
