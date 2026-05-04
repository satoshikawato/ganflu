import json
import re
import sys
import subprocess
import zipfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import ganflu
from ganflu.scripts import gff3togbk


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "ganflu" / "web"


def load_python_helpers_namespace():
    helper_text = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    match = re.search(r"export const PYTHON_HELPERS = `\n(.*)\n`;\n?\Z", helper_text, re.S)
    assert match is not None
    namespace = {}
    exec(match.group(1), namespace)
    return namespace


def load_cloudflare_pages_module():
    module_path = REPO_ROOT / "tools" / "prepare_cloudflare_pages.py"
    spec = spec_from_file_location("prepare_cloudflare_pages", module_path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
        WEB_ROOT / "_headers",
        WEB_ROOT / "open-source-notices.html",
        WEB_ROOT / "js" / "app.js",
        WEB_ROOT / "js" / "app" / "pyodide.js",
        WEB_ROOT / "js" / "app" / "python-helpers.js",
        WEB_ROOT / "samples" / "IAV_PR8.fasta",
        WEB_ROOT / "samples" / "IBV_B_Victoria_2_1987.fa",
        WEB_ROOT / "samples" / "ICV_Ann_Arbor_1_1950.fna",
        WEB_ROOT / "samples" / "IDV_swine_Oklahoma_1334_2011.fna",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "python_stdlib.zip",
        WEB_ROOT / "vendor" / "pyodide-wheels" / "numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl",
        WEB_ROOT / "vendor" / "pyodide-wheels" / "biopython-1.85-cp313-cp313-pyodide_2025_0_wasm32.whl",
        WEB_ROOT / "wasm" / "miniprot" / "miniprot-ganflu.js",
    ]
    missing = [str(path.relative_to(REPO_ROOT)) for path in required if not path.exists()]
    assert not missing
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert index_html.index('<option value="IAV">IAV</option>') < index_html.index(
        '<option value="auto">AUTO</option>'
    )
    assert index_html.index('data-sample-key="IDV"') < index_html.index('data-sample-key="auto"')
    assert 'placeholder="defaults to output stem"' in index_html
    assert 'id="run-summary"' in index_html
    assert "Download Results" in index_html
    assert '<script type="module" src="./js/app.js?v=20260504"></script>' in index_html
    assert 'id="download-select"' in index_html
    assert 'id="advanced-run-details"' in index_html
    assert 'id="output-tabs"' not in index_html
    assert 'id="summary"' not in index_html
    assert 'id="output-text"' not in index_html
    assert "target-tabs" in index_html
    assert "target-tab-panel" in index_html
    assert 'id="advanced-hit-filters"' in index_html
    assert "sequence-boxes" in index_html
    assert 'class="help-tip"' in index_html
    assert "Output stem is updated from the file name" in index_html
    assert "Load sample data" in index_html
    assert 'data-sample-key="auto" title="Load mixed IAV, IBV, ICV, and IDV sample data.">AUTO</button>' in index_html
    assert 'data-sample-key="IAV"' in index_html
    assert 'data-sample-key="IDV"' in index_html
    assert "Minimum normalized hit score" in index_html
    assert "Miniprot output filter for secondary hits" in index_html
    assert 'id="min-identity"' in index_html
    assert 'id="min-identity" type="number" min="0" max="1" step="0.01" value="0.70"' in index_html
    assert 'id="max-secondary-alignments"' in index_html
    app_js = (WEB_ROOT / "js" / "app.js").read_text(encoding="utf-8")
    headers_text = (WEB_ROOT / "_headers").read_text(encoding="utf-8")
    miniprot_js = (WEB_ROOT / "js" / "app" / "miniprot.js").read_text(encoding="utf-8")
    helpers_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    assert "runAutoGff3ToOutputs" in app_js
    assert "renderRunSummary" in app_js
    assert "createZipBlob" in app_js
    assert "downloadResultsZip" in app_js
    assert "downloadSelectedResult" in app_js
    assert "formatDownloadLabel" in app_js
    assert "getAutoContigTabKey" in app_js
    assert "groupAutoContigsByTab" in app_js
    assert "Unclassified" in app_js
    assert "formatSequenceFasta" in app_js
    assert "Copied!" in app_js
    assert "button.disabled = true" not in app_js
    assert "copySequence" in app_js
    assert "SAMPLE_PROFILES" in app_js
    assert "loadSampleProfile('auto')" not in app_js
    assert "cache: 'no-cache'" in app_js
    assert "Loading sample data" in app_js
    assert "Sample load failed" in app_js
    assert "elements.target.value = 'IAV'" in app_js
    assert "elements.target.value = 'auto'" not in app_js
    assert "elements.sampleTabs" in app_js
    assert "renderSampleTabs" in app_js
    assert "A/Puerto Rico/8/1934" in app_js
    assert "D/swine/Oklahoma/1334/2011" in app_js
    assert "./samples/IAV_PR8.fasta" in app_js
    assert "renderTabs" not in app_js
    assert "GANFLU_HIT_SETTINGS_JSON" in app_js
    assert "maxSecondaryAlignments" in app_js
    assert "elements.outputStem.value = makeSafeStem(file.name)" in app_js
    assert "bestN = 100" in miniprot_js
    assert "outputScoreRatio = 0.1" in miniprot_js
    assert "secondaryToPrimaryRatio = 0.1" in miniprot_js
    assert "/js/*" in headers_text
    assert "/samples/*" in headers_text
    assert "Cache-Control: public, max-age=0, must-revalidate" in headers_text
    assert "run_ganflu_auto_web" in helpers_js
    assert "_normalize_hit_settings" in helpers_js
    assert "_build_genbank_feature_data" in helpers_js
    assert 'f"{stem}.summary.json"' in helpers_js


def test_local_index_keeps_cloudflare_analytics_deploy_only():
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "https://static.cloudflareinsights.com/beacon.min.js" not in index_html
    assert "static.cloudflareinsights.com" not in index_html
    assert "cloudflareinsights.com" not in index_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_SCRIPT -->" in index_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_NOTICE -->" in index_html


def test_cloudflare_bundle_copies_required_static_assets(tmp_path):
    module = load_cloudflare_pages_module()
    output_root = module.build_cloudflare_pages_bundle(output_root=tmp_path / "cloudflare-pages")

    required = [
        output_root / "index.html",
        output_root / "_headers",
        output_root / "open-source-notices.html",
        output_root / "js" / "app.js",
        output_root / "js" / "config.js",
        output_root / "samples" / "IAV_PR8.fasta",
        output_root / "samples" / "IBV_B_Victoria_2_1987.fa",
        output_root / "samples" / "ICV_Ann_Arbor_1_1950.fna",
        output_root / "samples" / "IDV_swine_Oklahoma_1334_2011.fna",
        output_root / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        output_root / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        output_root / "vendor" / "pyodide" / "v0.29.0" / "full" / "python_stdlib.zip",
        output_root / "wasm" / "miniprot" / "miniprot-ganflu.js",
        output_root / "wasm" / "miniprot" / "dist" / "miniprot-ganflu.wasm",
    ]
    missing = [str(path.relative_to(output_root)) for path in required if not path.exists()]
    assert not missing
    assert list((output_root / "vendor" / "pyodide-wheels").glob("biopython-*.whl"))
    assert "Cache-Control: public, max-age=0, must-revalidate" in (
        output_root / "_headers"
    ).read_text(encoding="utf-8")

    config_text = (output_root / "js" / "config.js").read_text(encoding="utf-8")
    match = re.search(r'GANFLU_WHEEL_NAME\s*=\s*"([^"]+)"', config_text)
    assert match is not None
    assert sorted(path.name for path in output_root.glob("ganflu-*.whl")) == [match.group(1)]


def test_cloudflare_bundle_includes_analytics_only_when_requested(tmp_path):
    module = load_cloudflare_pages_module()

    no_analytics_root = module.build_cloudflare_pages_bundle(
        output_root=tmp_path / "cloudflare-pages-no-analytics"
    )
    no_analytics_html = (no_analytics_root / "index.html").read_text(encoding="utf-8")
    assert "static.cloudflareinsights.com" not in no_analytics_html
    assert "cloudflareinsights.com" not in no_analytics_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_SCRIPT -->" not in no_analytics_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_NOTICE -->" not in no_analytics_html

    analytics_root = module.build_cloudflare_pages_bundle(
        output_root=tmp_path / "cloudflare-pages-analytics",
        analytics_token="test-token",
    )
    analytics_html = (analytics_root / "index.html").read_text(encoding="utf-8")
    assert "https://static.cloudflareinsights.com/beacon.min.js" in analytics_html
    assert "test-token" in analytics_html
    assert "Uploaded FASTA files are processed locally in your browser" in analytics_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_SCRIPT -->" not in analytics_html
    assert "<!-- CLOUDFLARE_WEB_ANALYTICS_NOTICE -->" not in analytics_html

    source_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "test-token" not in source_html
    assert "https://static.cloudflareinsights.com/beacon.min.js" not in source_html


def test_wrangler_points_at_cloudflare_pages_bundle():
    wrangler_toml = (REPO_ROOT / "wrangler.toml").read_text(encoding="utf-8")
    assert 'name = "ganflu"' in wrangler_toml
    assert 'pages_build_output_dir = "./dist/cloudflare-pages"' in wrangler_toml


def test_single_target_web_helper_skips_no_hit_contigs():
    helpers = load_python_helpers_namespace()
    fasta = ">hit\nATGAAATAA\n>nohit\nATGAAATAA\n"
    gff3 = "\n".join(
        [
            "##gff-version 3",
            "hit\tminiprot\tmRNA\t1\t9\t1\t+\t.\tID=MP000001;Rank=1;Identity=0.9500;Positive=0.9500;Target=PB2 1 3",
            "hit\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP000001;Identity=0.9500;Target=PB2 1 3",
        ]
    ) + "\n"

    hit_settings = json.dumps(
        {
            "minIdentity": 0.9,
            "minAaCoverage": 0.001,
            "minScore": 0.001,
            "minMargin": 0.2,
            "completeAaCoverage": 0.001,
            "maxSecondaryAlignments": 7,
            "outputScoreRatio": 0.2,
            "secondaryToPrimaryRatio": 0.3,
        }
    )

    result = json.loads(
        helpers["run_ganflu_web"](fasta, gff3, "IAV", "sample", "sample", False, hit_settings)
    )

    assert "error" not in result
    assert sorted(result["outputs"]) == [
        "sample.cds.fna",
        "sample.faa",
        "sample.gbk",
        "sample.gff3",
        "sample.hits.tsv",
        "sample.log",
        "sample.summary.json",
    ]
    assert result["summary"]["input_contigs"] == 2
    assert result["summary"]["annotated_contigs"] == 1
    assert result["summary"]["skipped_contigs"] == 1
    assert result["summary"]["passed_contigs"] == 1
    assert result["summary"]["failed_contigs"] == 1
    assert result["summary"]["by_qc_result"] == {"fail": 1, "pass": 1}
    assert result["summary"]["by_segment"] == {"PB2": 1}
    assert result["summary"]["thresholds"]["min_identity"] == 0.9
    assert result["summary"]["thresholds"]["min_aa_coverage"] == 0.001
    assert result["summary"]["thresholds"]["complete_aa_coverage"] == 0.001
    assert result["summary"]["miniprot"]["max_secondary_alignments"] == 7
    assert result["summary"]["miniprot"]["output_score_ratio"] == 0.2
    assert result["summary"]["schema_version"] == 1
    assert result["summary"]["mode"] == "target"
    assert result["summary"]["run"]["counts"]["input_contigs"] == 2
    assert result["summary"]["run"]["counts"]["annotated_contigs"] == 1
    assert result["summary"]["downloads"][-1]["kind"] == "log"
    assert result["summary"]["run_log"]["thresholds"]["min_identity"] == 0.9

    summary_json = json.loads(result["outputs"]["sample.summary.json"])
    assert summary_json["schema_version"] == 1
    assert len(summary_json["contigs"]) == 2
    hit_contig = next(contig for contig in summary_json["contigs"] if contig["input_id"] == "hit")
    assert hit_contig["record_id"] == "sample_PB2"
    assert hit_contig["target"] == "IAV"
    assert hit_contig["segment"] == "PB2"
    assert hit_contig["qc_result"] == "pass"
    feature = hit_contig["features"][0]
    assert feature["gene"] == "PB2"
    assert feature["product"] == "polymerase PB2"
    assert feature["metrics"]["identity"] == 0.95
    assert feature["metrics"]["aa_coverage"] == 1.0
    assert feature["reference"]["product"] == "PB2"
    assert feature["sequences"]["cds_nt"] == "ATGAAATAA"
    assert feature["sequences"]["aa"] == "MK"
    assert feature["flags"] == []
    assert "sample.log" in result["outputs"]

    strict_result = json.loads(
        helpers["run_ganflu_web"](
            fasta,
            gff3,
            "IAV",
            "sample",
            "sample",
            False,
            json.dumps({"minIdentity": 0.99, "minAaCoverage": 0.001, "minScore": 0.001}),
        )
    )
    assert strict_result["error"]["type"] == "ValueError"
    assert "No segment-compatible miniprot hits" in strict_result["error"]["message"]


def test_web_summary_displays_missing_stop_on_feature_not_contig():
    helpers = load_python_helpers_namespace()
    fasta = ">hit\nATGAAAAAA\n"
    gff3 = "\n".join(
        [
            "##gff-version 3",
            "hit\tminiprot\tmRNA\t1\t9\t1\t+\t.\tID=MP000001;Rank=1;Identity=0.9500;Positive=0.9500;Target=PB2 1 3",
            "hit\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP000001;Identity=0.9500;Target=PB2 1 3",
        ]
    ) + "\n"
    hit_settings = json.dumps(
        {"minIdentity": 0.9, "minAaCoverage": 0.001, "minScore": 0.001, "completeAaCoverage": 0.001}
    )

    result = json.loads(
        helpers["run_ganflu_web"](fasta, gff3, "IAV", "sample", "sample", False, hit_settings)
    )

    summary = json.loads(result["outputs"]["sample.summary.json"])
    hit_contig = summary["contigs"][0]
    feature = hit_contig["features"][0]
    assert "missing_stop" not in hit_contig["flags"]
    assert feature["flags"] == ["missing_stop"]


def test_web_summary_does_not_mark_nonterminal_feature_fragment_missing_stop():
    helpers = load_python_helpers_namespace()
    fasta = ">hit\nATGAAA\n"
    gff3 = "\n".join(
        [
            "##gff-version 3",
            "##PAF\tPB2\t3\t0\t2\t+\thit\t6\t0\t6\t6\t6\t0\tAS:i:1",
            "hit\tminiprot\tmRNA\t1\t6\t1\t+\t.\tID=MP000001;Rank=1;Identity=0.9500;Positive=0.9500;Target=PB2 1 2",
            "hit\tminiprot\tCDS\t1\t6\t.\t+\t0\tParent=MP000001;Identity=0.9500;Target=PB2 1 2",
        ]
    ) + "\n"
    hit_settings = json.dumps(
        {"minIdentity": 0.9, "minAaCoverage": 0.001, "minScore": 0.001, "completeAaCoverage": 0.001}
    )

    result = json.loads(
        helpers["run_ganflu_web"](fasta, gff3, "IAV", "sample", "sample", False, hit_settings)
    )

    summary = json.loads(result["outputs"]["sample.summary.json"])
    hit_contig = summary["contigs"][0]
    feature = hit_contig["features"][0]
    assert "missing_stop" not in hit_contig["flags"]
    assert feature["flags"] == []
    assert not any("stop codon not found" in note for note in feature["notes"])


def test_auto_web_helper_summary_downloads_and_features():
    helpers = load_python_helpers_namespace()
    fasta = ">hit\nATGAAATAA\n>nohit\nATGAAATAA\n"
    gff3 = "\n".join(
        [
            "##gff-version 3",
            "hit\tminiprot\tmRNA\t1\t9\t1\t+\t.\tID=MP000001;Rank=1;Identity=0.9500;Positive=0.9500;Target=PB2 1 3",
            "hit\tminiprot\tCDS\t1\t9\t.\t+\t0\tParent=MP000001;Identity=0.9500;Target=PB2 1 3",
        ]
    ) + "\n"
    hit_settings = json.dumps(
        {
            "minIdentity": 0.9,
            "minAaCoverage": 0.001,
            "minScore": 0.001,
            "minMargin": 0.2,
            "completeAaCoverage": 0.001,
        }
    )

    result = json.loads(
        helpers["run_ganflu_auto_web"](
            fasta,
            json.dumps({"IAV": gff3}),
            "sample",
            "sample",
            False,
            hit_settings,
        )
    )

    assert "error" not in result
    assert sorted(result["outputs"]) == [
        "sample.IAV.cds.fna",
        "sample.IAV.faa",
        "sample.IAV.gbk",
        "sample.IAV.gff3",
        "sample.auto.log",
        "sample.auto.summary.json",
        "sample.auto.tsv",
    ]
    summary = json.loads(result["outputs"]["sample.auto.summary.json"])
    assert summary["schema_version"] == 1
    assert summary["mode"] == "auto"
    assert summary["run"]["target"] == "auto"
    assert summary["run"]["targets_scanned"] == ["IAV"]
    assert summary["run"]["counts"]["annotated_contigs"] == 1
    assert {download["kind"] for download in summary["downloads"]} >= {
        "summary_json",
        "gff3",
        "genbank",
        "cds_fasta",
        "protein_fasta",
        "contig_report",
        "log",
    }
    assert any(download["name"] == "sample.IAV.gff3" and download["target"] == "IAV" for download in summary["downloads"])
    hit_contig = next(contig for contig in summary["contigs"] if contig["input_id"] == "hit")
    assert hit_contig["features"][0]["sequences"]["cds_nt"] == "ATGAAATAA"
    assert hit_contig["features"][0]["sequences"]["aa"] == "MK"


def test_browser_wheel_contains_ganflu_db_without_web_assets():
    wheel_path = WEB_ROOT / f"ganflu-{ganflu.__version__}-py3-none-any.whl"
    assert wheel_path.exists()
    with zipfile.ZipFile(wheel_path) as zf:
        names = set(zf.namelist())
    assert "ganflu/db/IAV/IAV.toml" in names
    assert "ganflu/scripts/auto_mode.py" in names
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
    assert "ganflu/web/wasm/miniprot/miniprot-ganflu.js" in names
    assert "ganflu/web/samples/IAV_PR8.fasta" in names
    assert "ganflu/web/samples/IDV_swine_Oklahoma_1334_2011.fna" in names
