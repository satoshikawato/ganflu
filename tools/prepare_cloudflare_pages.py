#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
from html import escape
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "ganflu" / "web"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "dist" / "cloudflare-pages"
ANALYTICS_TOKEN_ENV = "CLOUDFLARE_WEB_ANALYTICS_TOKEN"
SCRIPT_MARKER = "<!-- CLOUDFLARE_WEB_ANALYTICS_SCRIPT -->"
NOTICE_MARKER = "<!-- CLOUDFLARE_WEB_ANALYTICS_NOTICE -->"
SCRIPT_SRC_BASE = "script-src 'self' 'unsafe-eval' 'wasm-unsafe-eval';"
SCRIPT_SRC_ANALYTICS = (
    "script-src 'self' 'unsafe-eval' 'wasm-unsafe-eval' "
    "https://static.cloudflareinsights.com;"
)
CONNECT_SRC_BASE = "connect-src 'self';"
CONNECT_SRC_ANALYTICS = "connect-src 'self' https://cloudflareinsights.com;"


def _load_prepare_browser_wheel_module():
    module_path = REPO_ROOT / "tools" / "prepare_browser_wheel.py"
    spec = spec_from_file_location("prepare_browser_wheel", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load browser wheel preparation module from {module_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _replace_once(source: str, old: str, new: str) -> str:
    if old not in source:
        raise RuntimeError(f"Expected marker not found while preparing Cloudflare bundle: {old}")
    return source.replace(old, new, 1)


def _normalize_analytics_token(token: str | None) -> str | None:
    if token is None:
        return None
    token = token.strip()
    return token or None


def _render_analytics_script(token: str) -> str:
    beacon_config = escape(json.dumps({"token": token}, separators=(",", ":")), quote=True)
    return (
        "    <!-- Cloudflare Web Analytics -->\n"
        "    <script defer src=\"https://static.cloudflareinsights.com/beacon.min.js\" "
        f"data-cf-beacon='{beacon_config}'></script>\n"
        "    <!-- End Cloudflare Web Analytics -->"
    )


def _render_analytics_notice() -> str:
    return (
        "        <span class=\"hosted-notice\">"
        "Hosted deployment notice: Cloudflare Web Analytics records aggregate page-visit "
        "and performance metrics. Uploaded FASTA files are processed locally in your browser "
        "and are not sent to Cloudflare by ganflu."
        "</span>"
    )


def build_cloudflare_pages_bundle(
    *,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    analytics_token: str | None = None,
) -> Path:
    output_root = Path(output_root)
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(WEB_ROOT, output_root)

    token = _normalize_analytics_token(analytics_token)
    index_path = output_root / "index.html"
    index_html = index_path.read_text(encoding="utf-8")
    if token:
        index_html = _replace_once(index_html, SCRIPT_SRC_BASE, SCRIPT_SRC_ANALYTICS)
        index_html = _replace_once(index_html, CONNECT_SRC_BASE, CONNECT_SRC_ANALYTICS)
    index_html = _replace_once(
        index_html,
        SCRIPT_MARKER,
        _render_analytics_script(token) if token else "",
    )
    index_html = _replace_once(
        index_html,
        NOTICE_MARKER,
        _render_analytics_notice() if token else "",
    )
    index_path.write_text(index_html, encoding="utf-8")
    return output_root


def prepare_cloudflare_pages(
    *,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    analytics_token: str | None = None,
    analytics_enabled: bool = True,
) -> Path:
    prepare_browser_wheel_module = _load_prepare_browser_wheel_module()
    prepare_browser_wheel_module.prepare_browser_wheel()

    token = None
    if analytics_enabled:
        token = analytics_token if analytics_token is not None else os.environ.get(ANALYTICS_TOKEN_ENV)
    return build_cloudflare_pages_bundle(output_root=output_root, analytics_token=token)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the static Cloudflare Pages bundle in dist/cloudflare-pages with the "
            "current ganflu browser wheel."
        )
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Output directory for the deployable static bundle.",
    )
    analytics_group = parser.add_mutually_exclusive_group()
    analytics_group.add_argument(
        "--analytics-token",
        help=(
            "Cloudflare Web Analytics token. Defaults to "
            f"{ANALYTICS_TOKEN_ENV} when set."
        ),
    )
    analytics_group.add_argument(
        "--no-analytics",
        action="store_true",
        help="Do not inject Cloudflare Web Analytics or the hosted analytics notice.",
    )
    args = parser.parse_args(argv)

    output_root = prepare_cloudflare_pages(
        output_root=args.output_root,
        analytics_token=args.analytics_token,
        analytics_enabled=not args.no_analytics,
    )
    print(f"Prepared Cloudflare Pages bundle: {output_root.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
