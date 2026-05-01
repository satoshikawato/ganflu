#!/usr/bin/env python3
from __future__ import annotations

import re
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "ganflu" / "web"
INIT_PATH = REPO_ROOT / "ganflu" / "__init__.py"
CONFIG_PATH = WEB_ROOT / "js" / "config.js"
BROWSER_WHEEL_BUILD_ENV = "GANFLU_BUILDING_BROWSER_WHEEL"

VERSION_RE = re.compile(r'^__version__ = ["\']([^"\']+)["\']', re.MULTILINE)
WHEEL_NAME_RE = re.compile(r'^(export const GANFLU_WHEEL_NAME\s*=\s*")[^"]+(";\s*)$', re.MULTILINE)


def read_version() -> str:
    match = VERSION_RE.search(INIT_PATH.read_text(encoding="utf-8"))
    if match is None:
        raise RuntimeError(f"Could not find __version__ in {INIT_PATH}")
    return match.group(1)


def expected_wheel_name() -> str:
    return f"ganflu-{read_version()}-py3-none-any.whl"


def update_config(wheel_name: str) -> None:
    text = CONFIG_PATH.read_text(encoding="utf-8")
    updated, replacements = WHEEL_NAME_RE.subn(
        lambda match: f"{match.group(1)}{wheel_name}{match.group(2)}",
        text,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError(f"Could not update GANFLU_WHEEL_NAME in {CONFIG_PATH}")
    if updated != text:
        CONFIG_PATH.write_text(updated, encoding="utf-8")


def prepare_browser_wheel() -> Path:
    WEB_ROOT.mkdir(parents=True, exist_ok=True)
    wheel_name = expected_wheel_name()
    with tempfile.TemporaryDirectory(prefix="ganflu-browser-wheel-") as tmpdir:
        dist_dir = Path(tmpdir) / "dist"
        dist_dir.mkdir()
        env = os.environ.copy()
        env[BROWSER_WHEEL_BUILD_ENV] = "1"
        subprocess.run(
            [sys.executable, "setup.py", "bdist_wheel", "--dist-dir", str(dist_dir)],
            cwd=REPO_ROOT,
            env=env,
            check=True,
        )
        wheel_path = dist_dir / wheel_name
        if not wheel_path.exists():
            available = ", ".join(path.name for path in sorted(dist_dir.glob("*.whl")))
            raise FileNotFoundError(f"Expected {wheel_name}; built wheels: {available or 'none'}")
        for old_wheel in WEB_ROOT.glob("ganflu-*.whl"):
            old_wheel.unlink()
        target_path = WEB_ROOT / wheel_name
        shutil.copy2(wheel_path, target_path)
    update_config(wheel_name)
    return target_path


def main() -> int:
    target_path = prepare_browser_wheel()
    print(f"Prepared browser wheel: {target_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
