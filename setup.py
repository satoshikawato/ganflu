import os
import re
from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent
BROWSER_WHEEL_BUILD_ENV = "GANFLU_BUILDING_BROWSER_WHEEL"

BASE_PACKAGE_DATA = [
    "db/*/*.toml",
    "db/*/prot/*.faa",
]

WEB_PACKAGE_DATA = [
    "web/*.html",
    "web/*.svg",
    "web/*.whl",
    "web/samples/*.fa",
    "web/samples/*.fasta",
    "web/samples/*.fna",
    "web/js/*.js",
    "web/js/app/*.js",
    "web/vendor/fonts/inter/*.woff2",
    "web/vendor/fonts/noto-sans-jp/*.woff2",
    "web/vendor/pyodide-wheels/*.whl",
    "web/vendor/pyodide/v0.29.0/full/*.js",
    "web/vendor/pyodide/v0.29.0/full/*.json",
    "web/vendor/pyodide/v0.29.0/full/*.mjs",
    "web/vendor/pyodide/v0.29.0/full/*.wasm",
    "web/vendor/pyodide/v0.29.0/full/*.whl",
    "web/vendor/pyodide/v0.29.0/full/*.zip",
    "web/wasm/miniprot/*.js",
    "web/wasm/miniprot/dist/*.mjs",
    "web/wasm/miniprot/dist/*.wasm",
]


def get_version():
    init_py = (ROOT / "ganflu" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__ = ["\']([^"\']+)["\']', init_py, re.MULTILINE)
    if not match:
        raise RuntimeError("Unable to find __version__ in ganflu/__init__.py")
    return match.group(1)


def get_package_data():
    if os.environ.get(BROWSER_WHEEL_BUILD_ENV) == "1":
        return BASE_PACKAGE_DATA
    return [*BASE_PACKAGE_DATA, *WEB_PACKAGE_DATA]


setup(
    name="ganflu",
    version=get_version(),
    packages=find_packages(),
    install_requires=[
        "biopython",
        "toml",
    ],
    include_package_data=False,
    package_data={"ganflu": get_package_data()},
    python_requires=">=3.10",
    author="Satoshi Kawato",
    author_email="kawato@kaiyodai.ac.jp",
    description="Genome ANnotation for inFLUenza viruses (GANFLU)",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url="https://github.com/satoshikawato/ganflu/",
    license="MIT",
    entry_points={"console_scripts": ["ganflu = ganflu.ganflu:main"]},
)
