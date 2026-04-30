import re
from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent


def get_version():
    init_py = (ROOT / "ganflu" / "__init__.py").read_text(encoding="utf-8")
    match = re.search(r'^__version__ = ["\']([^"\']+)["\']', init_py, re.MULTILINE)
    if not match:
        raise RuntimeError("Unable to find __version__ in ganflu/__init__.py")
    return match.group(1)


setup(
    name="ganflu",
    version=get_version(),
    packages=find_packages(),
    install_requires=[
        "biopython",
        "toml",
    ],
    include_package_data=True,
    package_data={
        "ganflu": [
            "db/*/*.toml",
            "db/*/prot/*.faa",
        ],
    },
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
