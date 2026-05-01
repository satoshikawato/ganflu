[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/ganflu)
# ganflu
Genome ANnotator for inFLUenza viruses

## Usage
```bash
ganflu
usage: ganflu [-h] -i INPUT [-o OUTPUT] -t {IAV,IBV,ICV,IDV} [-d DB_DIR] --isolate
              ISOLATE [--preserve_original_id] [--log-file LOG_FILE]
              [--verbose] [-v]

ganflu v0.1.0: Influenza virus genome annotation

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input FASTA file
  -o, --output OUTPUT   basename for Output GenBank file name (default: <input>)
  -t, --target {IAV,IBV,ICV,IDV}
                        Target virus
  -d, --db_dir DB_DIR   Data path (optional; default: ganflu/db)
  --isolate ISOLATE     isolate name (e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940", "C/Ann_Arbor/1/1950", "D/swine/Oklahoma/1334/2011")
  --preserve_original_id, --preserve-original-id
                        Preserve original FASTA record IDs in GenBank output
  --log-file LOG_FILE   Log file path (default: <output>.log)
  --verbose             Show debug logs in the terminal
  -v, --version         show program's version number and exit
```

Outputs are written using the output stem:

- `<output>.gbk`: GenBank annotation
- `<output>.gff3`: miniprot GFF3
- `<output>.cds.fna`: CDS nucleotide FASTA
- `<output>.faa`: amino acid FASTA

## Web app

The static browser app is in `ganflu/web/` and runs Miniprot WebAssembly plus
Pyodide without a backend server.

```bash
python tools/prepare_browser_wheel.py
cd ganflu/web
python -m http.server 8765 --bind 127.0.0.1
```

Then open <http://127.0.0.1:8765/>.

For a distributable package that includes the web app assets, build in this
order:

```bash
python tools/prepare_browser_wheel.py
python -m build
```
