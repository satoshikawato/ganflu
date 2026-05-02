[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/ganflu)
# ganflu
Genome ANnotator for inFLUenza viruses

## Usage
```bash
ganflu
usage: ganflu [-h] -i INPUT [-o OUTPUT] -t {IAV,IBV,ICV,IDV,auto} [-d DB_DIR]
              [--isolate ISOLATE] [--preserve_original_id] [--log-file LOG_FILE]
              [--verbose] [-v]

ganflu v0.1.0: Influenza virus genome annotation

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input FASTA file
  -o, --output OUTPUT   basename for Output GenBank file name (default: <input>)
  -t, --target {IAV,IBV,ICV,IDV,auto}
                        Target virus
  -d, --db_dir DB_DIR   Data path (optional; default: ganflu/db)
  --isolate ISOLATE     isolate name (default: output stem prefix; e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940", "C/Ann_Arbor/1/1950", "D/swine/Oklahoma/1334/2011")
  --preserve_original_id, --preserve-original-id
                        Preserve original FASTA record IDs in GenBank output
  --log-file LOG_FILE   Log file path (default: <output>.log; auto mode: <output>.auto.log)
  --verbose             Show debug logs in the terminal
  -v, --version         show program's version number and exit

Subcommands:
  ganflu gui              Launch the browser-based ganflu Web app
  ganflu gui --help       Show ganflu Web app options
```

Outputs are written using the output stem:

- `<output>.gbk`: GenBank annotation
- `<output>.gff3`: miniprot GFF3
- `<output>.cds.fna`: CDS nucleotide FASTA
- `<output>.faa`: amino acid FASTA

## Auto mode

Use `-t auto` to scan each input contig against the packaged IAV, IBV, ICV,
and IDV references, classify the best target/segment, and annotate accepted
contigs by target.

```bash
ganflu -i contigs.fa -o sample -t auto
```

Auto mode writes target-suffixed annotation outputs plus reports:

- `<output>.<target>.gbk`
- `<output>.<target>.gff3`
- `<output>.<target>.cds.fna`
- `<output>.<target>.faa`
- `<output>.auto.tsv`
- `<output>.auto.summary.json`

## Web app

The static browser app is in `ganflu/web/` and runs Miniprot WebAssembly plus
Pyodide without a backend server.

Launch it from the CLI:

```bash
ganflu gui
```

Then open the URL printed in the terminal. By default ganflu asks the OS for a
free local port. Use `--open-browser` to ask ganflu to open the URL with the
system browser.

```bash
ganflu gui --port 8888 --open-browser
```

For source checkouts, the same static app can also be served manually:

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
