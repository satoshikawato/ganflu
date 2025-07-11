# ganflu
Genome ANnotator for inFLUenza viruses

## Usage
```bash
ganflu
usage: ganflu [-h] -i INPUT [-o OUTPUT] -t {IAV,IBV} [-d DB_DIR] [--isolate ISOLATE] [-v]

ganflu v0.1.0: Influenza virus genome annotation

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input FASTA file
  -o, --output OUTPUT   basename for Output GenBank file name (default: <input>)
  -t, --target {IAV,IBV}
                        Target (IAV or IBV)
  -d, --db_dir DB_DIR   Data path (optional; default: ganflu/db)
  --isolate ISOLATE     isolate name (e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940")
  -v, --version         show program's version number and exit
```
