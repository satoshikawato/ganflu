[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/ganflu)
# ganflu
Genome ANnotator for inFLUenza viruses

## Usage
```bash
ganflu
usage: ganflu [-h] -i INPUT [-o OUTPUT] -t {IAV,IBV} [-d DB_DIR] [--isolate ISOLATE] [--preserve_original_id] [-v]

ganflu v0.1.0: Influenza virus genome annotation

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Input FASTA file
  -o, --output OUTPUT   basename for Output GenBank file name (default: <input>)
  -t, --target {IAV,IBV}
                        Target (IAV or IBV)
  -d, --db_dir DB_DIR   Data path (optional; default: ganflu/db)
  --isolate ISOLATE     isolate name (e.g. "A/Narita/1/2009", "A/goose/Guangdong/1/1996", "B/Lee/1940")
  --preserve_original_id, --preserve-original-id
                        Preserve original FASTA record IDs in GenBank output
  -v, --version         show program's version number and exit
```

Outputs are written using the output stem:

- `<output>.gbk`: GenBank annotation
- `<output>.gff3`: miniprot GFF3
- `<output>.cds.fna`: CDS nucleotide FASTA
- `<output>.faa`: amino acid FASTA
