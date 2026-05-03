# Mixed Input Implementation Plan

## Goal

Support mixed IAV FASTA inputs such as `H1N1_H5N1_mix.fa` and
`H1N1_H7N9_mix.fna` without losing valid secondary hits or mixing HA/NA subtype
state across contigs.

The output should keep the existing model:

- one annotation run produces one isolate-level GenBank set
- isolate-level serotype is a single global value such as `HXN1` or `HXNX`
- HA/NA subtype details are stored only on the corresponding HA/NA CDS features

Expected examples:

```text
H1N1 + H5N1 -> global serotype HXN1
H1N1 + H7N9 -> global serotype HXNX
```

## Current Findings

### `auto` Is Not a Subtype Splitter

`ganflu -t auto` decides the influenza target class (`IAV`, `IBV`, `ICV`,
`IDV`) per contig. It does not split H1N1, H5N1, or H7N9 isolates. For mixed
IAV inputs, accepted contigs should still be annotated as `IAV`.

### Default `miniprot` Output Drops Valid Hits in Mixed Inputs

The current launcher calls:

```bash
miniprot -P MP --gff -J 15 input.fa IAV_proteome_consensus.faa
```

With mixed inputs, the same reference protein can have more than one valid
target contig. For example, `PB2` can match both the H1N1 PB2 and the H5N1/H7N9
PB2. With default secondary-output filtering, one valid hit may be omitted.

Observed local reproduction:

```text
H1N1+H5N1 default:  mRNA 44, CDS 46, contigs with hit 13/16
H1N1+H5N1 relaxed: mRNA 80, CDS 84, contigs with hit 16/16

H1N1+H7N9 default:  mRNA 41, CDS 43, contigs with hit 12/16
H1N1+H7N9 relaxed: mRNA 76, CDS 80, contigs with hit 16/16
```

The relaxed command used for testing was:

```bash
miniprot -P MP --gff -J 15 -N 100 --outs=0.10 -p 0.10 input.fa ref.faa
```

### Relaxed `miniprot` Must Be Followed by Pruning

Relaxed output restores valid secondary hits, but it also exposes extra
candidates:

- low-identity HA/NA subtype alternatives
- duplicate alignments with the same product and coordinates
- secondary hits that should not become extra GenBank features

Therefore the correct default is not simply "emit more and annotate all of it".
The pipeline should:

1. run `miniprot` with relaxed secondary-output settings
2. prune the raw GFF3 into one annotation-ready GFF3
3. pass only the pruned GFF3 to `gff3togbk`

### HA/NA Subtype State Is Currently Too Global

`gff3togbk` currently chooses HA/NA subtype candidates using shared state across
the full GFF3 processing pass. That allows one contig's HA/NA candidate to
affect another contig.

Subtype candidate selection must be contig-local:

- choose the best HA subtype for that contig's HA feature
- choose the best NA subtype for that contig's NA feature
- add `/note="subtype: H1"` only to HA
- add `/note="subtype: N1"` only to NA
- do not add subtype notes to PB2/PB1/PA/NP/M/NS

After features are built, the existing global serotype aggregation should collect
those HA/NA notes across all contigs.

## Design

### 1. Add Configurable `miniprot` Secondary Output

Update `ganflu/launchers/miniprot.py` so `MiniprotCommandLine` can accept:

```text
max_secondary_alignments -> -N
secondary_to_primary_ratio -> -p
output_score_ratio -> --outs
```

Use these defaults for ganflu annotation runs after pruning is available:

```text
max_secondary_alignments = 100
secondary_to_primary_ratio = 0.10
output_score_ratio = 0.10
```

The launcher should still allow stricter values later if needed, but the first
implementation can keep these as internal defaults rather than adding user-facing
CLI options.

### 2. Add a GFF3 Pruning Module

Add an importable helper, probably in `ganflu/scripts/gff3_prune.py`.

Responsibilities:

- parse `mRNA`, `CDS`, `stop_codon`, and related child features by parent ID
- preserve one selected mRNA parent and all of its child features
- remove exact duplicate alignments
- remove unselected HA/NA subtype alternatives
- write a valid GFF3 suitable for `gff3togbk`

Selection rules:

```text
For HA:
  per contig, keep only the best HA subtype candidate.

For NA:
  per contig, keep only the best NA subtype candidate.

For non-HA/NA products:
  per contig and product, keep the best parent alignment.

For fragmented products such as PA-X_fragment01/PA-X_fragment02:
  do not collapse different fragments into one candidate.
  Keep the best parent for each full product name.

For spliced or multi-CDS products such as M2/NS2:
  select the parent mRNA once, then keep all child CDS rows for that parent.
```

Ranking should mirror the current auto logic where possible:

```text
identity
reference amino-acid coverage
raw score
longer query span
lower Rank value if present
```

Exact ranking can be implemented as a single sort key and tested with focused
fixtures.

### 3. Use Pruned GFF3 in Fixed-Target Mode

For `ganflu -t IAV`, `IBV`, `ICV`, or `IDV`:

1. run relaxed `miniprot` to a raw temporary GFF3
2. prune raw GFF3 to the user-facing `<output>.gff3`
3. call `gff3togbk` on `<output>.gff3`

The public `.gff3` should be the annotation-ready pruned GFF3. Raw GFF3 should
not be an additional default output unless debugging later requires it.

### 4. Use Pruned GFF3 in Auto Annotation Output

For `ganflu -t auto`:

1. keep using the raw scan GFF3 for target/segment classification
2. after contig calls are accepted, prune/filter the raw scan GFF3 for each
   target
3. pass the pruned target GFF3 to `gff3togbk`

This replaces the current segment-only filtering behavior for annotation output.
The classifier can still see broad candidate evidence, but GenBank output should
receive only the selected annotation features.

### 5. Make HA/NA Subtype Selection Contig-Local

Update `ganflu/scripts/gff3togbk.py` so `to_seqfeatures()` processes each seqid
with a fresh antigen candidate dictionary.

The result should be:

```text
NC_002017.1 HA -> /note="subtype: H1"
NC_007362.1 HA -> /note="subtype: H5"
NC_002018.1 NA -> /note="subtype: N1"
NC_026429.1 NA -> /note="subtype: N9"
```

Then global aggregation remains isolate-level:

```text
H1 + H5, N1 only -> HXN1
H1 + H7, N1 + N9 -> HXNX
```

### 6. Keep Subtype Notes Limited to HA/NA

Do not add subtype notes to non-antigen genes. The source of truth for subtype
aggregation should be HA and NA features only.

## Test Plan

### Unit Tests

Add focused tests for the pruning helper:

- HA candidates on one contig keep only the highest-identity subtype
- NA candidates on one contig keep only the highest-identity subtype
- duplicate same-product/same-coordinate parents are deduplicated
- M1 and M2 are both retained because they are distinct products
- NS1 and NS2 are both retained because they are distinct products
- PA-X fragments are not accidentally collapsed into a single retained fragment

Add focused tests for `gff3togbk`:

- two HA contigs with different best subtypes keep separate HA subtype notes
- HA/NA notes aggregate to `HXN1` or `HXNX`
- PB2/PB1/PA/NP/M/NS records do not receive subtype notes

### CLI Smoke Tests

Add smoke tests guarded by `shutil.which("miniprot")`:

```text
ganflu -i tests/data/H1N1_H5N1_mix.fa -t IAV
ganflu -i tests/data/H1N1_H7N9_mix.fna -t IAV
ganflu -i tests/data/H1N1_H5N1_mix.fa -t auto
ganflu -i tests/data/H1N1_H7N9_mix.fna -t auto
```

Assertions:

- command exits successfully
- GenBank exists and has 16 records for fixed-target runs
- auto accepts 16 IAV contigs for these curated mixed inputs
- H1N1+H5N1 global organism/serotype contains `HXN1`
- H1N1+H7N9 global organism/serotype contains `HXNX`
- HA features have only HA subtype notes
- NA features have only NA subtype notes
- non-HA/NA features do not have subtype notes
- no repeated identical location is produced from duplicate parent alignments

### Regression Tests

Keep current smoke tests for:

- PR8 IAV
- IBV
- ICV
- IDV
- PR8 auto

The PR8 single-isolate output should remain complete and should not gain extra
features after relaxed `miniprot` plus pruning.

## Implementation Order

1. Add `miniprot` secondary-output parameters to the launcher.
2. Implement GFF3 parsing and pruning helper with unit tests.
3. Wire pruning into fixed-target CLI runs.
4. Wire pruning into auto target annotation output.
5. Make `gff3togbk` HA/NA subtype selection contig-local.
6. Add mixed-input CLI smoke tests.
7. Run the focused test set, then the full test suite.
8. Update README or auto-mode docs if the user-facing behavior needs mention.

## Non-Goals

- Do not split mixed inputs into separate isolate files.
- Do not make `auto` classify H1/H5/H7 as separate target classes.
- Do not annotate every relaxed secondary alignment.
- Do not add subtype notes to non-HA/NA genes.
- Do not add raw unpruned GFF3 as a default public output.

## Open Questions

No blocking questions for implementation.

Two details can be decided during implementation:

- whether to keep raw GFF3 as a debug-only temporary file or delete it after
  pruning
- whether relaxed `miniprot` parameters should remain internal defaults or later
  become advanced CLI options
