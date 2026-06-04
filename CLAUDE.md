# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

The Kelly Lab's **Metabarcoding Taxonomic Assignment Pipeline** — takes de-multiplexed FASTQ files and produces taxon-by-sample count matrices. Written in R with external calls to shell tools and a remote BLAST cluster.

## How to Run

This is an interactive R script, not a package or CLI. Open `pipeline_code_5.0.R` in RStudio or another R session and run it section by section. There is no build step, test suite, or lint command.

**Section 1 (CONFIGURATION)** is the only section users need to edit for each run:

```r
RUN_NAME   # label for this sequencing run (used in output filenames)
PRIMERNAME # must match a name in primer.data.csv (e.g., "MV1", "MFU", "DL")
BLASTSCOPE # "vertebrate" or "eukaryote"
SP_THOLD   # % identity floor for species-level assignment

PARENT_DIR    # parent directory containing a "Fastq/" subfolder with raw reads
PROCESSED_DIR # final run folder is created here as <RUN_NAME>_<PRIMERNAME>/
DB_DIR        # directory holding per-primer <PRIMERNAME>_database.csv files

CUTADAPT      # "cutadapt" if on PATH, otherwise the full path to the binary
TAXONKIT_PATH # full path to the taxonkit binary

SKIP_DADA2    # TRUE to reload a previous ASV table and re-run only annotation
```

FASTQs must live in `PARENT_DIR/Fastq/` and be named `<PRIMERNAME>_<sample>_R1_*.fastq.gz` (hyphen separator also accepted). Sample names are extracted by stripping the primer prefix and the `_R1`/`_R2` suffix, so they may contain any number of underscores.

## Architecture

### Pipeline stages (in order)

1. **Primer trimming** — `make_primer_shell_script.R` is sourced to generate `trim_primers.sh`, which runs `cutadapt` on all FASTQs and writes trimmed reads to `for_dada2/`.

2. **DADA2 denoising** (`pipeline_code_*.R`) — quality filtering → error learning → denoising (`dada()`) → paired-end merging → chimera removal → size filtering by `MAX_AMPLICON_LENGTH`/`MIN_AMPLICON_LENGTH` from `primer.data.csv`.

3. **Hashing** — every unique ASV sequence is assigned a SHA1 hash (via `digest`). Outputs: `<RUN>_<PRIMER>_hash_key.csv` and `_ASV_table.csv`.

4. **BLAST on CEG cluster** (`CEG_BLAST_function_2.0.R`, `ceg_blast()`) — SSHes into `frustule.ocean.washington.edu:3004` using `~/.ssh/ceg_rsa`, uploads the FASTA of new sequences, runs blastn against `nt_euk` or `nt_vrt`, and downloads results. The shell templates `MURIblast_eukaryote_template*.sh` / `MURIblast_vertebrate_template*.sh` are the BLAST job scripts modified in-place before upload.

5. **LCA assignment** (`LCA_function_4.0.R`, `LCA()`) — applies a tiered percent-identity decision tree (100% → >99.3% → >98% → >96% → >94% → all) to select BLAST hits per ASV, then calls `taxonkit lca` + `taxonkit reformat` (via temp files in `tmp/`) to resolve Lowest Common Ancestor lineages. Species-level assignments are blanked when `Max_pident < SP_THOLD`. Each taxon's ASVs receive sequential `HaplotypeNumber`s. Results are written to/appended into a per-primer hash database CSV at `DATABASE_LOCATION`.

6. **Output tables** — three CSVs per run: `taxon_table.csv` (long), `taxon_table_wide.csv` (pivoted), `haplotype_table.csv` (one row per taxon+haplotype).

### Key files

| File | Purpose |
|------|---------|
| `pipeline_code_5.0.R` | Current main pipeline |
| `CEG_BLAST_function_2.0.R` | `ceg_blast()` — SSH upload → BLAST → download |
| `LCA_function_4.0.R` | `LCA()` — pident filtering, taxonkit LCA, haplotype assignment |
| `primer.data.csv` | Per-primer sequences, lengths, amplicon size bounds, overlap |
| `MURI_taxids_to_exclude.txt` | NCBI taxids excluded from BLAST (marine mammal species) |
| `MURIblast_eukaryote_template_2.0.sh` | CEG BLAST job template for eukaryote scope |
| `MURIblast_vertebrate_template_2.0.sh` | CEG BLAST job template for vertebrate scope |

### Directory layout

```
PARENT_DIR/                            # raw data root (temporary contents cleaned up)
  Fastq/                               # raw paired-end FASTQs (input, untouched)
  for_dada2/                           # primer-trimmed reads (removed on CLEANUP)
  filtered/                            # quality-filtered reads (removed on CLEANUP)

PROCESSED_DIR/<RUN_NAME>_<PRIMERNAME>/ # permanent run archive
  code_etc/                            # copies of scripts + primer data for this run
  outputs/                             # ASV table, hash key FASTA, taxon tables, BLAST results
  logs/                                # cutadapt report, read-tracking CSV, filtered-ASV CSV
```

Outputs and logs are written directly to the run archive — there is no intermediate move step. Cleanup only removes `for_dada2/` and `filtered/` under `PARENT_DIR`.

### Hash database

`DB_DIR/<PRIMERNAME>_database.csv` is a cumulative per-primer annotation cache. `LCA()` skips re-annotating hashes already in this database and appends new results. Set `SKIP_DADA2 = TRUE` to reload a previous ASV table and re-run only the annotation stage.

## External dependencies

- **R packages**: tidyverse, dada2, digest, seqinr, ssh, sys, ShortRead, here
- **System tools**: `cutadapt`, `taxonkit` (paths set in the script)
- **Remote cluster**: CEG server at `frustule.ocean.washington.edu:3004`, SSH key `~/.ssh/ceg_rsa`, login `KellyCEG`. BLAST databases live at `/mnt/nfs/home/KellyCEG/blastdb_euk/` (eukaryote) on the server.
- **taxonkit** uses temp files in `tmp/` at the repo root (`taxids.txt`, `taxids_lca.txt`, `taxids_lineages.txt`).
